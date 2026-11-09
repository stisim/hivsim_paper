"""
EMOD′'s side of the condensed cascade (common/cascade): `EmodHost` implements its `CascadeHost` on EmodHIV, the
demography and the risk module. ART arms are the campaign's ARTMortalityTable interventions (by arm and sex).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from common.cascade.interface import CascadeHost # noqa: E402

ARMS = dict(suppressive='EffectiveART', nonsuppressive='NonSuppressiveART') # The campaign's events that start each arm


def art_tables(events):
    """ {(arm, sex): ARTMortalityTable config} from the campaign """
    out = {}
    for ev in events:
        ic = ev['Event_Coordinator_Config'].get('Intervention_Config', {})
        act = ic.get('Actual_IndividualIntervention_Config', {})
        if act.get('class') == 'ARTMortalityTable':
            arm = next(a for a, e in ARMS.items() if e in ic['Trigger_Condition_List'])
            out[(arm, ic['Target_Gender'].lower()[0])] = act
    return out


class EmodHost(CascadeHost):
    def __init__(self, sim):
        super().__init__(sim)
        self.hiv = sim.diseases.hiv
        self.ppl = sim.people
        self.tables = art_tables(self.hiv.inp.events)
        self.prev_twelve = np.zeros(0, np.int64)
        return

    # People
    def age_years(self, uids):
        return self.ppl.emod_age.raw[uids]/365

    def female(self, uids):
        return self.ppl.female.raw[uids]

    def pregnant(self, uids):
        return ~np.isnan(self.sim.demographics.emoddemog.preg_days.raw[uids])

    def post_debut(self, uids):
        return self.ppl.emod_age.raw[uids]/365 >= self.sim.networks.emodnet.debut.raw[uids]

    def dying(self, uids):
        return np.isin(uids, self.hiv.dead_uids)

    def node(self, uids):
        return self.ppl.node.raw[uids] + 1

    def risk(self, uids):
        return self.ppl.risk.raw[uids]

    # Events
    def debuts(self):
        return self.sim.interventions.emodrisk.new_debuts

    def newly_symptomatic(self):
        return self.hiv.newly_symptomatic

    def twelve_weeks_pregnant(self):
        """ Those of the step before: EMOD's event comes after its interventions' update """
        out, self.prev_twelve = self.prev_twelve, self.sim.demographics.emoddemog.twelve_weeks
        return out

    def six_weeks_old(self):
        au, dt = self.ppl.auids, self.hiv.dt_days
        age = self.ppl.emod_age.raw[au]
        return au[(age - 2*dt < 42) & (age - dt >= 42)] # Six weeks old on the step before, likewise

    # HIV
    def infected(self, uids):
        return self.hiv.infected.raw[uids]

    def test(self, uids):
        h = self.hiv
        pos = h.infected.raw[uids].copy()
        h.ever_tested[uids] = True
        h.last_test[uids] = h.t.now('year')
        h.diagnosed[uids[pos]] = True
        return pos

    def cd4(self, uids):
        return np.where(self.hiv.infected.raw[uids], self.hiv.cd4(uids), np.inf)

    def record_cd4(self, uids):
        h = self.hiv
        h.lowest_cd4[uids] = np.minimum(h.lowest_cd4.raw[uids], self.cd4(uids))
        return h.lowest_cd4.raw[uids]

    def who_stage(self, uids):
        return np.where(self.hiv.infected.raw[uids], self.hiv.who_stage(uids), 1.0)

    def on_art(self, uids):
        return self.hiv.on_art(uids)

    def start_art(self, uids, arm):
        for sex, mask in (('f', self.ppl.female.raw[uids]), ('m', ~self.ppl.female.raw[uids])):
            t = self.tables[(arm, sex)]
            self.hiv.go_on_art(uids[mask], t, t['ART_Multiplier_On_Transmission_Prob_Per_Act'], t['Days_To_Achieve_Viral_Suppression'])
        return

    def stop_art(self, uids):
        self.hiv.go_off_art(uids)

    def start_pmtct(self, uids, efficacy, days):
        h, now = self.hiv, self.hiv.t.now('year')
        h.pmtct[uids] = np.maximum(np.where(h.pmtct_until.raw[uids] >= now, h.pmtct.raw[uids], 0), efficacy) # The larger, as EMOD's
        h.pmtct_until[uids] = now + days/365

    def circumcise(self, uids):
        men = uids[~self.ppl.female.raw[uids] & ~self.ppl.circumcised.raw[uids]]
        self.ppl.circumcised[men] = True
        return men

    def start_prep(self, uids, efficacy, days):
        self.hiv.prep_eff[uids] = efficacy
        self.hiv.prep_until[uids] = self.hiv.t.now('year') + days/365

    def cascade_stage(self, uids):
        return dict(diagnosed=self.hiv.diagnosed.raw[uids], on_art=self.hiv.on_art(uids), suppressive=self.ppl.arm.raw[uids] == 1) # people.arm: the cascade's arm of the last start
