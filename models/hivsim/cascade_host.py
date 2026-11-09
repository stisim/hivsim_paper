"""
HIVsim's side of the condensed cascade (common/cascade): `StisimHost` implements its `CascadeHost` on STIsim's
`sti.HIV` and `sti.ART`, the pregnancy module and the sexual network. Needs STIsim with HIV.who_stage(),
new_symptoms(), start_pmtct() and ART.start_art(effective=) (branch rs-cascade-host).

`sti.ART` must be in the sim, as the holder of ART's states (CD4 on ART, mortality, infectiousness), with no coverage
target and no stopping of its own (`art_for_cascade()`), so the cascade alone starts and stops ART.
"""
import sys
from pathlib import Path
import numpy as np
import starsim as ss
import stisim as sti

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.cascade.interface import CascadeHost # noqa: E402


def art_for_cascade(**kwargs):
    """ sti.ART holding the states of ART for the cascade: it starts no one, and stops no one of its own accord """
    return sti.ART(dur_on_art=ss.constant(ss.years(200)), **kwargs)


class StisimHost(CascadeHost):
    def __init__(self, sim):
        super().__init__(sim)
        self.hiv = sim.diseases.hiv
        self.art = sim.interventions.art
        self.ppl = sim.people
        self.preg = sim.demographics.pregnancy
        self.net = sim.networks.structuredsexual
        self.dt_years = float(sim.t.dt_year)
        return

    # People
    def age_years(self, uids):
        return self.ppl.age.raw[uids]

    def female(self, uids):
        return self.ppl.female.raw[uids]

    def pregnant(self, uids):
        return self.preg.P.raw[uids]

    def post_debut(self, uids):
        return self.ppl.age.raw[uids] >= self.net.debut.raw[uids]

    def dying(self, uids):
        return self.ppl.ti_dead.raw[uids] <= self.sim.ti

    def risk(self, uids):
        return self.net.risk_group.raw[uids]

    # Events
    def crossing(self, age):
        """ Those who reach `age` (a number, or one for each person alive) in this step """
        au = self.ppl.auids
        a = self.ppl.age.raw[au]
        return au[(a >= age) & (a - self.dt_years < age)]

    def debuts(self):
        return self.crossing(self.net.debut.raw[self.ppl.auids])

    def newly_symptomatic(self):
        return self.hiv.new_symptoms()

    def twelve_weeks_pregnant(self):
        u = self.preg.P.uids
        weeks = (self.sim.ti - self.preg.ti_P.raw[u])*self.dt_years*52.18 # Since conception
        return u[(weeks >= 12) & (weeks - self.dt_years*52.18 < 12)]

    def six_weeks_old(self):
        return self.crossing(6/52.18)

    # HIV
    def infected(self, uids):
        return self.hiv.infected.raw[uids]

    def test(self, uids):
        pos = uids[self.hiv.infected.raw[uids]]
        new = pos[~self.hiv.diagnosed.raw[pos]]
        self.hiv.diagnosed[new] = True
        self.hiv.ti_diagnosed[new] = self.sim.ti
        return self.hiv.infected.raw[uids]

    def cd4(self, uids):
        return np.where(self.hiv.infected.raw[uids], self.hiv.cd4.raw[uids], np.inf)

    def record_cd4(self, uids):
        """ The lowest CD4 count each has had (STIsim's nadir, the lowest off ART, rather than the lowest measured) """
        return np.where(self.hiv.infected.raw[uids], np.minimum(self.hiv.cd4_nadir.raw[uids], self.hiv.cd4.raw[uids]), np.inf)

    def who_stage(self, uids):
        return np.maximum(self.hiv.who_stage(uids), 1)

    def on_art(self, uids):
        return self.art.on_treatment.raw[uids]

    def start_art(self, uids, arm):
        uids = uids[self.hiv.infected.raw[uids] & ~self.art.on_treatment.raw[uids]]
        if len(uids):
            self.art.start_art(uids, effective=np.full(len(uids), arm == 'suppressive'))

    def stop_art(self, uids):
        uids = uids[self.art.on_treatment.raw[uids]]
        if len(uids):
            self.art.stop_art(uids)

    def start_pmtct(self, uids, efficacy, days):
        self.hiv.start_pmtct(uids, eff=efficacy, dur=ss.days(days))

    def circumcise(self, uids):
        men = uids[~self.ppl.female.raw[uids] & ~self.hiv.circumcised.raw[uids]]
        self.hiv.circumcise(men)
        return men

    def start_prep(self, uids, efficacy, days):
        self.hiv.start_prep(uids[~self.hiv.infected.raw[uids]], eff=efficacy, dur=ss.days(days))

    def cascade_stage(self, uids):
        return dict(diagnosed=self.hiv.diagnosed.raw[uids], on_art=self.on_art(uids), suppressive=self.hiv.on_effective_art.raw[uids])
