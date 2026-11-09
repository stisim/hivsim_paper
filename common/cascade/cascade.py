"""
A condensed HIV care cascade as Starsim deliveries: testing pathways, staging and linkage, pre-ART, ART, dropout and
return, VMMC, PrEP and PMTCT, driven by a table of parameters (params.py) and acting on the model through a
`CascadeHost` (interface.py). Needs Starsim 4.0.0rc3 or later (cycles through delays).

Each person has one place in care (`people.care`): none, uptake (offered routine testing now and then), testing (in
a routine testing loop), staging, pre-ART, linking (waiting to start ART), on ART, out (dropped out of ART and not
back in care) or lost (for good). Testing pathways reach those not in care (none, uptake or testing); staging those
not further on. The pathways:

- hct: at debut (and, at its start year, everyone past debut), into the testing loop with probability `p_debut`,
  else into uptake; from uptake, after `uptake_interval`, into the loop with probability `p_uptake`, else again; in the
  loop, a test every `test_interval`. A positive test starts ART at once with probability `p_immediate_art` (by sex),
  else goes to staging.
- symptomatic: at first symptoms, a test with probability `p_present` (positive: staging), else straight to the
  assessment of staging with probability `p_present_again`.
- anc: at 12 weeks of pregnancy, a test with probability `p_test`; a positive one leads to staging with probability
  `p_stage`, and to PMTCT (see pmtct).
- infant: at 6 weeks of age, a test with probability `p_test` (positive: staging).
  A further testing pathway is a row with a trigger (see `TESTING`) and its probability.
- staging: lost with probability `p_lost` (back to uptake); otherwise, after `delay`, a CD4 count and the
  eligibility rule (pathway `eligibility`: WHO stage, pregnancy, children's ages, CD4 thresholds, by year); eligible:
  linking with probability `p_link_art` (else uptake); not: pre-ART with probability `p_link_preart` (else uptake).
- preart: every `interval`, lost with probability `p_lost` (back to uptake), else eligibility again.
- art: a start `start_delay` after linking (or from pre-ART), in the suppressive arm with probability
  `p_suppressive` (by sex, under and over `age_split`), else the non-suppressive one; dropout after `dropout_delay`.
- dropout: lost for good with probability `p_lost`; else out of care with probability `p_out` (else back to uptake),
  and those out start ART again with probability `p_restart`, and are offered dropout and return again after
  `dropout_delay`. With `p_restart` 0 (as in EMOD Malawi from 2009), those out stay out of care, apart from `p_lost`.
- pmtct: after a positive ANC test, combination PMTCT with probability `p_combination` (Option B with probability
  `p_option_b`, else A) or single-dose nevirapine; with `blocked_if_staged` 1, women sent to staging then get none of
  the combination (as EMOD's campaign does, through its update order).
- vmmc, prep: numbers a year by group (the flows table), chosen at random among those eligible.
"""
import numpy as np
import pandas as pd
import starsim as ss
from starsim.internals.conditions import Cond

NONE, UPTAKE, TESTING, STAGING, PREART, LINKING, ON_ART, OUT, LOST = range(9)
CARE = ['none', 'uptake', 'testing', 'staging', 'preart', 'linking', 'on_art', 'out', 'lost']
ARMS = ['', 'suppressive', 'nonsuppressive'] # people.arm: the arm of the last ART start
TESTING_FROM = [NONE, UPTAKE, TESTING]           # Reached by testing pathways
STAGING_FROM = [NONE, UPTAKE, TESTING, STAGING]  # Reached by staging
RISKS = ['LOW', 'MEDIUM', 'HIGH'] # Risk groups a flow can be restricted to (see CascadeHost.risk())
TRIGGERS = dict(debut='debuts', symptoms='newly_symptomatic', anc='twelve_weeks_pregnant', infant='six_weeks_old') # Events of the host


def ev(name, outcome=None):
    return Cond('event', name, outcome)


def any_of(conds):
    out = None
    for c in conds:
        out = c if out is None else out | c
    return out


def care_in(codes):
    return any_of([ss.person.care == k for k in codes])


class Holder(ss.Intervention):
    """ The cascade's people (`care`), its host, parameters and flows """
    def __init__(self, params, host, name='cascade'):
        super().__init__(name=name)
        self.params, self.host_class = params, host
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.care = sim.people.add_state(ss.IntArr('care', default=NONE))
        self.cx_dying = sim.people.add_state(ss.BoolArr('cx_dying', default=False))
        self.cx_postdebut = sim.people.add_state(ss.BoolArr('cx_postdebut', default=False))
        self.vmmc = sim.people.add_state(ss.BoolArr('cx_circ', default=False)) # Circumcised by the programme
        self.arm = sim.people.add_state(ss.IntArr('arm', default=0))
        self.host = self.host_class(sim)
        self.rng = np.random.default_rng(sim.pars.rand_seed + 7)
        self.flow_plans = {}
        return

    def step(self):
        au = self.sim.people.auids
        self.cx_dying[au] = self.host.dying(au)
        self.cx_postdebut[au] = self.host.post_debut(au)
        self.year = self.t.now('year')
        return

    def p(self, pathway, parameter, sex='all'):
        return self.params.p(pathway, parameter, self.year, sex)

    def flow(self, intervention):
        """ Those chosen this step for a flow (VMMC, PrEP): each group's number for its period, spread evenly over the period's steps (remainder first), at random among the HIV-negative eligible """
        f, host, ppl = self.params.flows, self.host, self.sim.people
        out = []
        for r in f[f.intervention == intervention].itertuples():
            if not (r.start <= self.year < r.end):
                continue
            scale = self.params.value('flows', 'scale') if self.params.has('flows', 'scale') else 1/self.sim.pop_scale # Agents per person in the table
            total = int(r.n*scale)
            steps = int(round((r.end - r.start)/self.t.dt_year))
            k = int(round((self.year - r.start)/self.t.dt_year))
            n = total//steps + (k < total % steps)
            if not n:
                continue
            u = ppl.alive.uids
            u = u[~self.cx_dying.raw[u] & (host.female(u) == (r.sex == 'f'))]
            if r.node != 'all':
                u = u[host.node(u) == int(r.node)]
            age = host.age_years(u)
            u = u[(age >= r.age_lo) & (age < r.age_hi) & ~host.infected(u)]
            for res in str(r.restrict).split(';') if isinstance(r.restrict, str) else []:
                key, val = res.split(':')
                if key == 'risk':
                    u = u[np.isin(host.risk(u), [RISKS.index(v) for v in val.split('|')])]
                elif key == 'not' and val == 'vmmc':
                    u = u[~self.vmmc.raw[u]]
            out.append(u[self.rng.permutation(len(u))[:n]])
        return ss.uids(np.concatenate(out)) if out else ss.uids()


class CareProduct(ss.Product):
    """ A product acting on the cascade's people through the holder """
    def __init__(self, name, **kw):
        super().__init__(name=name)
        self.kw = kw

    @property
    def cx(self):
        return self.sim.interventions.cascade


class SetCare(CareProduct):
    def administer(self, uids):
        self.cx.care[uids] = self.kw['to']


class Test(CareProduct):
    """ An HIV test, through the host: outcomes positive and negative """
    def outcome_names(self):
        return ['positive', 'negative']

    def administer(self, uids):
        pos = self.cx.host.test(uids)
        return dict(positive=uids[pos], negative=uids[~pos])


class Assess(CareProduct):
    """ Staging: a CD4 count, then the eligibility rule of the year """
    def outcome_names(self):
        return ['treat', 'not_yet'] # Not 'eligible': a result every delivery has

    def administer(self, uids):
        cx = self.cx
        P, year, h = cx.params, cx.year, cx.host
        cd4 = h.record_cd4(uids)
        who = h.who_stage(uids)
        age = h.age_years(uids)
        preg = h.pregnant(uids)
        v = lambda name, default: P.value('eligibility', name, year) if P.has('eligibility', name) else default
        adult = age >= v('adult_age', 5)
        ok_adult = (who >= v('who_adult', 6)) | (preg & (v('pregnant', 0) > 0))
        ok_child = (age*365 <= v('child_age', -1)*365) | (who >= v('who_child', 6))
        ok_cd4 = (cd4 <= v('cd4', -np.inf)) | (preg & (cd4 <= v('cd4_pregnant', -np.inf)))
        ok = h.infected(uids) & (np.where(adult, ok_adult, ok_child) | ok_cd4)
        return dict(treat=uids[ok], not_yet=uids[~ok])


class ANC(CareProduct):
    """ After a positive ANC test: staging, and PMTCT (combination, blocked for those staged if `blocked_if_staged`) """
    def outcome_names(self):
        return ['staged']

    def administer(self, uids):
        cx = self.cx
        r = cx.rng.random((3, len(uids)))
        staged = r[0] < cx.p('anc', 'p_stage')
        comb = r[1] < cx.p('pmtct', 'p_combination')
        blocked = staged & (cx.params.value('pmtct', 'blocked_if_staged') > 0)
        b = r[2] < cx.p('pmtct', 'p_option_b')
        days = cx.params.value('pmtct', 'days')
        for mask, eff in ((comb & ~blocked & b, 'eff_option_b'), (comb & ~blocked & ~b, 'eff_option_a'), (~comb, 'eff_sdnvp')):
            if mask.any():
                cx.host.start_pmtct(uids[mask], cx.params.value('pmtct', eff), days)
        return dict(staged=uids[staged])


class StartART(CareProduct):
    """ ART, in the suppressive arm with probability by sex and age """
    def administer(self, uids):
        cx = self.cx
        h = cx.host
        female, age = h.female(uids), h.age_years(uids)
        split = cx.params.value('art', 'age_split')
        p = np.array([cx.p('art', 'p_suppressive_under' if a < split else 'p_suppressive_over', 'f' if f else 'm') for f, a in zip(female, age)])
        supp = cx.rng.random(len(uids)) < p
        for mask, arm in ((supp, 'suppressive'), (~supp, 'nonsuppressive')):
            if mask.any():
                h.start_art(uids[mask], arm)
                cx.arm[uids[mask]] = ARMS.index(arm)
        cx.care[uids] = ON_ART


class StopART(CareProduct):
    def administer(self, uids):
        self.cx.host.stop_art(uids[self.cx.host.on_art(uids)])


class Circumcise(CareProduct):
    def administer(self, uids):
        done = self.cx.host.circumcise(uids)
        self.cx.vmmc[done] = True


class PrEP(CareProduct):
    def administer(self, uids):
        P = self.cx.params
        self.cx.host.start_prep(uids, P.value('prep', 'efficacy'), P.value('prep', 'days'))


def each_step():
    """ Rounds at every step (the models' steps are months), for a delivery to those a function of the sim chooses """
    return ss.schedule(every=ss.months(1))


def delay(params, pathway, parameter):
    form, v = params.delay(pathway, parameter)
    return dict(exp_days=lambda: ss.expon(ss.days(v[0])), uniform_days=lambda: ss.uniform(ss.days(v[0]), ss.days(v[1])), days=lambda: ss.days(v[0]))[form]()


def condensed_cascade(params, host, testing=None):
    """
    The holder and the deliveries of the cascade

    Args:
        params (Params): the parameter table
        host (class): the model's `CascadeHost`
        testing (list): further testing pathways, as (name, trigger) with trigger a key of `TRIGGERS`, each with
            `p_test` (and optionally `start`) in the table; a positive test leads to staging
    """
    P = params
    yearfn = lambda pathway, parameter, sex='all': (lambda sim: sim.interventions.cascade.p(pathway, parameter, sex))
    start = lambda pathway: dict(start=P.value(pathway, 'start')) if P.has(pathway, 'start') else {}
    ok = ~ss.person.cx_dying
    testable, stageable = care_in(TESTING_FROM) & ok, care_in(STAGING_FROM) & ok
    D = []
    def add(name, when, who, what, delay=None, **kw):
        when = ss.on(when, delay=delay) if delay is not None else when
        D.append(ss.Delivery(what, name=name, when=when, who=who, disease=kw.pop('disease', None) or host_disease, **kw))
        return name
    host_disease = P.value('host', 'disease') if P.has('host', 'disease') else 'hiv'

    # Events of the host
    for key, method in TRIGGERS.items():
        add(f'cx_ev_{key}', each_step(), lambda sim, m=method: ss.uids(getattr(sim.interventions.cascade.host, m)()), None)
    to = lambda name, code: SetCare(name + '_p', to=code)

    # Routine testing (hct): uptake and the testing loop
    D.append(ss.Delivery(to('cx_hct_init', UPTAKE), name='cx_hct_init', disease=host_disease, when=[P.value('hct', 'start')],
                         who=(ss.person.care == NONE) & ok & ss.person.cx_postdebut))
    add('cx_hct_debut_yes', ev('cx_ev_debut'), testable, to('cx_hct_debut_yes', TESTING), coverage=yearfn('hct', 'p_debut'), **start('hct'))
    add('cx_hct_debut_no', ev('cx_ev_debut') & ~ev('cx_hct_debut_yes'), testable, to('cx_hct_debut_no', UPTAKE), **start('hct'))
    into_uptake = ['cx_hct_init', 'cx_hct_debut_no', 'cx_uptake_no', 'cx_stage_lost', 'cx_link_no', 'cx_prelink_no', 'cx_preart_lost', 'cx_drop_back']
    into_testing = ['cx_hct_debut_yes', 'cx_uptake_yes', 'cx_hct_retest']
    add('cx_uptake_wait', any_of([ev(n) for n in into_uptake]), (ss.person.care == UPTAKE) & ok, None, delay=delay(P, 'hct', 'uptake_interval'))
    add('cx_uptake_yes', ev('cx_uptake_wait'), ss.person.care == UPTAKE, to('cx_uptake_yes', TESTING), coverage=yearfn('hct', 'p_uptake'))
    add('cx_uptake_no', ev('cx_uptake_wait') & ~ev('cx_uptake_yes'), ss.person.care == UPTAKE, to('cx_uptake_no', UPTAKE))
    add('cx_test_wait', any_of([ev(n) for n in into_testing]), (ss.person.care == TESTING) & ok, None, delay=delay(P, 'hct', 'test_interval'))
    add('cx_hct_test', ev('cx_test_wait'), ss.person.care == TESTING, Test('cx_hct_test_p'))
    add('cx_hct_retest', ev('cx_hct_test', 'negative'), None, to('cx_hct_retest', TESTING))
    for sex, cond in (('f', ss.person.female), ('m', ~ss.person.female)):
        add(f'cx_immediate_{sex}', ev('cx_hct_test', 'positive'), cond, StartART(f'cx_immediate_{sex}_p'), coverage=yearfn('hct', 'p_immediate_art', sex))
    into_staging = [ev('cx_hct_test', 'positive') & ~ev('cx_immediate_f') & ~ev('cx_immediate_m')]

    # Testing triggered by events: symptoms, ANC, infants, and any others given
    add('cx_sym_test', ev('cx_ev_symptoms'), testable, Test('cx_sym_test_p'), coverage=yearfn('symptomatic', 'p_present'), **start('symptomatic'))
    into_staging.append(ev('cx_sym_test', 'positive'))
    add('cx_sym_again', ev('cx_ev_symptoms') & ~ev('cx_sym_test'), testable, to('cx_sym_again', STAGING), coverage=yearfn('symptomatic', 'p_present_again'), **start('symptomatic'))
    add('cx_anc_test', ev('cx_ev_anc'), testable, Test('cx_anc_test_p'), coverage=yearfn('anc', 'p_test'), **start('anc'))
    add('cx_anc', ev('cx_anc_test', 'positive'), None, ANC('cx_anc_p'))
    into_staging.append(ev('cx_anc', 'staged'))
    for name, trigger in [('infant', 'infant')] + list(testing or []):
        add(f'cx_{name}_test', ev(f'cx_ev_{trigger}'), testable, Test(f'cx_{name}_test_p'), coverage=yearfn(name, 'p_test'), **start(name))
        into_staging.append(ev(f'cx_{name}_test', 'positive'))

    # Staging, linkage, pre-ART
    staging = any_of(into_staging)
    add('cx_stage_lost', staging, stageable, to('cx_stage_lost', UPTAKE), coverage=yearfn('staging', 'p_lost'))
    add('cx_stage_in', staging & ~ev('cx_stage_lost'), stageable, to('cx_stage_in', STAGING))
    add('cx_assess', ev('cx_stage_in') | ev('cx_sym_again'), (ss.person.care == STAGING) & ok, Assess('cx_assess_p'), delay=delay(P, 'staging', 'delay'))
    for sex, cond in (('f', ss.person.female), ('m', ~ss.person.female)):
        add(f'cx_link_{sex}', ev('cx_assess', 'treat'), cond, to(f'cx_link_{sex}', LINKING), coverage=yearfn('staging', 'p_link_art', sex))
    add('cx_link_no', ev('cx_assess', 'treat') & ~ev('cx_link_f') & ~ev('cx_link_m'), None, to('cx_link_no', UPTAKE))
    add('cx_prelink', ev('cx_assess', 'not_yet'), None, to('cx_prelink', PREART), coverage=yearfn('staging', 'p_link_preart'))
    add('cx_prelink_no', ev('cx_assess', 'not_yet') & ~ev('cx_prelink'), None, to('cx_prelink_no', UPTAKE))
    add('cx_preart_wait', ev('cx_prelink') | ev('cx_preart_again'), (ss.person.care == PREART) & ok, None, delay=delay(P, 'preart', 'interval'))
    add('cx_preart_lost', ev('cx_preart_wait'), ss.person.care == PREART, to('cx_preart_lost', UPTAKE), coverage=yearfn('preart', 'p_lost'))
    add('cx_preart_check', ev('cx_preart_wait') & ~ev('cx_preart_lost'), ss.person.care == PREART, Assess('cx_preart_check_p'))
    add('cx_preart_again', ev('cx_preart_check', 'not_yet'), None, to('cx_preart_again', PREART))

    # ART, dropout and return
    add('cx_art_wait', ev('cx_link_f') | ev('cx_link_m') | ev('cx_preart_check', 'treat'), care_in([LINKING, PREART]) & ok, None, delay=delay(P, 'art', 'start_delay'))
    add('cx_art_start', ev('cx_art_wait'), care_in([LINKING, PREART]) & ok, StartART('cx_art_start_p'))
    add('cx_drop_wait', ev('cx_art_start') | ev('cx_immediate_f') | ev('cx_immediate_m') | ev('cx_drop_out') | ev('cx_restart'), ok, None, delay=delay(P, 'art', 'dropout_delay'))
    add('cx_drop', ev('cx_drop_wait'), care_in([ON_ART, OUT]) & ok, StopART('cx_drop_p'))
    add('cx_drop_lost', ev('cx_drop'), None, to('cx_drop_lost', LOST), coverage=yearfn('dropout', 'p_lost'))
    add('cx_drop_stay', ev('cx_drop') & ~ev('cx_drop_lost'), None, None, coverage=yearfn('dropout', 'p_out'))
    add('cx_restart', ev('cx_drop_stay'), None, StartART('cx_restart_p'), coverage=yearfn('dropout', 'p_restart'))
    add('cx_drop_out', ev('cx_drop_stay') & ~ev('cx_restart'), None, to('cx_drop_out', OUT))
    add('cx_drop_back', ev('cx_drop') & ~ev('cx_drop_lost') & ~ev('cx_drop_stay'), None, to('cx_drop_back', UPTAKE))

    # Prevention: flows of VMMC and PrEP
    add('cx_vmmc', each_step(), lambda sim: sim.interventions.cascade.flow('vmmc'), Circumcise('cx_vmmc_p'))
    add('cx_prep', each_step(), lambda sim: sim.interventions.cascade.flow('prep'), PrEP('cx_prep_p'))
    return [Holder(P, host)] + D
