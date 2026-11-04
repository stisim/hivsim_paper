"""
EMOD′ HIV: transmission per coital act and natural history without treatment, as in EMOD-HIV (docs/emod_inventory.md,
sections 4, 5, 10 (traditional circumcision), 11 (transmission at birth) and 12). In each step, in EMOD's order:

1. Campaign: traditional circumcision from 1975 (everyone then, and each newborn, with a probability by node; EMOD
   draws for both sexes and circumcises only men), and the seeding at SeedYr: a share of everyone at Risk HIGH is infected, partway through their
   course (a uniform fraction of their prognosis, InfectionHIV::SetParameters).
2. Transmission at birth: a child of a mother infected at the birth is infected with probability 0.3.
3. Coital acts (Relationship::Consummate, from the network's `acts`): each act in a relationship with one partner infected
   transmits with probability
   Base_Infectivity × stage (26 acute, 10 AIDS, else 1) × max of the STI multipliers of the two (5.5 if the infected
   partner, or the uninfected, has the STI flag) × (1 − 0.9) if a condom is used (drawn for each act, from the edge's
   `condoms`) × (1 − 0.6) for a circumcised man × the male-to-female multiplier for a woman (by her age, linear between
   0, 15 and 25 years). Stages are as at the start of the step. Each uninfected person's acts are taken in random order
   and the first that transmits infects them (RelationshipGroups::ExposeToContagion).
4. Progression of those infected before this step's exposure: the time since infection grows by the step; acute for
   2.9 months, AIDS for the last 9 months of the prognosis, death from HIV at the end of it (InfectionHIV::SetStageFromDuration).
5. The new infections of step 3. The prognosis is drawn at infection from the age then (FerrandAgeDependentDistribution):
   under 15, an exponential (rate 1.52/year) with probability 0.57, else a Weibull (median 16 years, shape 2.7); from 15,
   a Weibull of shape 2 and scale 21.182 − 0.2717 × age (age capped at 50); at least 7 days.

CD4 is deterministic (both heterogeneities 0): its square root falls linearly from √540.55 at infection to √31.63 at
death (SusceptibilityHIV). People: `people.sti` (from the risk module), `people.circumcised`.
"""
import numpy as np
import pandas as pd
import starsim as ss

ACUTE, LATENT, AIDS = 0, 1, 2


class EmodHIV(ss.Disease):
    """
    Args:
        inp (objdict): EMOD's inputs for one parameter set, from `inputs.load()`
        net (str): the name of the network
        record (bool): whether to keep a record of each infection and HIV death (`infections`, `deaths`)
    """
    def __init__(self, inp, net='emodnet', record=True, name='hiv', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp, self.netname, self.record = inp, net, record
        self.define_states(
            ss.BoolState('infected', label='Infected'),
            ss.FloatArr('days'),       # Days since infection
            ss.FloatArr('prognosis'),  # Days from infection to death from HIV
            ss.FloatArr('acute_days'),
            ss.FloatArr('aids_days'),
            ss.IntArr('stage', default=-1),
            ss.FloatArr('year_infected'),
            ss.IntArr('n_infections', default=0), # How many this person has infected
        )
        self.infections, self.deaths = [], []
        self.seeded = self.circ_started = False
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.sti = sim.people.add_state(ss.BoolArr('sti', default=False))
        self.circumcised = sim.people.add_state(ss.BoolArr('circumcised', default=False))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 5)
        self.dt_days = self.inp.dt
        h = self.inp.hiv
        self.acute_len, self.aids_len = h.Acute_Duration_In_Months*365/12, h.AIDS_Duration_In_Months*365/12
        return

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('new_infections', dtype=int, scale=True, label='New infections'),
                            ss.Result('hiv_deaths', dtype=int, scale=True, label='HIV deaths'),
                            ss.Result('n_infected', dtype=int, scale=True, label='People with HIV'))
        return

    # %% Natural history

    def draw_prognosis(self, age):
        """ Days from infection to death from HIV, by age at infection (years) """
        h, u = self.inp.hiv, self.rng.random(len(age))
        child = age < h.HIV_Age_Max_for_Child_Survival_Function
        a = h.HIV_Child_Survival_Rapid_Progressor_Fraction
        rapid = child & (u < a)
        slow = child & ~rapid
        years = np.zeros(len(age))
        years[rapid] = -np.log(u[rapid]/a)/h.HIV_Child_Survival_Rapid_Progressor_Rate
        years[slow] = h.HIV_Child_Survival_Slow_Progressor_Scale*(-np.log((u[slow] - a)/(1 - a))/np.log(2))**(1/h.HIV_Child_Survival_Slow_Progressor_Shape)
        adult = np.minimum(age[~child], h.HIV_Age_Max_for_Adult_Age_Dependent_Survival)
        lam = h.HIV_Adult_Survival_Scale_Parameter_Slope*adult + h.HIV_Adult_Survival_Scale_Parameter_Intercept
        years[~child] = lam*(-np.log(1 - u[~child]))**(1/h.HIV_Adult_Survival_Shape_Parameter)
        return np.maximum(years*365, 7)

    def infect(self, uids, fast_forward=False):
        """ Start infections; seeded ones partway through their course """
        uids = uids[~self.infected.raw[uids]]
        if not len(uids):
            return uids
        ppl = self.sim.people
        prog = self.draw_prognosis(ppl.emod_age.raw[uids]/365)
        acute = np.minimum(self.acute_len, prog)
        self.infected[uids] = True
        self.prognosis[uids] = prog
        self.acute_days[uids] = acute
        self.aids_days[uids] = np.maximum(np.minimum(self.aids_len, prog - acute), 0)
        self.days[uids] = self.rng.random(len(uids))*prog if fast_forward else 0.0
        self.stage[uids] = ACUTE
        self.set_stage(uids)
        self.year_infected[uids] = self.t.now('year') - self.days.raw[uids]/365
        return uids

    def set_stage(self, uids):
        """ EMOD's SetStageFromDuration: returns those whose prognosis has run out """
        d, st = self.days.raw[uids], self.stage.raw[uids]
        latent_end = self.prognosis.raw[uids] - self.aids_days.raw[uids]
        early = (st == ACUTE) | (st == LATENT)
        to_aids = early & (d >= latent_end)
        to_latent = (st == ACUTE) & ~to_aids & (d >= self.acute_days.raw[uids])
        self.stage[uids[to_aids]] = AIDS
        self.stage[uids[to_latent]] = LATENT
        return uids[d >= self.prognosis.raw[uids]]

    def cd4(self, uids):
        h = self.inp.hiv
        f = np.clip(self.days.raw[uids]/self.prognosis.raw[uids], 0, 1)
        return (np.sqrt(h.CD4_Post_Infection_Weibull_Scale) - f*(np.sqrt(h.CD4_Post_Infection_Weibull_Scale) - np.sqrt(h.CD4_At_Death_LogLogistic_Scale)))**2

    # %% Steps

    def campaign(self):
        """ Traditional circumcision and the seeding """
        ppl, c = self.sim.people, self.inp.campaign
        year = self.t.now('year')
        if year >= c.circ_year:
            if self.circ_started:
                uids = self.sim.demographics.emoddemog.born
            else:
                uids, self.circ_started = ppl.auids, True
            got = uids[self.rng.random(len(uids)) < c.p_circ[self.node.raw[uids]]] # EMOD draws for everyone, and circumcises only men
            self.circumcised[got[~ppl.female.raw[got]]] = True
        if year >= c.seed_year and not self.seeded:
            self.seeded = True
            high = ppl.auids[ppl.risk.raw[ppl.auids] == 2]
            seeds = self.infect(high[self.rng.random(len(high)) < c.seed_coverage], fast_forward=True)
            self.log_infections(seeds, src=np.full(len(seeds), -1), kind='seed')
        return

    def at_birth(self):
        dem = self.sim.demographics.emoddemog
        if not len(dem.born):
            return
        hit = self.infected.raw[dem.mothers] & (self.rng.random(len(dem.born)) < self.inp.hiv.Maternal_Infection_Transmission_Probability)
        kids = self.infect(dem.born[hit])
        self.log_infections(kids, src=dem.mothers[hit], kind='birth')
        return

    def exposures(self):
        """ The coital acts of this step: who is infected, the index of the edge, and the partner who infected them """
        net, h, ppl = self.sim.networks[self.netname], self.inp.hiv, self.sim.people
        e = net.edges
        none = ss.uids(), np.zeros(0, int), np.zeros(0, int)
        if not len(net):
            return none
        inf = self.infected.raw
        disc = inf[e.p1] != inf[e.p2]
        idx = np.flatnonzero(disc & (e.acts > 0))
        p1, p2 = e.p1[idx], e.p2[idx]
        src = np.where(inf[p1], p1, p2)
        dst = np.where(inf[p1], p2, p1)
        stage_mult = np.array([h.Acute_Stage_Infectivity_Multiplier, 1, h.AIDS_Stage_Infectivity_Multiplier])[self.stage.raw[src]]
        sti = np.maximum(np.where(self.sti.raw[src], h.STI_Coinfection_Transmission_Multiplier, 1), np.where(self.sti.raw[dst], h.STI_Coinfection_Acquisition_Multiplier, 1))
        female = ppl.female.raw[dst]
        acq = np.where(female, np.interp(ppl.emod_age.raw[dst]/365, h.Male_To_Female_Relative_Infectivity_Ages, h.Male_To_Female_Relative_Infectivity_Multipliers),
                       np.where(self.circumcised.raw[dst], 1 - self.inp.campaign.circ_efficacy, 1))
        p_edge = h.Base_Infectivity*stage_mult*sti*acq
        n = e.acts[idx].astype(int)
        act = np.repeat(np.arange(len(idx)), n) # One row per act
        condom = self.rng.random(len(act)) < e.condoms[idx][act]
        p = p_edge[act]*np.where(condom, 1 - h.Condom_Transmission_Blocking_Probability, 1)
        hit = np.flatnonzero(self.rng.random(len(act)) < p)
        if not len(hit):
            return none
        order = hit[np.argsort(self.rng.random(len(hit)))] # Acts in random order; the first for each person infects them
        who = dst[act[order]]
        _, first = np.unique(who, return_index=True)
        k = idx[act[order[first]]]
        return ss.uids(who[first]), k, src[act[order[first]]]

    def step(self):
        ppl = self.sim.people
        self.campaign()
        self.at_birth()
        new, k, src = self.exposures()
        src_stage, src_days = self.stage.raw[src].copy(), self.days.raw[src].copy()
        inf = ppl.auids[self.infected.raw[ppl.auids]]
        self.days[inf] = self.days.raw[inf] + self.dt_days
        dead = self.set_stage(inf)
        if len(dead):
            ppl.request_death(dead, cause=self)
            if self.record:
                self.deaths.append(dict(year=np.full(len(dead), self.t.now('year')), age=ppl.emod_age.raw[dead]/365, female=ppl.female.raw[dead],
                                        node=self.node.raw[dead], years_infected=self.days.raw[dead]/365, cd4=self.cd4(dead)))
        new = self.infect(new)
        self.log_infections(new, src=src, kind=self.sim.networks[self.netname].edges.type[k], src_stage=src_stage, src_days=src_days)
        self.n_new, self.n_dead = len(new), len(dead)
        return

    def log_infections(self, uids, src, kind, src_stage=None, src_days=None):
        np.add.at(self.n_infections.raw, src[src >= 0], 1)
        if not self.record or not len(uids):
            return
        ppl = self.sim.people
        kinds = {'seed': -1, 'birth': -2}
        typ = np.full(len(uids), kinds[kind]) if isinstance(kind, str) else kind
        nan = np.full(len(uids), np.nan)
        self.infections.append(dict(year=np.full(len(uids), self.t.now('year')), type=typ, node=self.node.raw[uids], female=ppl.female.raw[uids],
                                    age=ppl.emod_age.raw[uids]/365, src_stage=nan if src_stage is None else src_stage, src_days=nan if src_days is None else src_days,
                                    src_age=np.where(src >= 0, ppl.emod_age.raw[np.maximum(src, 0)]/365, np.nan)))
        return

    def step_die(self, uids):
        super().step_die(uids)
        self.infected[uids] = False
        return

    def update_results(self):
        super().update_results()
        self.results.new_infections[self.ti] = getattr(self, 'n_new', 0)
        self.results.hiv_deaths[self.ti] = getattr(self, 'n_dead', 0)
        self.results.n_infected[self.ti] = np.count_nonzero(self.infected.raw[self.sim.people.auids])
        return

    def finalize(self):
        super().finalize()
        cat = lambda recs: pd.DataFrame({k: np.concatenate([r[k] for r in recs]) for k in recs[0]}) if recs else pd.DataFrame()
        self.infections, self.deaths = cat(self.infections), cat(self.deaths)
        return
