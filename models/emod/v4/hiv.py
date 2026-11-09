"""
EMOD′ HIV: transmission per coital act, natural history, and the effects of ART, as in EMOD-HIV (docs/emod_inventory.md,
sections 4, 5, 9 (ART), 10 (traditional circumcision), 11 (transmission at birth) and 12). In each step, in EMOD's order:

1. Campaign basics: traditional circumcision from 1975 (everyone then, and each newborn, with a probability by node;
   EMOD draws for both sexes and circumcises only men), and the seeding at SeedYr: a share of everyone at Risk HIGH is
   infected, partway through their course (a uniform fraction of the prognosis, InfectionHIV::SetParameters).
2. Transmission at birth from a mother infected at the birth: probability 0.3, × 0.03334 if she is on ART, else
   × (1 − her PMTCT efficacy) (IndividualHumanHIV::getProbMaternalTransmission).
3. The infectiousness of each act, as at the start of the step (Relationship::Consummate): Base_Infectivity × stage (26
   acute, 10 AIDS, else 1) × ART (from 1 down to its multiplier, linearly over 183 days, then the multiplier; 1 once ART
   fails) × max of the STI multipliers of the two partners × (1 − 0.9) with a condom.
4. Progression (InfectionHIV, SusceptibilityHIV, HIVInterventionsContainer::InfectiousLoopUpdate): time since infection
   (or since starting ART) grows by the step; acute for 2.9 months, AIDS for the last 9 months; death from HIV at the
   end of the prognosis. √CD4 falls linearly from √540.55 to √31.63 over the prognosis; on ART it recovers
   (+15.584 m − 0.2113 m² after m months, at most +287.34, up to its value at infection). People become symptomatic a
   Weibull time (scale 618 days, shape 2) before death, which fires NewlySymptomatic. On ART, failure comes 9 months
   before the death drawn at ART start: the AIDS stage, with no ART effect, and CD4 falling to √31.63 by the death.
5. The campaign's care cascade (cascade.py), after progression, as EMOD updates interventions after infections.
   Steps 1-4 are `step_state()`; the cascade is an intervention, which Starsim runs after it (and after the network);
   step 6 is `step()`.
6. Exposure: each uninfected person's acts, in random order, until one transmits (RelationshipGroups::ExposeToContagion),
   with × (1 − 0.6) for a circumcised man, × the male-to-female multiplier for a woman by her age, and × (1 − PrEP).
   New infections draw their prognosis from the age then (FerrandAgeDependentDistribution): under 15, exponential
   (rate 1.52/year) with probability 0.57, else Weibull (median 16 years, shape 2.7); from 15, Weibull of shape 2 and
   scale 21.182 − 0.2717 × age (age capped at 50); at least 7 days.

ART (`go_on_art`, ARTMortalityTable): the time from starting to death is drawn from annual rates by ART duration bin,
age and CD4 (interpolated), each later bin's draw timed from ART start, as EMOD does. Dropout (`go_off_art`,
InfectionHIV::ApplySuppressionDropout): a new prognosis is drawn from the age then and shortened by the share of the
CD4 decline already done; its AIDS stage starts at the new prognosis's latent end, so someone who recovered little
CD4 on ART can die before reaching it, as in EMOD. People: `people.sti` (risk module), `people.circumcised`.
"""
import numpy as np
import pandas as pd
import starsim as ss

ACUTE, LATENT, AIDS, ON_ART = 0, 1, 2, 3
OFF, RAMP, SUPPRESSED, FAILING = 0, 1, 2, 3 # ART status (HIVInterventionsContainer)
WHO_KAPPA, WHO_LAMBDA = np.array([0.9664, 0.9916, 0.9356]), np.array([0.26596, 0.19729, 0.34721])
CD4_SQ, CD4_LIN, CD4_SAT, CD4_MONTHS = -0.2113, 15.584, 287.3415, 3.073*12 # CD4 reconstitution on ART
IDEAL_MONTH = 30 # IDEALDAYSPERMONTH


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
        F = ss.FloatArr
        self.define_states(
            ss.BoolState('infected', label='Infected'),
            F('days'), F('prognosis'), F('acute_days'), F('aids_start'), # Days since infection (or ART start); days to death; stage ends
            ss.IntArr('stage', default=-1), F('year_infected'), ss.IntArr('n_infections', default=0),
            F('sqrt_cd4'), F('cd4_rate'),                          # √CD4 and its change per day off ART
            F('symp_days'), ss.BoolArr('symptomatic', default=False), # Days before death of becoming symptomatic
            F('who1'), F('who2'), F('who3'),                       # Shares of the prognosis in WHO stages 1-3
            ss.IntArr('art', default=OFF), F('art_mult'), F('supp_timer'), F('fail_timer'), F('art_days', default=-1.0), F('cd4_at_art'),
            ss.BoolArr('ever_tested', default=False), ss.BoolArr('diagnosed', default=False), F('last_test'), F('lowest_cd4', default=np.inf),
            ss.BoolArr('ever_art', default=False), F('pmtct', default=0.0), F('pmtct_until', default=np.inf), F('prep_until', default=-np.inf), F('prep_eff', default=0.0),
            F('p_act'), ss.IntArr('stage0', default=-1), F('days0'), # Per-act infectiousness, stage and days since infection at the start of the step
        )
        self.infections, self.deaths = [], []
        self.seeded = self.circ_started = False
        self.newly_symptomatic = self.new_uids = self.dead_uids = ss.uids()
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.sti = sim.people.add_state(ss.BoolArr('sti', default=False))
        self.circumcised = sim.people.add_state(ss.BoolArr('circumcised', default=False))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 5)
        self.dt_days = self.t.dt_year*365 # The step in days (EMOD's Simulation_Timestep)
        h = self.inp.hiv
        self.acute_len, self.aids_len = h.Acute_Duration_In_Months*365/12, h.AIDS_Duration_In_Months*365/12
        self.sqrt_post, self.sqrt_death = np.sqrt(h.CD4_Post_Infection_Weibull_Scale), np.sqrt(h.CD4_At_Death_LogLogistic_Scale)
        return

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('new_infections', dtype=int, scale=True, label='New infections'),
                            ss.Result('hiv_deaths', dtype=int, scale=True, kind='flow', label='HIV deaths'),
                            ss.Result('n_infected', dtype=int, scale=True, label='People with HIV'),
                            ss.Result('n_on_art', dtype=int, scale=True, label='On ART'))
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

    def draw_symptoms(self, uids):
        h = self.inp.hiv
        het = h.Days_Between_Symptomatic_And_Death_Weibull_Heterogeneity
        self.symp_days[uids] = h.Days_Between_Symptomatic_And_Death_Weibull_Scale*(-np.log(self.rng.random(len(uids))))**het
        self.symptomatic[uids] = False
        return

    def untreated_timers(self, uids):
        """ InfectionHIV::SetupNonSuppressedDiseaseTimers: a prognosis from the age now, its stages, and WHO-stage shares; returns the prognosis """
        prog = self.draw_prognosis(self.sim.people.emod_age.raw[uids]/365)
        acute = np.minimum(self.acute_len, prog)
        aids = np.maximum(np.minimum(self.aids_len, prog - acute), 0)
        self.prognosis[uids] = prog
        self.acute_days[uids] = acute
        self.aids_start[uids] = prog - aids # The end of acute + latent
        self.days[uids] = 0.0
        for k, arr in enumerate((self.who1, self.who2, self.who3)):
            arr[uids] = WHO_LAMBDA[k]*(-np.log(self.rng.random(len(uids))))**(1/WHO_KAPPA[k])
        return prog

    def infect(self, uids, fast_forward=False):
        """ Start infections; seeded ones partway through their course """
        uids = uids[~self.infected.raw[uids]]
        if not len(uids):
            return uids
        self.infected[uids] = True
        prog = self.untreated_timers(uids)
        self.stage[uids] = ACUTE
        self.sqrt_cd4[uids] = self.sqrt_post
        self.cd4_rate[uids] = (self.sqrt_death - self.sqrt_post)/prog
        self.draw_symptoms(uids)
        if fast_forward:
            ff = self.rng.random(len(uids))*prog
            self.days[uids] = ff
            self.sqrt_cd4[uids] = self.sqrt_post + self.cd4_rate.raw[uids]*ff
        self.set_stage(uids)
        self.year_infected[uids] = self.t.now('year') - self.days.raw[uids]/365
        return uids

    def set_stage(self, uids):
        """ EMOD's SetStageFromDuration: returns those whose prognosis has run out """
        d, st = self.days.raw[uids], self.stage.raw[uids]
        early = (st == ACUTE) | (st == LATENT)
        to_aids = early & (d >= self.aids_start.raw[uids])
        to_latent = (st == ACUTE) & ~to_aids & (d >= self.acute_days.raw[uids])
        self.stage[uids[to_aids]] = AIDS
        self.stage[uids[to_latent]] = LATENT
        return uids[d >= self.prognosis.raw[uids]]

    def cd4(self, uids):
        return self.sqrt_cd4.raw[uids]**2

    def prognosis_done(self, uids):
        """ SusceptibilityHIV::GetPrognosisCompletedFraction, from CD4 """
        return np.clip((self.sqrt_post - self.sqrt_cd4.raw[uids])/(self.sqrt_post - self.sqrt_death), 0, 1)

    def who_stage(self, uids):
        """ InfectionHIV::GetWHOStage: 1-4, with the fraction through the stage """
        f = self.prognosis_done(uids)
        w = np.stack([self.who1.raw[uids], self.who2.raw[uids], self.who3.raw[uids]], axis=1)
        cum = np.cumsum(w, axis=1)
        out = np.full(len(uids), 4.0)
        done = np.zeros(len(uids), bool)
        for k in range(3):
            hit = ~done & (cum[:, k] >= f)
            out[hit] = k + 1 + np.minimum((f[hit] - (cum[hit, k] - w[hit, k]))/w[hit, k], 1)
            done |= hit
        return out # 4 for those past the first three stages (MAX_WHO_HIV_STAGE - 1), with no fraction

    def progress(self, uids):
        """ One step of each infection, its CD4 and its ART timers; returns those who die of HIV and those newly symptomatic """
        dt = self.dt_days
        self.days[uids] = self.days.raw[uids] + dt
        dead = self.set_stage(uids)
        art = self.art.raw[uids]
        rec = (art == RAMP) | (art == SUPPRESSED)
        u = uids[rec]
        months = np.minimum((self.art_days.raw[u] + dt)/IDEAL_MONTH, CD4_MONTHS)
        inc = np.minimum(CD4_SQ*months**2 + CD4_LIN*months, CD4_SAT)
        low = self.sqrt_cd4.raw[u] < self.sqrt_post
        self.sqrt_cd4[u[low]] = np.minimum(np.sqrt(self.cd4_at_art.raw[u[low]] + inc[low]), self.sqrt_post)
        u = uids[~rec]
        self.sqrt_cd4[u] = np.maximum(self.sqrt_cd4.raw[u] + self.cd4_rate.raw[u]*dt, self.sqrt_death)
        new_symp = u[~self.symptomatic.raw[u] & (self.prognosis.raw[u] - self.days.raw[u] <= self.symp_days.raw[u])]
        self.symptomatic[new_symp] = True
        self.symp_days[new_symp] = -np.inf
        # ART timers
        u = uids[art == RAMP]
        self.supp_timer[u] = self.supp_timer.raw[u] - dt
        self.art[u[self.supp_timer.raw[u] <= 0]] = SUPPRESSED
        u = uids[self.fail_timer.raw[uids] >= 0]
        self.fail_timer[u] = self.fail_timer.raw[u] - dt
        fail = u[self.fail_timer.raw[u] <= 0]
        if len(fail):
            failing = np.minimum(self.aids_len, self.prognosis.raw[fail])
            self.art[fail] = FAILING
            self.supp_timer[fail] = -1.0
            self.stage[fail] = AIDS
            self.cd4_rate[fail] = np.where(failing > 0, (self.sqrt_death - self.sqrt_cd4.raw[fail])/np.maximum(failing, 1e-9), 0)
            self.draw_symptoms(fail)
        u = uids[self.art_days.raw[uids] >= 0]
        self.art_days[u] = self.art_days.raw[u] + dt
        return dead, new_symp[~np.isin(new_symp, dead)]

    # %% ART, as the cascade gives it

    def art_death_days(self, uids, table):
        """ ARTMortalityTable::ComputeDurationFromEnrollmentToArtAidsDeath """
        dur_bins, age_bins, cd4_bins = (np.array(table[k], float) for k in ('ART_Duration_Days_Bins', 'Age_Years_Bins', 'CD4_Count_Bins'))
        mort = np.array(table['MortalityTable'], float) # Duration × age × CD4
        bin_ = lambda v, b: np.minimum(np.searchsorted(b, v, side='left'), len(b) - 1) # ReportUtilities::GetBinIndex
        on = self.art_days.raw[uids]
        di = bin_(on, dur_bins)
        ai = bin_(self.sim.people.emod_age.raw[uids]/365, age_bins)
        cd4 = self.cd4(uids)
        ci = bin_(cd4, cd4_bins)
        days = np.full(len(uids), np.inf)
        todo = np.ones(len(uids), bool)
        while todo.any():
            k = np.flatnonzero(todo)
            rows = mort[di[k], ai[k]]
            lo = rows[np.arange(len(k)), np.maximum(ci[k] - 1, 0)]
            hi = rows[np.arange(len(k)), ci[k]]
            frac = np.where(ci[k] > 0, (cd4[k] - cd4_bins[np.maximum(ci[k] - 1, 0)])/np.maximum(cd4_bins[ci[k]] - cd4_bins[np.maximum(ci[k] - 1, 0)], 1e-9), 0)
            rate = np.where(cd4[k] >= cd4_bins[-1], rows[:, -1], np.where(ci[k] == 0, rows[:, 0], lo + frac*(hi - lo)))
            d = np.where(rate > 0, -np.log(self.rng.random(len(k)))/np.where(rate > 0, rate, 1)*365, np.inf)
            days[k] = d
            di[k] += 1
            todo[k] = (di[k] < len(dur_bins)) & (d > dur_bins[di[k] - 1] - on[k])
        return days

    def go_on_art(self, uids, table, mult, days_to_suppress):
        """ HIVInterventionsContainer::GoOnART, for those infected """
        uids = uids[self.infected.raw[uids]]
        if not len(uids):
            return
        death = self.art_death_days(uids, table)
        start = uids[self.art.raw[uids] == OFF]
        self.art_days[start] = 0.0
        self.cd4_at_art[start] = self.cd4(start)
        self.ever_art[start] = True
        self.art_mult[uids] = mult
        self.art[uids] = RAMP
        self.prognosis[uids] = death
        self.days[uids] = 0.0
        self.stage[uids] = ON_ART
        self.supp_timer[uids] = days_to_suppress
        failing = np.minimum(self.aids_len, death)
        early = np.minimum(days_to_suppress, death - failing)
        self.fail_timer[uids] = np.maximum(death - early - failing, 0) + early
        return

    def go_off_art(self, uids):
        """ HIVInterventionsContainer::GoOffART and InfectionHIV::ApplySuppressionDropout """
        uids = uids[self.infected.raw[uids] & (self.art.raw[uids] != OFF)]
        if not len(uids):
            return
        self.art[uids] = OFF
        self.supp_timer[uids] = -1.0
        self.fail_timer[uids] = -1.0
        self.art_days[uids] = -1.0
        p0 = self.untreated_timers(uids)
        done = self.prognosis_done(uids)
        left = p0*(1 - done)
        self.cd4_rate[uids] = np.where(left > 0, (self.sqrt_death - self.sqrt_cd4.raw[uids])/np.maximum(left, 1e-9), 0)
        self.draw_symptoms(uids)
        self.prognosis[uids] = left
        self.aids_start[uids] = self.aids_start.raw[uids] - self.acute_days.raw[uids] # Acute is skipped: AIDS starts after the new prognosis's latent length
        self.acute_days[uids] = 0.0
        self.stage[uids] = LATENT
        self.set_stage(uids)
        return

    def on_art(self, uids):
        return self.art.raw[uids] != OFF

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
            self.step_new.append(seeds)
            self.log_infections(seeds, src=np.full(len(seeds), -1), kind='seed')
        return

    def at_birth(self):
        dem, h = self.sim.demographics.emoddemog, self.inp.hiv
        if not len(dem.born):
            return
        m = dem.mothers
        p = np.where(self.on_art(m), h.Maternal_Infection_Transmission_Probability*h.Maternal_Transmission_ART_Multiplier,
                     h.Maternal_Infection_Transmission_Probability*(1 - np.where(self.pmtct_until.raw[m] >= self.t.now('year'), self.pmtct.raw[m], 0)))
        hit = self.infected.raw[m] & (self.rng.random(len(m)) < p)
        kids = self.infect(dem.born[hit])
        self.step_new.append(kids)
        self.log_infections(kids, src=m[hit], kind='birth')
        return

    def infectiousness(self, uids):
        """ Per-act transmission probability of infected people before acquisition factors (UpdateInfectiousnessSTI) """
        h = self.inp.hiv
        st, art = self.stage.raw[uids], self.art.raw[uids]
        stage_mult = np.select([st == ACUTE, st == AIDS], [h.Acute_Stage_Infectivity_Multiplier, h.AIDS_Stage_Infectivity_Multiplier], 1.0)
        m = self.art_mult.raw[uids]
        supp = np.select([art == RAMP, art == SUPPRESSED], [1 - (1 - m)/183*(183 - self.supp_timer.raw[uids]), m], 1.0)
        return h.Base_Infectivity*stage_mult*supp

    def acts(self):
        """ The coital acts of this step between an infected and an uninfected partner, with the infected side's part of their probability """
        net, h = self.sim.networks[self.netname], self.inp.hiv
        e = net.edges
        if not len(net):
            return None
        inf = self.infected.raw
        idx = np.flatnonzero((inf[e.p1] != inf[e.p2]) & (e.acts > 0))
        p1, p2 = e.p1[idx], e.p2[idx]
        src = np.where(inf[p1], p1, p2)
        dst = np.where(inf[p1], p2, p1)
        sti = np.maximum(np.where(self.sti.raw[src], h.STI_Coinfection_Transmission_Multiplier, 1), np.where(self.sti.raw[dst], h.STI_Coinfection_Acquisition_Multiplier, 1))
        return dict(idx=idx, src=src, dst=dst, p=self.p_act.raw[src]*sti, stage=self.stage0.raw[src], days=self.days0.raw[src], type=e.type[idx])

    def exposures(self, a):
        """ Who is infected by the acts, and by which """
        none = ss.uids(), np.zeros(0, int)
        if a is None or not len(a['idx']):
            return none
        h, ppl, e = self.inp.hiv, self.sim.people, self.sim.networks[self.netname].edges
        dst = a['dst']
        ok = ppl.alive.raw[dst] & ~self.infected.raw[dst]
        female = ppl.female.raw[dst]
        acq = np.where(female, np.interp(ppl.emod_age.raw[dst]/365, h.Male_To_Female_Relative_Infectivity_Ages, h.Male_To_Female_Relative_Infectivity_Multipliers),
                       np.where(self.circumcised.raw[dst], 1 - self.inp.campaign.circ_efficacy, 1))
        acq = acq*np.where(self.prep_until.raw[dst] >= self.t.now('year'), 1 - self.prep_eff.raw[dst], 1)*ok
        n = e.acts[a['idx']].astype(int)
        act = np.repeat(np.arange(len(n)), n)
        condom = self.rng.random(len(act)) < e.condoms[a['idx']][act]
        p = (a['p']*acq)[act]*np.where(condom, 1 - h.Condom_Transmission_Blocking_Probability, 1)
        hit = np.flatnonzero(self.rng.random(len(act)) < p)
        if not len(hit):
            return none
        order = hit[np.argsort(self.rng.random(len(hit)))] # Acts in random order; the first for each person infects them
        who = dst[act[order]]
        _, first = np.unique(who, return_index=True)
        return ss.uids(who[first]), act[order[first]]

    def step_state(self):
        """ Steps 1-4: campaign basics, births, the infectiousness of the step's acts, and progression """
        super().step_state()
        ppl = self.sim.people
        self.step_new = []
        self.campaign()
        self.at_birth()
        inf = ppl.auids[self.infected.raw[ppl.auids]]
        self.p_act[inf] = self.infectiousness(inf) # As at the start of the step, for transmission in step()
        self.stage0[inf] = self.stage.raw[inf]
        self.days0[inf] = self.days.raw[inf]
        dead, self.newly_symptomatic = self.progress(inf)
        if len(dead):
            ppl.request_death(dead, cause=self)
            if self.record:
                self.deaths.append(dict(year=np.full(len(dead), self.t.now('year')), age=ppl.emod_age.raw[dead]/365, female=ppl.female.raw[dead],
                                        node=self.node.raw[dead], years_infected=(self.t.now('year') - self.year_infected.raw[dead]), cd4=self.cd4(dead),
                                        on_art=self.on_art(dead)))
        self.dead_uids = dead
        return

    def step(self):
        """ Step 6: exposure, after the cascade (an intervention, so it runs between step_state and step) """
        a = self.acts()
        dead = self.dead_uids
        new, k = self.exposures(a)
        new = self.infect(new)
        self.sexual = new, (a['src'][k] if a is not None else np.zeros(0, int)) # This step's sexual transmissions: recipients and sources
        if a is not None:
            self.log_infections(new, src=a['src'][k], kind=a['type'][k], src_stage=a['stage'][k], src_days=a['days'][k])
        self.n_new, self.n_dead = len(new), len(dead)
        self.new_uids = ss.uids(np.concatenate(self.step_new + [new])) # Every infection of the step, as EMOD's NewInfectionEvent
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

    def on_death(self, uids):
        super().on_death(uids)
        self.infected[uids] = False
        self.art[uids] = OFF
        return

    def record_flows(self):
        super().record_flows()
        self.results.new_infections[self.ti] = getattr(self, 'n_new', 0)
        self.results.hiv_deaths[self.ti] = getattr(self, 'n_dead', 0)
        return

    def record_stocks(self):
        super().record_stocks()
        au = self.sim.people.auids
        self.results.n_infected[self.ti] = np.count_nonzero(self.infected.raw[au])
        self.results.n_on_art[self.ti] = np.count_nonzero(self.art.raw[au] != OFF)
        return

    def finalize(self):
        super().finalize()
        cat = lambda recs: pd.DataFrame({k: np.concatenate([r[k] for r in recs]) for k in recs[0]}) if recs else pd.DataFrame()
        self.infections, self.deaths = cat(self.infections), cat(self.deaths)
        return
