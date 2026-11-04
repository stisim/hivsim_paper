"""
EMOD′ demography: EMOD-HIV Malawi's population, births and background deaths (docs/emod_inventory.md, section 1).

- Initial population: InitialPopulation × x_Base_Population agents in each of the 3 nodes, sex 50/50, ages drawn
  from the AgeDistribution CDF by linear interpolation (Node::PopulateFromDemographics, calculateInitialAge). Women of
  childbearing age start pregnant with probability rate × 280 days, with the time left uniform
  (conditionallyInitializePregnancy).
- Pregnancies (INDIVIDUAL_PREGNANCIES_BY_AGE_AND_YEAR, Node::considerPregnancyForIndividual): a woman aged over 14 and
  under 45 who is not pregnant conceives with probability r/(1 − 280r) × dt, r the FertilityDistribution rate per day
  at her age and the year; the birth is 280 − U×dt days later. A woman who gives birth does not conceive in the same step.
- Background deaths (NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER, IndividualHuman::CheckVitalDynamics):
  probability rate × dt, the rate from the table for her or his sex, age and the year. The tables hold deaths per
  person per year, but EMOD reads them with ResultScaleFactor 1/365,000, so every rate is 1000 times too low; that is
  copied, and `fix_mortality=True` reads them as intended.
- Tables are EMOD's 2-axis distributions: linear between groups, constant beyond them (pickDistributionValue).
- No migration: each person stays in their node (`people.node`, 0-2), and newborns are in their mother's.
- Ages as EMOD keeps them: in days, in single precision (IndividualHuman::m_age is a float), increased by the step
  once in each step, after the pair formation and before deaths, pregnancies and reports. Adding 30.4167 days in single
  precision drifts, so after whole steps people are just below some 5-year edges (5, 15, 20, 25, 55-65) and just
  above others (10, 30-50); a monthly cohort falls on the other side of them from exact arithmetic. Newborns are aged 0
  until the step after their birth. `people.age` holds the age before this step's increase (what EMOD's pair
  formation sees), and `people.emod_age` the age in days after it (what its deaths, pregnancies, reports and debut
  events see).

The births of each step are in `born` (newborns) and `mothers`, for diseases that pass to the child at birth.
"""
import numpy as np
import starsim as ss

GESTATION = 280 # Days (DAYSPERWEEK * WEEKS_FOR_GESTATION)
f32 = np.float32


def lookup(tab, age, year):
    """ An EMOD 2-axis demographics distribution at a year, for each age (years) """
    col = np.array([np.interp(year, tab.year, row) for row in tab.vals])
    return np.interp(age, tab.age, col)


class EmodDemography(ss.Demographics):
    """
    Args:
        inp (objdict): EMOD's inputs for one parameter set, from `inputs.load()`
        fix_mortality (bool): read the mortality tables as deaths per person per year, as they are labelled (default: as EMOD does)
    """
    def __init__(self, inp, fix_mortality=False, name='emoddemog', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp, self.fix_mortality = inp, fix_mortality
        self.define_states(ss.FloatArr('preg_days')) # Days left in a pregnancy
        self.born = self.mothers = ss.uids()
        self.n_births = self.n_deaths = 0
        self.twelve_weeks = ss.uids()
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.emod_age = sim.people.add_state(ss.FloatArr('emod_age', default=0.0)) # Days, after this step's ageing
        self.rng = np.random.default_rng(sim.pars.rand_seed + 4)
        self.dt_days = self.t.dt_year*365 # The step in days (EMOD's Simulation_Timestep)
        return

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('births', dtype=int, scale=True, label='Births'),
                            ss.Result('deaths', dtype=int, scale=True, label='Background deaths'))
        return

    def init_post(self):
        super().init_post()
        inp, ppl = self.inp, self.sim.people
        uids = ppl.auids
        counts = inp.init_pop.astype(int)
        assert counts.sum() == len(uids), f'the sim needs n_agents={counts.sum()}'
        self.node[uids] = np.repeat(np.arange(3), counts)
        days = np.interp(self.rng.random(len(uids)), inp.age_cdf.cdf, inp.age_cdf.age)*365
        self.emod_age[uids] = days.astype(f32)
        ppl.age[uids] = self.emod_age.raw[uids]/365
        moms = uids[self.can_conceive(uids)]
        rate = lookup(inp.fertility, ppl.age.raw[moms], self.t.now('year'))/365
        preg = moms[self.rng.random(len(moms)) < rate*GESTATION]
        self.preg_days[preg] = self.rng.random(len(preg))*GESTATION
        return

    def can_conceive(self, uids, age=None):
        ppl = self.sim.people
        age = ppl.age.raw[uids] if age is None else age
        return ppl.female.raw[uids] & (age > 14) & (age < 45) & np.isnan(self.preg_days.raw[uids])

    def step(self):
        inp, ppl, dt = self.inp, self.sim.people, self.dt_days
        year = self.t.now('year')
        uids = ppl.auids
        female = ppl.female.raw[uids]

        # Ageing: what the rest of this step sees before it (people.age) and after it (people.emod_age, used here)
        before = self.emod_age.raw[uids]
        ppl.age[uids] = before/365
        self.emod_age[uids] = (before.astype(f32) + f32(dt)).astype(f32)
        age = self.emod_age.raw[uids]/365

        # Background deaths
        rate = np.where(female, lookup(inp.mortality[1], age, year), lookup(inp.mortality[0], age, year))/365*(1000 if self.fix_mortality else 1)
        dead = uids[self.rng.random(len(uids)) < rate*dt]
        ppl.request_death(dead, cause=self)

        # Births to those whose pregnancy ends, then new pregnancies
        preg = uids[~np.isnan(self.preg_days.raw[uids])]
        self.preg_days[preg] = self.preg_days.raw[preg] - dt
        left = self.preg_days.raw[preg]
        self.twelve_weeks = preg[(left - dt < GESTATION - 84) & (left >= GESTATION - 84)] # IndividualHumanHIV::UpdatePregnancy's TwelveWeeksPregnant
        born = preg[self.preg_days.raw[preg] <= 0]
        self.preg_days[born] = np.nan
        ok = self.can_conceive(uids, age) & ~np.isin(uids, born)
        moms = uids[ok]
        r = lookup(inp.fertility, age[ok], year)/365
        r = r/(1 - r*GESTATION) # So that the birth rate is r, allowing for the time pregnant
        new = moms[self.rng.random(len(moms)) < r*dt]
        self.preg_days[new] = GESTATION - self.rng.random(len(new))*dt
        self.mothers = born
        self.born = ppl.grow(len(born), parent=born, age=0) if len(born) else ss.uids()
        self.node[self.born] = self.node.raw[born]
        self.emod_age[self.born] = 0.0
        self.n_births, self.n_deaths = len(born), len(dead)
        return

    def update_results(self):
        super().update_results()
        self.results.births[self.ti] = self.n_births
        self.results.deaths[self.ti] = self.n_deaths
        return
