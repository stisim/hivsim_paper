"""
A standalone harness for the EMOD′ network: the network and risk modules on an EMOD-like population, with no HIV.

The demography here is minimal, only enough to give the network EMOD's population: EMOD's initial age distribution
and node sizes (1% scale), and its pregnancy and background mortality tables, read as EMOD reads them (so with the
mortality rates 1000 times too low). It is not the EMOD′ demography component: no HIV deaths, no migration, and
no check beyond the population size and age structure.

    sim = make_sim(pid=0, seed=1); sim.run()
"""
import numpy as np
import pandas as pd
import starsim as ss

import inputs
from network import EmodNetwork, KEYS
from risk import EmodRisk

GESTATION = 280 # Days (DAYSPERWEEK * WEEKS_FOR_GESTATION)
AGES = np.arange(0, 70, 5) # Lower edges of the report's age bins, as in the EMOD run


def lookup(tab, age, year):
    """ EMOD's 2-axis demographics distribution at a year, for each age: linear between groups, constant beyond them """
    col = np.array([np.interp(year, tab.year, row) for row in tab.vals])
    return np.interp(age, tab.age, col)


class EmodDemography(ss.Demographics):
    """ EMOD's initial population, pregnancies by age and year, and background deaths by age, sex and year, in 3 nodes """
    def __init__(self, inp, name='emoddemog', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp = inp
        self.define_states(ss.FloatArr('preg_days')) # Days left in a pregnancy
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 4)
        return

    def init_post(self):
        super().init_post()
        inp, ppl, dt = self.inp, self.sim.people, inputs.DT
        uids = ppl.auids
        counts = inp.init_pop.astype(int)
        assert counts.sum() == len(uids)
        self.node[uids] = np.repeat(np.arange(3), counts)
        ppl.age[uids] = np.interp(self.rng.random(len(uids)), inp.age_cdf.cdf, inp.age_cdf.age)
        moms = uids[self.mothers(uids)] # Pregnant at the start: for a share of the gestation, with the time left uniform
        rate = lookup(inp.fertility, ppl.age.raw[moms], self.t.now('year'))/365
        preg = moms[self.rng.random(len(moms)) < rate*GESTATION]
        self.preg_days[preg] = self.rng.random(len(preg))*GESTATION
        return

    def mothers(self, uids):
        ppl = self.sim.people
        age = ppl.age.raw[uids]
        return ppl.female.raw[uids] & (age > 14) & (age < 45) & np.isnan(self.preg_days.raw[uids])

    def step(self):
        inp, ppl, dt = self.inp, self.sim.people, inputs.DT
        year = self.t.now('year')
        uids = ppl.auids
        age, female = ppl.age.raw[uids], ppl.female.raw[uids]

        # Deaths (EMOD: rate per day times the step)
        rate = np.where(female, lookup(inp.mortality[1], age, year), lookup(inp.mortality[0], age, year))/365
        ppl.request_death(uids[self.rng.random(len(uids)) < rate*dt], cause=self)

        # Births to those whose pregnancy ends, then new pregnancies
        preg = uids[~np.isnan(self.preg_days.raw[uids])]
        self.preg_days[preg] = self.preg_days.raw[preg] - dt
        born = preg[self.preg_days.raw[preg] <= 0]
        self.preg_days[born] = np.nan
        moms = uids[self.mothers(uids)]
        moms = moms[~np.isin(moms, born)]
        r = lookup(inp.fertility, ppl.age.raw[moms], year)/365
        r = r/(1 - r*GESTATION) # So that the birth rate is r, allowing for the time pregnant
        new = moms[self.rng.random(len(moms)) < r*dt]
        self.preg_days[new] = GESTATION - self.rng.random(len(new))*dt
        if len(born):
            kids = ppl.grow(len(born), parent=born, age=0)
            self.node[kids] = self.node.raw[born]
        return


class Snapshots(ss.Analyzer):
    """ Once a year, at the times of EMOD's ReportHIVByAgeAndGender: people and their relationships by node, sex, Risk and age """
    def __init__(self, net='emodnet', name='snapshots', **kwargs):
        super().__init__(name=name, **kwargs)
        self.netname = net
        self.rows = []
        return

    def step(self):
        if (self.ti - 5) % 12: # EMOD reports at 1960.92, 1961.92, ...
            return
        ppl, net = self.sim.people, self.sim.networks[self.netname]
        uids = ppl.auids
        n = net.stack('n', uids)
        q = net.stack('q', uids)
        mx = net.stack('max', uids)
        debuted = ppl.age.raw[uids] >= net.debut.raw[uids]
        avail = available(n, q, mx, net.flags.raw[uids].astype(np.int64)) & debuted[:, None]
        df = pd.DataFrame(dict(node=net.node.raw[uids] + 1, female=ppl.female.raw[uids].astype(int), risk=net.risk.raw[uids],
                               age=AGES[np.searchsorted(AGES, ppl.age.raw[uids], side='right') - 1], pop=1, debuted=debuted,
                               concurrent=n.sum(1) > 1, partners=n.sum(1)))
        for t, k in enumerate(KEYS):
            df[f'cur_{k}'] = n[:, t] > 0
            df[f'num_{k}'] = n[:, t]
            df[f'que_{k}'] = q[:, t]
            df[f'avail_{k}'] = avail[:, t]
        out = df.groupby(['node', 'female', 'risk', 'age']).sum().reset_index()
        out.insert(0, 'year', round(self.t.now('year') + 1/12, 2))
        self.rows.append(out)
        return

    def finalize(self):
        super().finalize()
        self.df = pd.concat(self.rows, ignore_index=True)
        return


def available(n, q, mx, flags):
    from network import available_all
    return available_all(n, q, mx, flags)


def make_sim(pid=0, seed=1, stop=2025, record=True, **kwargs):
    """ The network on EMOD's population for parameter set `pid` """
    inp = inputs.load(pid)
    return ss.Sim(n_agents=int(inp.init_pop.sum()), start=inp.start, stop=stop, dt=1/12, rand_seed=seed, verbose=0,
                  demographics=EmodDemography(inp), networks=EmodNetwork(inp, record=record), interventions=EmodRisk(inp),
                  analyzers=Snapshots(), **kwargs)
