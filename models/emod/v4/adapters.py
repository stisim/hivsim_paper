"""
EMOD′ components on HIVsim's slots (docs/swap_design.md, section 5, steps 6 and 7): STIsim's demography (`ss.Pregnancy`,
`ss.Deaths` and the 1985 age structure, from HIVsim Malawi's data) and, for EmodHIV alone, its network and risk
(`sti.StructuredSexual`).

- `EmodPeople` (D bridge, named `emoddemog` so EMOD′ modules find it where they find `EmodDemography`): the people
  fields EMOD′ modules read, made from STIsim's demography. `people.node`: drawn at the start from EMOD's node shares
  (its InitialPopulation by node), and a newborn's is its mother's. `people.emod_age`: EMOD's age in days after the
  step's ageing, `(age + dt)*365` (Starsim ages people at the end of the step), 0 in the step of birth. `born` and
  `mothers`: the step's live births and their mothers, from the pregnancy's deliveries (EmodHIV's transmission at birth
  and traditional circumcision of newborns read them). `twelve_weeks`: women reaching 12 weeks of gestation this step.
  With `risk_net`, it also reads EMOD's R from HIVsim's (step 6): `people.risk` HIGH for current FSW and clients, LOW
  for `risk_group` 0 and MEDIUM otherwise, each step; `people.sti` drawn once at the network's debut (on day 1 for those
  past it), with EMOD's probabilities by that Risk.
- `EmodHIVOnSti` (step 6): EmodHIV reading `StructuredSexual`'s edges. STIsim's `acts` are already whole acts per step
  (drawn once per partnership), taken as they are; `condoms` is the chance an act has a condom, as EMOD′'s; types are
  mapped stable → MARITAL, casual → INFORMAL, one-time → TRANSITORY, sex work → COMMERCIAL (only the infection log
  reads them).
- `EmodNetworkAt` and `EmodRiskAt` (step 7): EmodNetwork and EmodRisk starting in a population that has no
  sexual history. At its first step the network draws debut ages, then EmodRisk gives those already past debut the
  sex work and STI flag they would have had (`EmodRiskAt.start_history()`), then the network runs a burn-in as a
  warm-up: 5 years of its own steps (end, form, match, acts) with time, ages and the rest of the sim frozen, so the
  queue-entry rates are set from the starting population and frozen before the hybrid's first step, as EMOD's are
  after its burn-in. The warm-up's relationships carry on into the run. EMOD's own burn-in length (15 years) is too
  long with ages frozen: long relationships pile up in the frozen age groups and the rates set from that population
  stay wrong after it (analysis 12: relationships by type, sex and age within 8% of EMOD′'s at 5 years, 17% at 15). The network reseeds its Numba stream from
  its own generator at each step (`own_stream`), so other Numba code drawing from the shared global stream can't move
  its draws.

`make_sim(pid, seed, hiv=True, network='sti'|'emod')` builds the hybrids; see `models/emod/v4/README.md`.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.data import loaders as stidata

import inputs
from hiv import EmodHIV
from network import EmodNetwork, nb_seed
from risk import EmodRisk, LOW, MEDIUM, HIGH

REPO = Path(__file__).resolve().parents[3]
DATADIR = REPO / 'models' / 'hivsim' / 'malawi' / 'data'
TYPE_MAP = np.array([2, 1, 0, 3]) # STIsim's stable, casual, onetime, sw → EMOD's MARITAL, INFORMAL, TRANSITORY, COMMERCIAL


def sti_demography(n_agents, start=1985, location='malawi', datafolder=DATADIR):
    """ STIsim's demography for a location, as `sti.Sim(demographics=location, datafolder=...)` builds it, as modules
    (so other demographic modules can follow them): returns `demographics`, `people` and `total_pop` for the sim """
    columns = dict(year='Time', age='AgeStart', value='Value')
    preg = ss.Pregnancy(fertility_rate=stidata.get_rates(location, 'asfr', datafolder), columns=columns)
    deaths = ss.Deaths(death_rate=stidata.get_rates(location, 'death', datafolder), rate_per=1, columns=dict(year='Time', sex='Sex', age='AgeStart', value='Value'))
    age = stidata.get_age_distribution(location, year=start, datafolder=datafolder)
    total_pop = int(age.value.sum())
    age['value'] /= age.value.sum()
    return dict(demographics=[preg, deaths], people=ss.People(n_agents, age_data=age), total_pop=total_pop)


class EmodPeople(ss.Demographics):
    """
    EMOD′'s people fields on STIsim's demography (see the module docstring); after the pregnancy and deaths modules

    Args:
        inp (objdict): EMOD's inputs, from `inputs.load()` (node shares, STI probabilities)
        preg (str): the name of the pregnancy module
        risk_net (str): the name of a `StructuredSexual` network to read EMOD's Risk and STI flag from (None if `EmodRisk` sets them)
        warm (str): the name of an `EmodNetworkAt` to warm up in the first step
    """
    def __init__(self, inp, preg='pregnancy', risk_net=None, warm=None, name='emoddemog', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp, self.pregname, self.risk_net, self.warm = inp, preg, risk_net, warm
        self.define_states(ss.BoolArr('sti_drawn', default=False))
        self.born = self.mothers = self.twelve_weeks = ss.uids()
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.emod_age = sim.people.add_state(ss.FloatArr('emod_age', default=0.0))
        if self.risk_net:
            self.risk = sim.people.add_state(ss.IntArr('risk', default=-1))
            self.sti = sim.people.add_state(ss.BoolArr('sti', default=False))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 6)
        self.dt_years = float(self.t.dt_year)
        return

    def init_post(self):
        super().init_post()
        uids = self.sim.people.auids
        share = self.inp.init_pop/self.inp.init_pop.sum()
        self.node[uids] = self.rng.choice(3, size=len(uids), p=share)
        return

    def step(self):
        ppl, preg = self.sim.people, self.sim.demographics[self.pregname]
        uids = ppl.auids
        age = ppl.age.raw[uids]
        born = uids[(age == 0) & (ppl.parent.raw[uids] >= 0)]
        moms = ss.uids(ppl.parent.raw[born])
        new = preg.ti_delivered.raw[moms] == self.ti # Born in this step's deliveries
        self.born, self.mothers = born[new], moms[new]
        self.node[self.born] = self.node.raw[self.mothers]
        self.emod_age[uids] = np.where(age == 0, 0, (age + self.dt_years)*365)
        p = preg.P.uids
        weeks = (self.ti - preg.ti_P.raw[p])*self.dt_years*52.18 # Since conception, at the end of the step
        self.twelve_weeks = p[(weeks >= 12) & (weeks - self.dt_years*52.18 < 12)]
        if self.risk_net:
            self.read_risk(uids)
        if self.warm is not None and self.ti == 0:
            self.sim.networks[self.warm].warm_up()
        return

    def read_risk(self, uids):
        """ EMOD's Risk and STI flag from StructuredSexual's risk group and sex work """
        net = self.sim.networks[self.risk_net]
        sw = net.fsw.raw[uids] | net.client.raw[uids]
        self.risk[uids] = np.where(sw, HIGH, np.where(net.risk_group.raw[uids] == 0, LOW, MEDIUM))
        new = uids[~self.sti_drawn.raw[uids] & (self.sim.people.age.raw[uids] >= net.debut.raw[uids])]
        self.sti_drawn[new] = True
        self.sti[new[self.rng.random(len(new)) < self.inp.campaign.p_sti[self.risk.raw[new]]]] = True
        return


class EmodHIVOnSti(EmodHIV):
    """ EmodHIV on `StructuredSexual` (see the module docstring) """
    def __init__(self, inp, net='structuredsexual', **kwargs):
        super().__init__(inp, net=net, **kwargs)
        return

    def contacts(self):
        e = self.sim.networks[self.netname].edges
        if not len(e.p1):
            return None
        return sc.objdict(p1=e.p1, p2=e.p2, acts=np.asarray(e.acts).astype(int), condoms=e.condoms, type=TYPE_MAP[np.asarray(e.edge_type).astype(int)])


class EmodRiskAt(EmodRisk):
    """ EmodRisk starting in a population that has no history: see `start_history()` """
    def start_history(self):
        """ For those past debut at the start (as the network has drawn it): the STI flag (EMOD's day 1), and the sex
        work of those whose debut was in or after the year the campaign's listeners start, as if from their debut """
        ppl, sw, net = self.sim.people, self.inp.sexwork, self.sim.networks[self.netname]
        year = self.t.now('year')
        au = ppl.auids
        since = ppl.emod_age.raw[au]/365 - net.debut.raw[au]
        past = au[since >= 0]
        since = since[since >= 0]
        self.debuted[past] = True
        self.sti[past[self.rng.random(len(past)) < self.inp.campaign.p_sti[self.risk.raw[past]]]] = True
        debut_year = year - since
        female = ppl.female.raw[past]
        p = sw.p_ever[self.node.raw[past], female.astype(int)]
        ever = (debut_year >= sw.start_year) & (self.rng.random(len(past)) < p)
        u, f, dy = past[ever], female[ever], debut_year[ever]
        uptake = dy + self.rng.random(len(u))*sw.delay_max[f.astype(int)]/365
        days = np.where(f, sw.dropout[1]['Lambda']*(-np.log(self.rng.random(len(u))))**(1/sw.dropout[1]['Kappa']),
                        sw.dropout[0]['Min'] + self.rng.random(len(u))*(sw.dropout[0]['Max'] - sw.dropout[0]['Min']))
        dropout = uptake + days/365
        later, now, done = uptake > year, (uptake <= year) & (dropout > year), dropout <= year
        self.t_uptake[u[later]] = uptake[later]
        self.risk[u[now]] = HIGH
        self.t_dropout[u[now]] = dropout[now]
        self.risk[u[done]] = MEDIUM
        return


class EmodNetworkAt(EmodNetwork):
    """
    EmodNetwork starting in a population with no relationships: its burn-in as a warm-up before the first step (see
    the module docstring)

    Args:
        risk (str): the name of the `EmodRiskAt` module that gives the starting population its sex work
        warmup (float): days of warm-up (default 5 years; EMOD's own burn-in is `inp.burnin`)
        own_stream (bool): reseed the Numba stream from the module's own generator at each step
    """
    def __init__(self, inp, risk='emodrisk', warmup=5*365, own_stream=True, **kwargs):
        super().__init__(inp, **kwargs)
        self.riskname, self.warmup, self.own_stream = risk, warmup, own_stream
        self.warmed = False
        return

    def warm_up(self):
        """ The starting population's history and the warm-up; `EmodPeople` calls it in the first step, before the diseases' step (EmodHIV seeds by Risk then) """
        if self.warmed:
            return
        self.warmed = True
        if self.own_stream:
            nb_seed(int(self.rng.integers(2**31)))
        self.set_people() # Debut ages, so the risk module can give the starting population its history
        self.sim.interventions[self.riskname].start_history()
        n = int(np.ceil((self.inp.burnin if self.warmup is None else self.warmup)/self.dt_days))
        record, self.record = self.record, False
        for _ in range(n):
            super().step()
        self.record = record
        self.burnin_left = 0 # Rates frozen from here, as after EMOD's burn-in
        for k in self.results.keys(): # The warm-up's counts aren't the first step's
            if k.startswith('new_'):
                self.results[k][self.ti] = 0
        return

    def step(self):
        if self.own_stream:
            nb_seed(int(self.rng.integers(2**31)))
        super().step()
        return


def make_sim(pid=0, seed=1, n_agents=10_000, start=1985, stop=2025, hiv=True, network='sti', record=True, condoms=True, **kwargs):
    """
    A hybrid of EMOD′ components on HIVsim's demography (STIsim's, with HIVsim Malawi's data), monthly from `start`

    Args:
        pid (int): EMOD′ parameter set
        network (str): 'sti' for HIVsim's network and risk (`StructuredSexual`, with Malawi's condom data if `condoms`; step 6), 'emod' for EMOD′'s (step 7)
        hiv (bool): with EmodHIV (no care cascade)
    """
    inp = inputs.load(pid)
    dem = sti_demography(n_agents, start=start)
    mods = {}
    if network == 'sti':
        cd = pd.read_csv(DATADIR/'condom_use.csv') if condoms else None
        mods['networks'] = [sti.StructuredSexual(condom_data=cd)]
        people = EmodPeople(inp, risk_net='structuredsexual')
        if hiv:
            mods['diseases'] = EmodHIVOnSti(inp, record=record)
    elif network == 'emod':
        mods['networks'] = [EmodNetworkAt(inp, record=record)]
        mods['interventions'] = [EmodRiskAt(inp)]
        people = EmodPeople(inp, warm='emodnet')
        if hiv:
            mods['diseases'] = EmodHIV(inp, record=record)
    else:
        raise ValueError(f'network must be "sti" or "emod", not {network!r}')
    return sti.Sim(n_agents=n_agents, start=start, stop=stop, dt=1/12, rand_seed=seed, verbose=0,
                   demographics=dem['demographics'] + [people], people=dem['people'], total_pop=dem['total_pop'], **mods, **kwargs)
