"""
Optima′'s side of the condensed cascade and of the outputs (common/cascade/interface.py, common/outputs.py), as agents:
`OptimaHost` implements `CascadeHost` on Optima′'s declared states, moving each agent with `ss.move(uids=)`, and
`OptimaTracker` (an intervention, after Optima's own rules: `agents.make_sim(interventions=[OptimaTracker()])`) keeps
what the states don't say: who has ever been on ART, who has had first symptoms, the step's events, and the sources of
the step's sexual infections.

How the interface maps onto Optima′ (each agent has a population, `people.pop`, and one of Optima's 38 states):
- Age and sex come from the population (poplabel.POPS): sex as the population's; age fixed within its band for each
  person (FSW, Clients and MSM: 15-49), so it jumps between bands when Optima ages people. Optima has no ages.
- Risk: FSW and Clients are HIGH (as EMOD's), everyone else MEDIUM; Optima has no risk within the general population.
- CD4 is a value in the band of the person's CD4 state, fixed for each person within each band (acute: as > 500).
  WHO stage from it (1: ≥ 500, 2: 350-500, 3: 200-350, 4: < 200), and first symptoms on first reaching CD4 < 200
  (Optima's AIDS: the gt50 and lt50 bands). Optima has no nadir: `record_cd4()` keeps the lowest one recorded.
- A test diagnoses the undiagnosed (undx → dx, same CD4); ART starts move from any state with HIV not on ART into svl
  (arm 'suppressive') or usvl (otherwise); stops into lost; circumcision moves uncircumcised susceptibles (susreg) into
  progcirc (men with HIV have no circumcision state in Optima: they are left as they are).
- Not representable per agent in Optima′ (raise if asked for anyone): pregnancy (Optima's births are a rate per woman,
  with no mother), so `pregnant()` is False and `twelve_weeks_pregnant()` and `six_weeks_old()` are empty; PMTCT and
  PrEP (in Optima, coverages by population inside births and the force of infection).
- Places in care (common/cascade/pathways.py) are Optima's states: `care_place()`, `cd4_band()` and `population()`
  read them, and `rate_scale()` is the survival of Optima's 0.2-year step, among whom Optima moves people; `link_to_care` moves dx and lost into care, `lose_from_care` care into lost, `switch_arm` between usvl
  and svl, each in its CD4 band.
- `cascade_stage()`: `diagnosed` (not undx), `on_art` (usvl or svl), `suppressive` (svl), `ever_treated`, and
  `category`, Figure 3A's (docs/core_question.md): undiagnosed; diagnosed, never treated (dx, and care and lost
  never on ART); on ART; interrupted (care and lost who have been on ART), analysis 19's mapping D8.
- Sources of sexual infections: Optima's force of infection on each population is a sum over sources (partner population
  × infected state, `Transmission.parts`), with no pairs. Each infection's source is drawn from those parts for the
  recipient's population and circumcision status in that step, as analysis 19 attributes Optima's infections, and then
  an agent in that state and population (at the end of the step) is drawn as the source. The draws use their own
  random numbers, so they don't change the run.
"""
import numpy as np
import starsim as ss
import optima_prime as opp
import poplabel as pl
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from common.cascade.interface import CascadeHost # noqa: E402
from common.cascade.pathways import PLACES # noqa: E402

CD4_RANGE = dict(acute=(500, 1000), gt500=(500, 1000), gt350=(350, 500), gt200=(200, 350), gt50=(50, 200), lt50=(0, 50))
INFECTED = [f'{c}_{k}' for c in opp.UNTREATED + opp.TREATED for k in opp.CD4]
DEFS = pl.definitions()
AGE_LO = DEFS.age_lo.to_numpy(float)
AGE_W = np.where(np.isfinite(DEFS.age_hi), DEFS.age_hi - DEFS.age_lo, 15).astype(float) # 85+: to 100
FEMALE = (DEFS.sex == 'f').to_numpy()
CHILD = np.array([p[1:] == '10-14' for p in pl.POPS])


def fixed_fraction(uids, k):
    """ A number in [0, 1) fixed for each agent (and each k), without random numbers """
    return (np.asarray(uids, float)*(0.6180339887498949 + 0.1*k)) % 1


class OptimaTracker(ss.Intervention):
    """ What the host needs and the states don't keep: see the module docstring """
    def __init__(self, name='optima_host', seed=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.define_states(ss.BoolArr('ever_treated', default=False), ss.BoolArr('symptomatic', default=False),
                           ss.FloatArr('lowest_cd4', default=np.inf), ss.IntArr('last_pop', default=-1))
        self.seed = seed
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        dis = sim.diseases.optima
        self.dis = dis
        self.rng = np.random.default_rng((sim.pars.rand_seed if self.seed is None else self.seed) + 2203)
        self.names = list(dis.foi.names) if dis.foi is not None else []
        self.hivnames = [n for n, h in zip(self.names, dis.foi.hiv) if h] if dis.foi is not None else INFECTED
        self._sexual, self._births, self._deaths = [], [], []
        for s in ('susreg', 'progcirc'):
            if f'infection_{s}' in dis.graph.transitions:
                dis.add_callback(f'infection_{s}', lambda u, c=int(s == 'progcirc'): self._sexual.append((ss.uids(u), c)))
        for name in dis.graph.transitions:
            if name.startswith('birth_'):
                dis.add_callback(name, lambda u: self._births.append(ss.uids(u)))
            elif name.startswith('hivdeath_'):
                dis.add_callback(name, lambda u: self._deaths.append(ss.uids(u)))
        self.pairs = (ss.uids(), ss.uids())
        self.new_sexual, self.new_births, self.dead = ss.uids(), ss.uids(), ss.uids()
        self.debuts = self.newly_symptomatic = ss.uids()
        self._done = None
        return

    def init_post(self):
        super().init_post()
        self.ever_treated[self.on_art(self.sim.people.auids)] = True # As analysis 19's tracer: those on ART at the start, and nobody in care or lost
        self.last_pop[self.sim.people.auids] = self.sim.people.pop.raw[self.sim.people.auids]
        return

    def on_art(self, uids):
        return np.any([self.dis[f'{c}_{k}'].raw[uids] for c in opp.TREATED for k in opp.CD4], axis=0)

    def step(self):
        """ After Optima's rules: everyone on ART now has been (only these rules and the host start ART); the step's debuts and first symptoms """
        au = self.sim.people.auids
        self.ever_treated[au[self.on_art(au)]] = True
        pop = self.sim.people.pop.raw[au]
        last = self.last_pop.raw[au]
        self.debuts = au[CHILD[np.maximum(last, 0)] & (last >= 0) & ~CHILD[pop] & (AGE_LO[pop] >= 15)]
        self.last_pop[au] = pop
        aids = np.any([self.dis[f'{c}_{k}'].raw[au] for c in opp.UNTREATED + opp.TREATED for k in opp.CD4[opp.AIDS:]], axis=0)
        self.newly_symptomatic = au[aids & ~self.symptomatic.raw[au]]
        self.symptomatic[self.newly_symptomatic] = True
        return

    def collect(self):
        """ The events of the step just finished (from the disease's callbacks), and the sources of its sexual infections; once a step """
        if self._done == self.ti:
            return
        self._done = self.ti
        cat = lambda lst: ss.uids(np.concatenate(lst)) if lst else ss.uids()
        sexual = self._sexual
        self.new_sexual = cat([u for u, _ in sexual])
        self.new_births, self.dead = cat(self._births), cat(self._deaths)
        self._sexual, self._births, self._deaths = [], [], []
        self.pairs = self.draw_sources(sexual)
        return

    def draw_sources(self, sexual):
        """ (recipients, sources) for the step's sexual infections: a source (partner population × infected state) from Optima's parts of the force of infection, then an agent in it """
        foi = self.dis.foi
        if not sexual or foi is None:
            return ss.uids(), ss.uids()
        ppl = self.sim.people
        au = ppl.auids
        state = np.full(len(au), -1)
        for k, n in enumerate(self.hivnames):
            state[self.dis[n].raw[au]] = k
        has = state >= 0
        cand, cstate, cpop = au[has], state[has], ppl.pop.raw[au[has]]
        npop, nst = len(pl.POPS), len(self.hivnames)
        key = cpop*nst + cstate
        order = np.argsort(key, kind='stable')
        cand, key = cand[order], key[order]
        start = np.searchsorted(key, np.arange(npop*nst))
        count = np.searchsorted(key, np.arange(npop*nst), side='right') - start
        parts = foi.parts # Acquiring pop × partner pop × circ × infected state, this step's
        dst, src = [], []
        for uids, c in sexual:
            pops = ppl.pop.raw[uids]
            for i in np.unique(pops):
                who = uids[pops == i]
                w = parts[i, :, c, :]*(count.reshape(npop, nst) > 0) # Only sources with agents in them
                if w.sum() <= 0:
                    w = (count.reshape(npop, nst) > 0).astype(float)
                cell = self.rng.choice(w.size, size=len(who), p=(w/w.sum()).ravel())
                pick = start[cell] + (self.rng.random(len(who))*count[cell]).astype(int)
                dst.append(who); src.append(cand[pick])
        return ss.uids(np.concatenate(dst)), ss.uids(np.concatenate(src))


class OptimaHost(CascadeHost):
    def __init__(self, sim):
        super().__init__(sim)
        self.dis = sim.diseases.optima
        self.ppl = sim.people
        self.tr = sim.interventions.optima_host
        self.sources = [] # (year, cascade states, categories) of each step's sources
        return

    def pop(self, uids):
        return self.ppl.pop.raw[uids]

    def state_in(self, uids, cascade, cd4=opp.CD4):
        return np.any([self.dis[f'{c}_{k}'].raw[uids] for c in cascade for k in cd4], axis=0)

    # People
    def age_years(self, uids):
        p = self.pop(uids)
        return AGE_LO[p] + fixed_fraction(uids, 0)*AGE_W[p]

    def female(self, uids):
        return FEMALE[self.pop(uids)]

    def pregnant(self, uids):
        return np.zeros(len(uids), bool)

    def post_debut(self, uids):
        return AGE_LO[self.pop(uids)] >= 15

    def dying(self, uids):
        return self.ppl.ti_dead.raw[uids] <= self.sim.ti

    def risk(self, uids):
        return np.where(self.pop(uids) <= pl.CLIENTS, 2, 1)

    # Events
    def debuts(self):
        return self.tr.debuts

    def newly_symptomatic(self):
        return self.tr.newly_symptomatic

    def twelve_weeks_pregnant(self):
        return ss.uids()

    def six_weeks_old(self):
        return ss.uids()

    # HIV
    def infected(self, uids):
        return self.state_in(uids, opp.UNTREATED + opp.TREATED)

    def move(self, uids, src, dst):
        """ Move these agents from cascade state `src` to `dst`, each in its CD4 band; returns those moved """
        moved = []
        for k in opp.CD4:
            u = uids[self.dis[f'{src}_{k}'].raw[uids]]
            if len(u):
                ss.move(self.dis[f'{src}_{k}'], self.dis[f'{dst}_{k}'], uids=u)
                moved.append(u)
        opp_stale()
        return ss.uids(np.concatenate(moved)) if moved else ss.uids()

    def test(self, uids):
        pos = self.infected(uids)
        self.move(uids, 'undx', 'dx')
        return pos

    def cd4(self, uids):
        out = np.full(len(uids), np.inf)
        f = fixed_fraction(uids, 1)
        for k, (lo, hi) in CD4_RANGE.items():
            m = self.state_in(uids, opp.UNTREATED + opp.TREATED, [k])
            out[m] = lo + f[m]*(hi - lo)
        return out

    def record_cd4(self, uids):
        lc = self.tr.lowest_cd4
        lc[uids] = np.minimum(lc.raw[uids], self.cd4(uids))
        return lc.raw[uids]

    def who_stage(self, uids):
        c = self.cd4(uids)
        return np.select([c >= 500, c >= 350, c >= 200], [1, 2, 3], 4).astype(float)

    def on_art(self, uids):
        return self.state_in(uids, opp.TREATED)

    def start_art(self, uids, arm):
        dst = 'svl' if arm == 'suppressive' else 'usvl'
        for c in opp.UNTREATED:
            self.tr.ever_treated[self.move(uids, c, dst)] = True
        return

    def stop_art(self, uids):
        for c in opp.TREATED:
            self.move(uids, c, 'lost')
        return

    def start_pmtct(self, uids, efficacy, days):
        if len(uids):
            raise NotImplementedError('Optima′ has no PMTCT per mother (its births have no mothers)')

    def circumcise(self, uids):
        men = uids[~self.female(uids) & self.dis.susreg.raw[uids]]
        if len(men):
            ss.move(self.dis.susreg, self.dis.progcirc, uids=men)
            opp_stale()
        return men

    def start_prep(self, uids, efficacy, days):
        if len(uids):
            raise NotImplementedError('Optima′ has no PrEP per person (a coverage by population in its force of infection)')

    # Places in care (common/cascade/pathways.py): Optima′'s own states
    def care_place(self, uids):
        out = np.zeros(len(uids), int)
        for i, c in enumerate(PLACES[1:], 1):
            out[self.state_in(uids, [c])] = i
        return out

    def cd4_band(self, uids):
        out = np.zeros(len(uids), int)
        for i, k in enumerate(opp.CD4):
            out[self.state_in(uids, opp.UNTREATED + opp.TREATED, [k])] = i
        return out

    def population(self, uids):
        return np.array(pl.POPS)[self.pop(uids)]

    def rate_scale(self, uids):
        """ Optima′ takes its cascade's moves among those who survive Optima's step of 0.2 years (HIV death by CD4 band, and on ART, then background death by population and year) """
        if not hasattr(self, '_death'):
            sp = opp.load().simpars
            self._tvec, self._bg = sp['tvec'], np.asarray(sp['death'], float)*opp.REF_DT
            d = np.array([sp[f'death{k}'] for k in opp.CD4])*opp.REF_DT
            self._death = dict(untreated=d, usvl=d*sp['deathusvl'], svl=d*sp['deathsvl'])
        band, place = self.cd4_band(uids), self.care_place(uids)
        d = np.where(place == PLACES.index('svl'), self._death['svl'][band], np.where(place == PLACES.index('usvl'), self._death['usvl'][band], self._death['untreated'][band]))
        bg = np.array([np.interp(self.sim.t.now('year'), self._tvec, row) for row in self._bg])[self.pop(uids)]
        return (1 - d)*(1 - bg)

    def link_to_care(self, uids):
        for c in ('dx', 'lost'):
            self.move(uids, c, 'care')

    def lose_from_care(self, uids):
        self.move(uids, 'care', 'lost')

    def switch_arm(self, uids, arm):
        src, dst = ('usvl', 'svl') if arm == 'suppressive' else ('svl', 'usvl')
        self.move(uids, src, dst)

    # Measurement
    def cascade_stage(self, uids):
        on = self.on_art(uids)
        diagnosed = self.infected(uids) & ~self.state_in(uids, ['undx'])
        ever = self.tr.ever_treated.raw[uids] | on
        category = np.select([~diagnosed, on, ever], ['undiagnosed', 'on_art', 'interrupted'], 'diagnosed_never_treated')
        category = np.where(self.infected(uids), category, '')
        return dict(diagnosed=diagnosed, on_art=on, suppressive=self.state_in(uids, ['svl']), ever_treated=ever, category=category)

    # For common/outputs.py: the step's sexual transmissions, new infections and HIV deaths
    def sexual_pairs(self):
        self.tr.collect()
        dst, src = self.tr.pairs
        if len(src): # For checks: each source's cascade state and Figure 3A category at the time
            st = np.full(len(src), '', dtype=object)
            for c in opp.UNTREATED + opp.TREATED:
                st[self.state_in(src, [c])] = c
            self.sources.append((self.sim.t.now('year'), st, self.cascade_stage(src)['category']))
        return dst, src

    def new_infections(self):
        self.tr.collect()
        return ss.uids(np.concatenate([self.tr.new_sexual, self.tr.new_births]))

    def hiv_deaths(self):
        self.tr.collect()
        return self.tr.dead


def opp_stale():
    """ The counts Optima′'s rules read are out of date after a move """
    import agents
    agents.stale()
