"""
Optima's population label (docs/swap_design.md §1.3): each agent's Optima population, from its age, sex and sex work,
for hybrids where Optima′'s components sit on another model's people.

Optima's Malawi project (models/optima/reference/optima_malawi_ref.obj, `simpars['popkeys']`) has 39 populations:
- FSW (women), Clients (men) and MSM (men), all 15-49. In Optima they have no ages: nobody ages into or out of them;
  people move between them and the general populations of the same sex aged 15-49 by `risktransit` (FSW ↔ each of
  F15-19 ... F45-49 at a mean of 32 years each, so 4.6 years in sex work; Clients ↔ M15-19 ... M45-49 at 47 years each,
  6.7 years), into each of the seven bands alike, whatever the person's age. MSM have no transitions at all.
- F0-4, M0-4, ..., F80-84, M80-84, F85+, M85+: everyone else, by sex and 5-year age band.

`label()` is that definition read the other way: a woman aged 15-49 who sells sex now is FSW, a man aged 15-49 who
buys it is a client, a man aged 15-49 who has sex with men is MSM (first that applies), and everyone else is in their
sex and age band. Neither EMOD′ nor HIVsim has MSM, so that population is empty on their people.

- `stisim_label(sim)`: HIVsim's people (STIsim's `StructuredSexual`: `fsw` and `client` are whether each sells or
  buys sex now, from their lifetime fate and the ages it starts and stops). Its risk groups 0-2 have no counterpart
  in Optima's populations (Optima has no risk within the general population), so they are not used.
- `emod_label(sim)`: EMOD′'s people (`emod_age` in days; `risk` HIGH is a current FSW or client, by sex). EMOD's
  nodes and its LOW and MEDIUM risk have no counterpart in Optima's populations either.
- `PopLabel`: an analyzer that keeps each agent's label up to date (`optima_pop`, from the step's end), and records the
  number of agents in each population at each whole year (`sizes`, years × populations, in people).
"""
import numpy as np
import pandas as pd
import starsim as ss

POPS = ['FSW', 'Clients', 'MSM'] + [f'{s}{a}-{a+4}' for a in range(0, 85, 5) for s in 'FM'] + ['F85+', 'M85+']
FSW, CLIENTS, MSM = 0, 1, 2
KEY = slice(0, 3)
NBANDS = 18 # 0-4 ... 80-84, 85+


def definitions():
    """ The 39 populations: sex, ages and who is in them """
    rows = [dict(pop='FSW', sex='f', age_lo=15, age_hi=50, who='women who sell sex now'),
            dict(pop='Clients', sex='m', age_lo=15, age_hi=50, who='men who buy sex now'),
            dict(pop='MSM', sex='m', age_lo=15, age_hi=50, who='men who have sex with men')]
    for p in POPS[3:]:
        lo = int(p[1:].split('-')[0].rstrip('+'))
        rows.append(dict(pop=p, sex=p[0].lower(), age_lo=lo, age_hi=lo + 5 if lo < 85 else np.inf, who='everyone else'))
    return pd.DataFrame(rows)


def check(ref):
    """ That POPS are the project's populations, in its order """
    assert list(ref.simpars['popkeys']) == POPS, 'Populations differ from the Optima project'


def label(age, female, fsw=None, client=None, msm=None):
    """ Each person's Optima population (index into POPS), from age in years, sex and current sex work (bool arrays, or None for none) """
    age, female = np.asarray(age, float), np.asarray(female, bool)
    band = np.clip(np.floor(age/5).astype(int), 0, NBANDS - 1)
    out = 3 + 2*band + ~female
    adult = (age >= 15) & (age < 50)
    no = np.zeros(len(age), bool)
    as_bool = lambda x: no if x is None else np.asarray(x, bool)
    for idx, flag, sex in ((MSM, msm, ~female), (CLIENTS, client, ~female), (FSW, fsw, female)): # The first that applies last
        out = np.where(adult & sex & as_bool(flag), idx, out)
    return out


def stisim_label(sim, uids=None):
    """ HIVsim's (STIsim's) people """
    ppl = sim.people
    uids = ppl.auids if uids is None else uids
    net = sim.networks.structuredsexual
    return label(ppl.age.raw[uids], ppl.female.raw[uids], net.fsw.raw[uids], net.client.raw[uids])


def emod_label(sim, uids=None):
    """ EMOD′'s people: Risk HIGH is a current FSW (women) or client (men) """
    ppl = sim.people
    uids = ppl.auids if uids is None else uids
    high = ppl.risk.raw[uids] == 2
    female = ppl.female.raw[uids]
    return label(ppl.emod_age.raw[uids]/365, female, high & female, high & ~female)


class PopLabel(ss.Analyzer):
    """
    Each agent's Optima population, kept up to date each step (`optima_pop`), and the number in each at each whole year
    (`sizes`: a DataFrame, years × POPS, in people, i.e. agents × `scale`)

    Args:
        labeller (func): `labeller(sim, uids)`, e.g. stisim_label or emod_label
        scale (float): people per agent (default: the sim's `pop_scale`)
    """
    def __init__(self, labeller, scale=None, name='poplabel', **kwargs):
        super().__init__(name=name, **kwargs)
        self.labeller, self.scale = labeller, scale
        self.define_states(ss.IntArr('optima_pop', default=-1))
        self.rows = {}
        return

    def init_post(self):
        super().init_post()
        if self.scale is None:
            self.scale = float(self.sim.pop_scale or 1)
        return

    def step(self):
        au = self.sim.people.auids
        lab = self.labeller(self.sim, au)
        self.optima_pop[au] = lab
        year = self.t.now('year')
        y = int(np.round(year))
        if abs(year - y) < 0.5*self.t.dt_year - 1e-9 and y not in self.rows: # The step nearest each whole year
            self.rows[y] = np.bincount(lab, minlength=len(POPS))*self.scale
        return

    def finalize(self):
        super().finalize()
        self.sizes = pd.DataFrame.from_dict(self.rows, orient='index', columns=POPS).sort_index()
        self.sizes.index.name = 'year'
        return
