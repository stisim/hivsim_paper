"""
The policy outputs of docs/swap_design.md 2.2 (all but the intervention ranking), and the common targets
(common/targets.csv), from any agent-based model with a `CascadeHost` that has `cascade_stage()`.

`Outputs` is an analyzer (a `TransmissionLog`, common/transmission.py, so it also logs each sexual transmission with
the source's cascade stage). Each step it counts, by sex, 5-year age band and risk HIGH: people, infected,
diagnosed, on ART, suppressed (on the suppressive arm); and the step's new infections and HIV deaths, and the HIV
deaths of people on ART (at the end of the step before, as the log's ART history). After the run:

    out = sim.analyzers.outputs
    out.strata                 # the counts, one row per step and stratum
    out.df                     # the transmission log
    yearly(out)                # one row per year: the outputs of 2.2 and their numerators
    measure(out.strata, targets, scale)   # {tid: value} for common/targets.csv

What it reads from the disease (`sim.diseases.hiv`): EmodHIV's `sexual`, `new_uids` and `dead_uids`; for a Starsim
disease (STIsim's `sti.HIV`), its infection log (so the sim needs `ss.infection_log()`), `ti_infected` and `ti_dead`.
"""
import numpy as np
import pandas as pd
import starsim as ss
from common.transmission import TransmissionLog, classify, COARSE, COARSE_ORDER

BANDS = np.arange(0, 100, 5) # Lower edges of the age bands (the last is 95+)
COLS = ['pop', 'infected', 'diagnosed', 'on_art', 'suppressed', 'new_infections', 'hiv_deaths', 'hiv_deaths_on_art']
SEXUAL = ('structuredsexual',) # Routes in a Starsim infection log that count as sexual transmissions


def sexual_pairs(sim):
    """ This step's sexual transmissions: (recipients, sources) """
    hiv = sim.diseases.hiv
    if hasattr(hiv, 'sexual'): # EmodHIV
        return hiv.sexual
    if hiv.infection_log is None:
        raise ValueError('Outputs needs the sources of transmissions: add ss.infection_log() to the sim')
    d = hiv.infection_log.to_df(ti=sim.ti)
    if not len(d):
        return np.zeros(0, int), np.zeros(0, int)
    d = d[d.network.isin(SEXUAL)]
    return d.target.to_numpy(int), d.source.to_numpy(int)


def new_infections(sim):
    hiv = sim.diseases.hiv
    return hiv.new_uids if hasattr(hiv, 'new_uids') else (hiv.ti_infected == sim.ti).uids


def hiv_deaths(sim):
    hiv = sim.diseases.hiv
    return hiv.dead_uids if hasattr(hiv, 'dead_uids') else (hiv.ti_dead == sim.ti).uids


class Outputs(TransmissionLog):
    """
    Args:
        host (class): the model's `CascadeHost`, with `cascade_stage()`
        scale (float): people per agent, for counts in measure() (default: `sim.pop_scale`)
    """
    def __init__(self, host, scale=None, name='outputs', **kwargs):
        super().__init__(host, sexual_pairs, name=name, **kwargs)
        self.scale = scale
        self.counts, self.times = [], []
        return

    def init_post(self):
        super().init_post()
        if self.scale is None:
            self.scale = float(self.sim.pop_scale or 1)
        try:
            self.host.risk(ss.uids([0]))
            self.has_risk = True
        except NotImplementedError:
            self.has_risk = False
        return

    def key(self, uids):
        """ Stratum of each: (female, risk HIGH, age band) as one index """
        age = np.clip(np.searchsorted(BANDS, self.host.age_years(uids), side='right') - 1, 0, len(BANDS) - 1) # Unborn (infected in pregnancy) in the first
        high = self.host.risk(uids) == 2 if self.has_risk else np.zeros(len(uids), bool)
        return (self.host.female(uids).astype(int)*2 + high)*len(BANDS) + age

    def step(self):
        super().step() # The log, and self.on: who is on ART now (as before for those dying)
        sim, ppl = self.sim, self.sim.people
        au = ppl.auids
        au = au[(ppl.ti_dead.raw[au] != self.ti) & (self.host.age_years(au) >= 0)] # At the end of the step, and born
        inf = self.host.infected(au)
        st = self.host.cascade_stage(au)
        new, dead = ss.uids(new_infections(sim)), ss.uids(hiv_deaths(sim))
        n = 4*len(BANDS)
        c = lambda k, w=None: np.bincount(k, weights=w, minlength=n)
        k = self.key(au)
        out = [c(k), c(k, inf), c(k, inf & st['diagnosed']), c(k, st['on_art']), c(k, st['on_art'] & st['suppressive']),
               c(self.key(new)), c(self.key(dead)), c(self.key(dead), self.on.raw[dead])]
        self.counts.append(np.array(out).T)
        self.times.append(self.t.now('year'))
        return

    def finalize(self):
        super().finalize()
        nb = len(BANDS)
        idx = np.arange(4*nb)
        strat = pd.DataFrame(dict(female=idx//(2*nb) == 1, high=(idx//nb) % 2 == 1, age=BANDS[idx % nb]))
        n = len(self.times)
        self.strata = pd.concat([strat]*n, ignore_index=True)
        self.strata.insert(0, 'year', np.repeat(self.times, len(idx)))
        self.strata[COLS] = np.concatenate(self.counts) if n else np.zeros((0, len(COLS)))
        self.dt_year = float(self.t.dt_year)
        return


def yearly(out):
    """
    One row per calendar year: the outputs of docs/swap_design.md 2.2 and their numerators. Transmissions by the
    source's cascade stage (`tx_<coarse category>` shares of `tx_n`, sexual transmissions); deaths on ART
    (`share_deaths_on_art` of `hiv_deaths`); HIV incidence at 15-49 by sex per 1000 HIV-negative person-years
    (`incidence_15_49_f`, `_m`); and new infections at all ages (`new_infections`, agents, not scaled)
    """
    s = out.strata.assign(yr=np.floor(out.strata.year + 1e-9).astype(int))
    g = s.groupby('yr')
    y = pd.DataFrame(dict(new_infections=g.new_infections.sum(), hiv_deaths=g.hiv_deaths.sum(), hiv_deaths_on_art=g.hiv_deaths_on_art.sum()))
    y['share_deaths_on_art'] = y.hiv_deaths_on_art / y.hiv_deaths
    a = s[(s.age >= 15) & (s.age < 50)]
    for sex, f in (('f', True), ('m', False)):
        b = a[a.female == f].groupby('yr')
        py = (b['pop'].sum() - b.infected.sum())*out.dt_year
        y[f'new_15_49_{sex}'] = b.new_infections.sum()
        y[f'py_15_49_{sex}'] = py
        y[f'incidence_15_49_{sex}'] = 1000*y[f'new_15_49_{sex}']/py
    log = out.df
    if len(log):
        cat = classify(log).map(COARSE)
        tx = pd.crosstab(np.floor(log.year + 1e-9).astype(int).rename('yr'), pd.Categorical(cat, COARSE_ORDER), dropna=False)
        tx = tx.reindex(columns=COARSE_ORDER, fill_value=0)
        y['tx_n'] = tx.sum(1)
        for c in COARSE_ORDER:
            y[f'tx_{c}'] = tx[c] / y.tx_n
        y = y.fillna({'tx_n': 0})
    y.index.name = 'year'
    return y.reset_index()


def measure(strata, targets, scale=1, dt=None):
    """ {tid: value} for the rows of common/targets.csv: stocks at the step nearest each target's year; counts
    (`plhiv`, `on_art`, `population`, `new_infections`, `aids_deaths`) times `scale`; flows over the calendar year;
    `incidence` the new infections in the year to the survey per HIV-negative person at the survey (as
    models/emod/reference/emod_targets.py); group `fsw` is women of risk HIGH; NaN for years outside the run """
    s = strata
    times = np.unique(s.year)
    out = {}
    for tid, r in targets.iterrows():
        if not times[0] <= r.year <= times[-1] + 1e-6:
            out[tid] = np.nan # Outside the run
            continue
        m = (s.age >= r.age_lo) & (s.age < r.age_hi)
        if r.sex in ('f', 'm'):
            m &= s.female == (r.sex == 'f')
        if r.group == 'fsw':
            m &= s.female & s.high
        elif r.group != 'all':
            raise ValueError(r.group)
        if r.target in ('new_infections', 'aids_deaths'):
            col = 'new_infections' if r.target == 'new_infections' else 'hiv_deaths'
            out[tid] = s.loc[m & (np.floor(s.year + 1e-9) == np.floor(r.year)), col].sum()*scale
            continue
        at = s[m & (s.year == times[np.abs(times - r.year).argmin()])]
        pop, inf, dx, art, supp = (at[c].sum() for c in ('pop', 'infected', 'diagnosed', 'on_art', 'suppressed'))
        if r.target == 'incidence':
            new = s.loc[m & (s.year > r.year - 1) & (s.year <= r.year), 'new_infections'].sum()
            out[tid] = new/(pop - inf)
            continue
        out[tid] = dict(prevalence=lambda: inf/pop, plhiv=lambda: inf*scale, on_art=lambda: art*scale, population=lambda: pop*scale,
                        art_coverage=lambda: art/inf, diagnosed=lambda: dx/inf, on_ART_among_diagnosed=lambda: art/dx,
                        suppressed_among_on_ART=lambda: supp/art, vls_among_plhiv=lambda: supp/inf)[r.target]()
    return {int(k): float(v) for k, v in out.items()}
