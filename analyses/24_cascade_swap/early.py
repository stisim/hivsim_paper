"""
Analysis 24 addendum: why PLHIV and new infections are higher early with Optima's cascade on EMOD′-fixed (README,
"Addendum"). Run with Starsim and STIsim on PYTHONPATH, as run.py.

    python early.py track [CORES]   # EMOD′-fixed own, + Optima's table, + with its PMTCT per person; set 0, seeds 1-8, to 2012, with Track
    python early.py tables          # outputs/early_*.csv and outputs/early_report.txt, from run.py emod_pm's store and the tracked runs
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
import run # noqa: E402

OUT = HERE / 'outputs'
TRACK = OUT / 'track'
ARMS = {'emod_fixed': 'own', 'optima': 'optima', 'optima_pm': 'optima_pm', 'optima_swp_pm': 'optima_swp_pm'}


class Track(ss.Analyzer):
    """ EMOD′'s prevention and births, at the middle of each year (counts are agents) """
    def init_pre(self, sim):
        super().init_pre(sim)
        self.rows = []

    def step(self):
        sim, yr = self.sim, self.t.now('year')
        if abs((yr % 1) - 0.5) > 1/24 - 1e-9 or (yr % 1) < 0.5 - 1e-9:
            return
        p, h, dem = sim.people, sim.diseases.hiv, sim.demographics.emoddemog
        au = p.auids
        age = p.emod_age.raw[au]/365
        f = p.female.raw[au]
        inf = h.infected.raw[au]
        art = h.on_art(au)
        m1549 = ~f & (age >= 15) & (age < 50)
        preg = ~np.isnan(dem.preg_days.raw[au])
        kid = age < 15
        self.rows.append(dict(year=int(yr), circumcised_15_49=p.circumcised.raw[au][m1549].mean(), sti=p.sti.raw[au].mean(),
                              pregnant_hiv=(preg & inf).sum(), pregnant_hiv_dx=(preg & inf & h.diagnosed.raw[au]).sum(),
                              pregnant_hiv_art=(preg & inf & art).sum(), pmtct=(preg & inf & (h.pmtct_until.raw[au] >= yr) & (h.pmtct.raw[au] > 0)).sum(),
                              prep=(h.prep_until.raw[au] >= yr).sum(), plhiv_kids=(kid & inf).sum(), art_kids=(kid & art).sum(),
                              diagnosed=(inf & h.diagnosed.raw[au]).sum(), on_art=(inf & art).sum()))

    def finalize(self):
        super().finalize()
        self.df = pd.DataFrame(self.rows)


def _track(job):
    from common import swaps
    arm, seed = job
    f = TRACK / f'{ARMS[arm]}_{seed}.csv'
    if f.exists():
        return arm, seed, 'stored'
    c = run.emod_pm_config()
    c['stop'] = 2012
    cell = next(x for x in swaps.cells(c) if x['pair'] == 0 and x['seed'] == seed)
    ch = next(h for h in c['hybrids'] if h['C'] == arm)
    sim = swaps.build(ch, cell, c)
    sim.pars.analyzers = sc.tolist(sim.pars.analyzers) + [Track(name='track')]
    sim.run()
    from common.outputs import yearly
    y = yearly(sim.analyzers.outputs)[['year', 'new_infections', 'tx_n', 'hiv_deaths', 'plhiv']]
    sim.analyzers.track.df.merge(y, on='year').to_csv(f, index=False)
    return arm, seed, 'ran'


def track(cores=24):
    from multiprocessing import Pool
    TRACK.mkdir(parents=True, exist_ok=True)
    jobs = [(a, s) for s in range(1, 9) for a in ('emod_fixed', 'optima', 'optima_pm')]
    with Pool(cores) as pool:
        for r in pool.imap_unordered(_track, jobs):
            print(*r, flush=True)


# %% Tables

def stored():
    """ run.py emod_pm's yearly tables, per arm, indexed by cell and year (agents) """
    from common import swaps
    df = swaps.load(run.emod_pm_config())
    df['arm'] = df.name.str.split('C=').str[1].map(ARMS)
    df['mtct'] = df.new_infections - df.tx_n
    return {a: d.set_index(['cell', 'year']).sort_index() for a, d in df.groupby('arm')}


def paired(d, arm, cols, years, scale=100):
    """ Mean paired Δ (arm − own) per year, people """
    return ((d[arm][cols] - d['own'][cols]).groupby('year').mean()*scale).loc[years]


def tables():
    import tables as T
    d = stored()
    cols = ['new_infections', 'tx_n', 'mtct', 'hiv_deaths', 'plhiv', 'diagnosed', 'on_art']
    lines = []
    yrs = list(range(1985, 2011))
    own = (d['own'][cols].groupby('year').mean()*100).loc[yrs]
    by_year = pd.concat({'own': own} | {f'Δ {a}': paired(d, a, cols, yrs) for a in d if a != 'own'}, axis=1).round(0)
    by_year.to_csv(OUT / 'early_by_year.csv')
    lines += ['Year by year, 1985-2010, own and paired Δ (people; 48 cells)', by_year.to_string(), '']
    # Cumulative 1990-2020 and by period
    rows = []
    for a in [a for a in d if a != 'own']:
        for lo, hi in [(1990, 2020), (1990, 2000), (2001, 2004), (2005, 2010), (2011, 2015), (2016, 2020)]:
            s = lambda z: z.loc[(slice(None), slice(lo, hi)), ['new_infections', 'tx_n', 'mtct', 'hiv_deaths']].groupby('cell').sum()
            x = (s(d[a]) - s(d['own']))*100
            x['plhiv_end'] = (d[a].xs(hi, level='year').plhiv - d['own'].xs(hi, level='year').plhiv)*100
            for q in x:
                m, l, h, n = T.ci(x[q])
                rows.append(dict(arm=a, period=f'{lo}-{hi}', q=q, mean=m, lo=l, hi=h))
    cum = pd.DataFrame(rows)
    cum.to_csv(OUT / 'early_cumulative.csv', index=False)
    lines += ['Cumulative paired Δ (people; 95% CI): infections, sexual (tx_n), mother-to-child, HIV deaths; PLHIV at the end of the period']
    for (a, p), g in cum.groupby(['arm', 'period'], sort=False):
        lines.append(f'{a:14s} {p}: ' + '; '.join(f'{r.q} {r.mean/1e3:+.1f}k ({r.lo/1e3:+.1f}, {r.hi/1e3:+.1f})' for r in g.itertuples()))
    lines.append('')
    # Headline (2020, 2024) as tables.main's direction 1
    rows = []
    for a, z in d.items():
        for cell, g in z.groupby('cell'):
            q = T.from_yearly(g.reset_index('cell', drop=True).reset_index(), 100)
            for yr in T.YEARS:
                rows.append(dict(arm=a, cell=cell, year=yr, **q.loc[yr]))
    em = pd.DataFrame(rows)
    em.to_csv(OUT / 'early_cells.csv', index=False)
    gap = pd.read_csv(OUT / 'gap.csv') if (OUT / 'gap.csv').exists() else None
    lines.append('Direction 1, paired Δ against own (95% CI)')
    for yr in T.YEARS:
        a0 = em[(em.year == yr) & (em.arm == 'own')].set_index('cell')
        lines.append(f'{yr}: own ' + ', '.join(f'{q} {T.level(a0[q].mean(), q)}' for q in T.Q))
        for a in [a for a in d if a != 'own']:
            b = em[(em.year == yr) & (em.arm == a)].set_index('cell').reindex(a0.index)
            lines.append(f'  {a:14s} ' + '; '.join(f'{q} {T.fmt(*T.ci(b[q] - a0[q])[:3], q)}' for q in T.Q))
    lines.append('')
    # Tracked runs
    if TRACK.exists() and any(TRACK.glob('*.csv')):
        tr = pd.concat([pd.read_csv(f).assign(arm=f.stem.rsplit('_', 1)[0], seed=int(f.stem.rsplit('_', 1)[1])) for f in sorted(TRACK.glob('*.csv'))])
        t = tr.groupby(['arm', 'year']).mean(numeric_only=True).drop(columns='seed')
        cnt = [c for c in t if c not in ('circumcised_15_49', 'sti')]
        t[cnt] = t[cnt]*100
        t = t.loc[(slice(None), slice(1995, 2012)), :]
        t.to_csv(OUT / 'early_track.csv')
        lines += ['Tracked runs (set 0, seeds 1-8, mid-year; counts in people): VMMC and STI flag as shares', t.round(3).to_string()]
    (OUT / 'early_report.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    if sys.argv[1] == 'track':
        track(int(sys.argv[2]) if len(sys.argv) > 2 else 24)
    elif sys.argv[1] == 'tables':
        tables()
