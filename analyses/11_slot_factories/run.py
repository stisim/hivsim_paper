"""
Analysis 11: HIVsim on rc3, EMOD′'s 1985 burn-in saved and resumed, and HIVsim at 10,000 agents against 1%. See README.md.

    python run.py hivsim --tag rc3                                 # HIVsim (own care), 10,000 agents, seeds 1-20, to outputs/hivsim_rc3.csv
    python run.py hivsim --tag before --repo ~/hivsim_paper        # the same, in env hivsim_paper (before the migration)
    python run.py scale                                            # HIVsim at 1% scale, seeds 1-10, to outputs/hivsim_scale1pct.csv
    python run.py burnin                                           # EMOD′ continuous vs burnt in to 1985, saved, loaded and run on (set 0, seeds 1-3)
    python run.py compare                                          # the comparison tables, to outputs/
    python run.py --plot                                           # figures/

Run with env hp_rc3 and PYTHONPATH=<STIsim rs-v4-art>, except `hivsim --tag before`.
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
import time
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
out = here / 'outputs'
figs = here / 'figures'
raw = repo / 'raw_results' / '11_slot_factories'
SEEDS = range(1, 21)
SEEDS_1PCT = range(1, 11)
N_1PCT = 76_205 # 1% of Malawi's 1985 population (malawi/data/malawi_age_1985.csv: 7,620,508)
MEASURES = ['alive', 'hiv.prevalence_15_49', 'hiv.new_infections', 'hiv.new_deaths', 'hiv.on_art', 'hiv.p_on_art', 'hiv.diagnosed', 'hiv.infected']
RATES = ['hiv.prevalence_15_49', 'hiv.p_on_art', 'incidence', 'death_per_infected'] # Independent of the number of agents (derived in yearly())
YEARS = range(2000, 2021)


# %% HIVsim

def hivsim_one(job):
    root, seed, n_agents = job
    sys.path.insert(0, str(Path(root) / 'models' / 'hivsim'))
    import swap
    t = time.time()
    sim = swap.make_sim('own', seed=seed, n_agents=n_agents, stop=2022)
    sim.run()
    secs = time.time() - t
    df = sim.to_df(resample='year', use_years=True, sep='.')
    df = df[['time'] + MEASURES].rename(columns={'time': 'year'}).assign(seed=seed, n_agents=n_agents, secs=secs)
    return df


def run_hivsim(root, seeds, n_agents, fname):
    with Pool(min(len(seeds), 40)) as pool:
        dfs = pool.map(hivsim_one, [(root, s, n_agents) for s in seeds], chunksize=1)
    df = pd.concat(dfs, ignore_index=True)
    out.mkdir(exist_ok=True)
    df.to_csv(out / fname, index=False)
    print(f'{fname}: {len(dfs)} runs, median {df.groupby("seed").secs.first().median():.1f} s')
    return df


# %% EMOD′ burn-in

def emod_setup():
    for p in ('common', 'models/emod/reference', 'models/emod/v4'):
        sys.path.insert(0, str(repo / p))
    sys.path.insert(0, str(repo))
    return


def emod_sim(seed, stop):
    emod_setup()
    import harness
    return harness.make_sim(pid=0, seed=seed, stop=stop, cascade=True, condensed=True)


def burnin_one(job):
    """ One seed: continuous; then burnt in to 1985, saved, loaded and run on in the same process. Saves the burn-in and the continuous run's arrays for resume_one() """
    seed, stop = job
    emod_setup()
    import starsim as ss
    import sciris as sc
    from common.identity import arrays, differences
    times = dict(seed=seed, stop=stop)

    t = time.time()
    cont = arrays(emod_sim(seed, stop).run())
    times['continuous'] = time.time() - t
    raw.mkdir(parents=True, exist_ok=True)
    sc.save(raw / f'emod_continuous_set0_seed{seed}.arrays', cont)

    t = time.time()
    sim = emod_sim(seed, stop)
    sim.init()
    years = np.asarray(sim.t.yearvec, dtype=float)
    last = float(years[np.argmin(np.abs(years - 1985)) - 1]) # The step that ends at 1985 is the last of the burn-in
    sim.run(until=last)
    times['burnin'] = time.time() - t
    times['burnin_last_step'] = last
    fn = raw / f'emod_burnin_set0_seed{seed}.sim'
    t = time.time()
    sim.save(fn)
    times['save'] = time.time() - t
    times['file_mb'] = fn.stat().st_size / 1e6
    del sim
    t = time.time()
    sim = ss.load(fn)
    times['load'] = time.time() - t
    t = time.time()
    sim.run()
    times['resume'] = time.time() - t
    d = differences(cont, arrays(sim)).assign(seed=seed, process='same')
    times['n_arrays'] = len(cont)
    times['n_differ_same_process'] = len(d)
    return times, d


def resume_one(seed):
    """ In a new process (as every hybrid will be): load the burn-in, run it on, and compare with the continuous run """
    emod_setup()
    import starsim as ss
    import sciris as sc
    import harness # noqa: F401 (so that the modules the pickle names are importable as they were)
    from common.identity import arrays, differences
    t = time.time()
    sim = ss.load(raw / f'emod_burnin_set0_seed{seed}.sim')
    load = time.time() - t
    t = time.time()
    sim.run()
    d = differences(sc.load(raw / f'emod_continuous_set0_seed{seed}.arrays'), arrays(sim)).assign(seed=seed, process='new')
    return dict(seed=seed, load_new_process=load, resume_new_process=time.time() - t, n_differ_new_process=len(d)), d


def run_burnin(seeds=(1, 2, 3), stop=2025):
    import multiprocessing as mp
    with Pool(len(seeds)) as pool:
        res = pool.map(burnin_one, [(s, stop) for s in seeds], chunksize=1)
    with mp.get_context('spawn').Pool(len(seeds)) as pool: # Fresh processes, which have run nothing before
        res2 = pool.map(resume_one, seeds, chunksize=1)
    times = pd.DataFrame([r[0] for r in res]).merge(pd.DataFrame([r[0] for r in res2]), on='seed')
    diffs = pd.concat([r[1] for r in res + res2], ignore_index=True)
    out.mkdir(exist_ok=True)
    times.to_csv(out / 'burnin_times.csv', index=False)
    diffs.to_csv(out / 'burnin_differences.csv', index=False)
    print(times.round(2).to_string())
    print(f'{len(diffs)} array differences' + (':\n' + diffs.groupby(['process', 'seed']).head(10).to_string() if len(diffs) else ''))
    return times, diffs


# %% Comparisons

def yearly(df):
    """ Measures by seed and year, with rates that don't depend on the number of agents """
    df = df.assign(year=df.year.astype(float).round().astype(int))
    return df.assign(incidence=df['hiv.new_infections'] / (df['alive'] - df['hiv.infected']), death_per_infected=df['hiv.new_deaths'] / df['hiv.infected'])


def z_table(a, b, measures):
    """ Per measure and year: medians, 5-95% of a, and z = (median b - median a)/SE (SE from both sets, as analysis 06) """
    rows = []
    for m in measures:
        for y in YEARS:
            x, w = a.loc[a.year == y, m].dropna(), b.loc[b.year == y, m].dropna()
            rows.append(dict(measure=m, year=y, a=x.median(), b=w.median(), a_lo=x.quantile(0.05), a_hi=x.quantile(0.95), se=np.sqrt(x.var()/len(x) + w.var()/len(w))))
    d = pd.DataFrame(rows)
    return d.assign(ratio=d.b / d.a, z=(d.b - d.a) / d.se)


def by_measure(d):
    return d.assign(within2=d.z.abs() <= 2).groupby('measure', sort=False).agg(n=('year', 'size'), within2=('within2', 'sum'), z_mean=('z', 'mean'), z_max=('z', lambda z: z.abs().max()),
                                                                               ratio_median=('ratio', 'median'), ratio_lo=('ratio', 'min'), ratio_hi=('ratio', 'max')).round(3)


def compare():
    rc3 = yearly(pd.read_csv(out / 'hivsim_rc3.csv'))
    before = yearly(pd.read_csv(out / 'hivsim_before.csv'))
    measures = MEASURES[:-1] + ['incidence', 'death_per_infected']
    z = z_table(before, rc3, measures)
    z.round(5).to_csv(out / 'hivsim_rc3_vs_before.csv', index=False)
    s = by_measure(z)
    s.to_csv(out / 'hivsim_rc3_vs_before_by_measure.csv')
    print('rc3 (b) vs before (a):\n' + s.to_string())

    pct = out / 'hivsim_scale1pct.csv'
    if pct.exists():
        one = yearly(pd.read_csv(pct))
        z = z_table(rc3, one, measures)
        z.round(5).to_csv(out / 'hivsim_1pct_vs_10k.csv', index=False)
        s = by_measure(z)
        s.to_csv(out / 'hivsim_1pct_vs_10k_by_measure.csv')
        print('1% (b) vs 10,000 agents (a):\n' + s.to_string())
    rt = pd.concat([d.groupby('seed').agg(secs=('secs', 'first'), n_agents=('n_agents', 'first')).assign(run=k) for k, d in
                    (('before', before), ('rc3', rc3)) + ((('1pct', one),) if pct.exists() else ())])
    rt.groupby('run').secs.describe().round(1).to_csv(out / 'hivsim_runtime.csv')
    print(rt.groupby('run').secs.describe().round(1))
    return


def plot():
    import matplotlib.pyplot as plt
    figs.mkdir(exist_ok=True)
    sets = [('before (rc1)', yearly(pd.read_csv(out / 'hivsim_before.csv')), '#888888'), ('rc3, 10,000 agents', yearly(pd.read_csv(out / 'hivsim_rc3.csv')), '#1f5fa8')]
    if (out / 'hivsim_scale1pct.csv').exists():
        sets.append(('rc3, 1% (76,205 agents)', yearly(pd.read_csv(out / 'hivsim_scale1pct.csv')), '#d9822b'))
    panels = [('hiv.prevalence_15_49', 'Prevalence 15-49'), ('incidence', 'Incidence (per susceptible per year)'), ('hiv.new_infections', 'New infections'),
              ('hiv.on_art', 'On ART'), ('hiv.p_on_art', 'Share of infected on ART'), ('hiv.new_deaths', 'HIV deaths'), ('alive', 'Population')]
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    for ax, (m, title) in zip(axes.flat, panels):
        for label, df, c in sets:
            q = df[df.year.between(1986, 2021)].groupby('year')[m].quantile([0.05, 0.5, 0.95]).unstack()
            ax.fill_between(q.index, q[0.05], q[0.95], color=c, alpha=0.2, lw=0)
            ax.plot(q.index, q[0.5], color=c, label=label)
        ax.set_title(title, fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)
    axes.flat[-1].axis('off')
    axes.flat[0].legend(fontsize=8, frameon=False)
    fig.suptitle('HIVsim Malawi (experiment 03 best set): median and 5-95% over seeds', fontsize=11)
    fig.tight_layout()
    fig.savefig(figs / 'hivsim_outputs.png', dpi=120)
    print(f'Wrote {figs / "hivsim_outputs.png"}')
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('what', nargs='?', choices=['hivsim', 'scale', 'burnin', 'compare'])
    ap.add_argument('--tag', default='rc3')
    ap.add_argument('--repo', default=str(repo))
    ap.add_argument('--stop', type=float, default=2025)
    ap.add_argument('--plot', action='store_true')
    a = ap.parse_args()
    if a.what == 'hivsim':
        run_hivsim(Path(a.repo).expanduser(), list(SEEDS), 10e3, f'hivsim_{a.tag}.csv')
    elif a.what == 'scale':
        run_hivsim(repo, list(SEEDS_1PCT), N_1PCT, 'hivsim_scale1pct.csv')
    elif a.what == 'burnin':
        run_burnin(stop=a.stop)
    elif a.what == 'compare':
        compare()
    if a.plot:
        plot()
