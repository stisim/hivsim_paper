"""
Analysis 10: EMOD′ and Optima′ on Starsim 4.0.0rc3, against the Starsim ebd52c1a runs saved in analyses 04, 06 and 08. See README.md.

    python run.py emod                             # analysis 06's and 04's EMOD′ runs on rc3 (90 runs), to outputs/
    python run.py optima                           # analysis 08's three Optima′ runs on rc3, to outputs/
    python run.py arrays --tag rc3 [--repo PATH]   # Optima′ as acceptance.py runs it: stocks and flows by year, to outputs/optima_arrays_<tag>.npz
    python run.py time --tag rc3 [--repo PATH]     # time EMOD′ runs (set 0, seed 1: condensed, compiled), to outputs/time_<tag>.csv
    python run.py compare                          # the comparison tables, to outputs/
    python run.py --plot                           # figures/

`--repo` is the checkout whose models are run (default: this one); with the ebd52c1a environment, give one at the commit
before the migration (e.g. ~/hivsim_paper at 36a4fe7), since the models' code changed with the Starsim API.
"""
import sys
import json
import time
import argparse
import importlib.util
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
A04, A06, A08 = [repo / 'analyses' / n for n in ('04_emod_full', '06_condensed_cascade', '08_optima_condensed_cascade')]
SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
SEEDS = [1, 2, 3]
VARIANTS = ['compiled', 'condensed', 'interpreted'] # 06's two, and 04's (the campaign interpreted, cascade.py)
BASE = dict(compiled=(A06, 'compiled'), condensed=(A06, 'condensed'), interpreted=(A04, 'v4')) # The ebd52c1a runs: analysis folder and model name


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def setup(root):
    for p in ('common', 'analyses/04_emod_full', 'models/emod/reference', 'models/emod/v4'): # v4 first: common/cascade would hide its cascade.py
        sys.path.insert(0, str(Path(root) / p))
    return


# %% EMOD′

def emod_one(job):
    setup(repo)
    import harness # noqa: F401 (first, so that `cascade` is EMOD′'s: emod_targets puts common/ first on the path)
    pid, seed, variant = job
    if variant == 'interpreted':
        run04 = load_module(A04 / 'run.py', 'run04')
        t = time.time()
        _, _, values, series = run04.run_v4((pid, seed))
        return dict(tag=dict(model=variant, pid=pid, seed=seed), secs=time.time() - t, values=values, series=series.assign(model=variant, pid=pid, seed=seed), infants=None)
    run06 = load_module(A06 / 'run.py', 'run06')
    r = run06.run_one((pid, seed, variant))
    return dict(tag=r['tag'], secs=r['secs'], values=r['values'], series=r['series'], infants=r['infants'])


def emod(cores):
    jobs = [(p, s, v) for p in SETS for s in SEEDS for v in VARIANTS]
    with Pool(cores) as pool:
        res = pool.map(emod_one, jobs, chunksize=1)
    out = here / 'outputs'
    out.mkdir(exist_ok=True)
    for v in VARIANTS:
        rows = [r for r in res if r['tag']['model'] == v]
        (out / f'{v}_results.jsonl').write_text(''.join(json.dumps(dict(i=i, kind='sample', pars=dict(set=r['tag']['pid'], seed=r['tag']['seed']), values=r['values'])) + '\n' for i, r in enumerate(rows)))
    pd.concat([r['series'] for r in res], ignore_index=True).to_csv(out / 'series.csv.gz', index=False)
    pd.concat([r['infants'] for r in res if r['infants'] is not None], ignore_index=True).to_csv(out / 'infants.csv.gz', index=False)
    pd.DataFrame([dict(**r['tag'], secs=r['secs']) for r in res]).to_csv(out / 'runtime.csv', index=False)
    return


def timing(tag, root):
    """ Run time of one condensed and one compiled run (set 0, seed 1), on whichever Starsim is imported """
    setup(root)
    import starsim as ss, harness, inputs
    rows = []
    for variant in ('condensed', 'compiled'):
        sim = harness.make_sim(0, seed=inputs.load(0).run_number*10 + 1, cascade=True, compiled=variant == 'compiled', condensed=variant == 'condensed', record=False)
        t = time.time()
        sim.run()
        rows.append(dict(tag=tag, starsim=ss.__version__, variant=variant, secs=round(time.time() - t, 1)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(here / 'outputs' / f'time_{tag}.csv', index=False)
    return


# %% Optima′

def arrays(tag, root):
    """ Optima′'s whole model as acceptance.py runs it: people by state, population and year, new infections and HIV deaths by year """
    sys.path.insert(0, str(Path(root) / 'models' / 'optima' / 'v4'))
    import sciris as sc, starsim as ss
    import optima_prime as opp, compare as cmp
    refs = sc.load(cmp.REFS)
    years = refs['E']['original']['years']
    sim = opp.make_sim(components=cmp.COMPONENTS['E'], dt=refs['E']['dt'], end=int(years[-1])).run()
    np.savez_compressed(here / 'outputs' / f'optima_arrays_{tag}.npz', years=years, people=cmp.prime_people(sim, years), inci=cmp.prime_inci(sim, years), deaths=cmp.prime_deaths(sim, years), starsim=ss.__version__)
    return


def optima_one(variant):
    sys.path.insert(0, str(repo))
    return load_module(A08 / 'run.py', 'run08').optima_one(variant)


def optima(cores):
    with Pool(cores) as pool:
        res = pool.map(optima_one, ['own', 'emodcx', 'emodcx_noflows']) # Analysis 08's VARIANTS
    out = here / 'outputs'
    for r in res:
        (out / f'optima_{r["variant"]}_results.jsonl').write_text(json.dumps(dict(i=0, kind='sample', pars=dict(variant=r['variant']), values=r['values'])) + '\n')
    pd.concat([r['series'] for r in res]).round(5).to_csv(out / 'optima_series.csv', index=False)
    return


# %% Comparisons

def z_table(a, b, t):
    """ Target points in both sets of runs: medians, 5-95% of a, and z = (median b - median a)/SE (SE from both, as analysis 06) """
    rows = []
    for tid in t.index:
        if tid in a and tid in b:
            x, y = a[tid].dropna(), b[tid].dropna()
            if len(x) and len(y):
                rows.append(dict(tid=tid, target=t.loc[tid, 'target'], year=t.loc[tid, 'year'], a=x.median(), b=y.median(), a_lo=x.quantile(0.05), a_hi=x.quantile(0.95),
                                 se=np.sqrt(x.var()/len(x) + y.var()/len(y))))
    d = pd.DataFrame(rows)
    return d.assign(ratio=d.b/d.a, z=(d.b - d.a)/d.se)


def by_kind(d):
    return d.assign(within2=d.z.abs() <= 2).groupby('target').agg(n=('tid', 'size'), within2=('within2', 'sum'), z_mean=('z', 'mean'), z_max=('z', lambda z: z.abs().max()),
                                                                  ratio_median=('ratio', 'median'), ratio_lo=('ratio', lambda r: r.quantile(0.05)), ratio_hi=('ratio', lambda r: r.quantile(0.95))).round(3)


def compare():
    sys.path.insert(0, str(repo / 'common'))
    import targetfit
    t = targetfit.load_targets()
    out = here / 'outputs'
    runs = lambda folder, model: targetfit.load_runs(folder, model)[1]
    rc3 = {v: runs(here, v) for v in VARIANTS}

    # Each variant: rc3 against its ebd52c1a runs
    tabs, kinds, summary = [], [], []
    for v in VARIANTS:
        d = z_table(runs(*BASE[v]), rc3[v], t).assign(variant=v)
        tabs.append(d)
        kinds.append(by_kind(d).assign(variant=v))
        (pa, va, _), (pb, vb, _) = targetfit.load_runs(*BASE[v]), targetfit.load_runs(here, v) # Run by run (same set and seed): identical?
        key = lambda p: list(zip(p['set'], p['seed']))
        va.index, vb.index = key(pa), key(pb)
        same = (va.loc[vb.index, vb.columns].fillna(-1) == vb.fillna(-1)).all(axis=1)
        summary.append(dict(comparison=f'{v}: rc3 vs ebd52c1a', points=len(d), within2=int((d.z.abs() <= 2).sum()), within3=int((d.z.abs() <= 3).sum()), z_mean=round(d.z.mean(), 3), z_sd=round(d.z.std(), 3), z_max=round(d.z.abs().max(), 2),
                            runs_identical=f'{int(same.sum())}/{len(same)}'))
    pd.concat(tabs).round(5).to_csv(out / 'targets_vs_ebd52c1a.csv', index=False)
    pd.concat(kinds).to_csv(out / 'targets_vs_ebd52c1a_by_kind.csv')

    # Analysis 06 again: condensed against compiled, both on rc3
    d = z_table(rc3['compiled'], rc3['condensed'], t)
    by_kind(d).to_csv(out / 'condensed_vs_compiled_by_kind.csv')
    summary.append(dict(comparison='condensed vs compiled, rc3 (06)', points=len(d), within2=int((d.z.abs() <= 2).sum()), within3=int((d.z.abs() <= 3).sum()), z_mean=round(d.z.mean(), 3), z_sd=round(d.z.std(), 3), z_max=round(d.z.abs().max(), 2)))
    pd.DataFrame(summary).to_csv(out / 'z_summary.csv', index=False)

    # Analysis 04 again: EMOD′ (interpreted) / EMOD, median ratio by target, on each version
    emod = runs(A04, 'emod')
    r = pd.concat([by_kind(z_table(emod, runs(A04, 'v4'), t))[['n', 'ratio_median', 'ratio_lo', 'ratio_hi']].add_prefix('ebd52c1a_'),
                   by_kind(z_table(emod, rc3['interpreted'], t))[['ratio_median', 'ratio_lo', 'ratio_hi']].add_prefix('rc3_')], axis=1)
    r.to_csv(out / 'emod_ratio_by_kind.csv')

    # The cascade over time: 2020 medians, each variant on each version
    s = pd.read_csv(out / 'series.csv.gz')
    base = pd.concat([pd.read_csv(A06 / 'outputs' / 'series.csv.gz'), pd.read_csv(A04 / 'outputs' / 'series.csv.gz').query('model == "v4"').assign(model='interpreted')])
    cols = ['plhiv', 'on_art', 'Newly Infected', 'Died_from_HIV', 'diagnosed', 'art_OffART', 'art_OnEffectiveART']
    m = pd.concat([x.assign(version=ver) for ver, x in (('ebd52c1a', base), ('rc3', s))])
    m[m.year.isin([2005, 2010, 2015, 2020])].groupby(['model', 'year', 'version'])[cols].median().round(1).unstack('version').to_csv(out / 'series_medians.csv')

    # Run time: the 90 runs here (busy machine) and the paired timing
    rt = pd.read_csv(out / 'runtime.csv').groupby('model').secs.median().rename('rc3_median_secs').to_frame()
    rt['ebd52c1a_median_secs'] = pd.read_csv(A06 / 'outputs' / 'runtime.csv').groupby('model').secs.median()
    rt.round(1).to_csv(out / 'runtime_medians.csv')
    pd.concat([pd.read_csv(f) for f in sorted(out.glob('time_*.csv'))]).pivot(index='variant', columns='tag', values='secs').to_csv(out / 'time_paired.csv')
    for name in ('z_summary', 'targets_vs_ebd52c1a_by_kind', 'condensed_vs_compiled_by_kind', 'emod_ratio_by_kind', 'runtime_medians', 'time_paired'):
        print(f'\n{name}\n', pd.read_csv(out / f'{name}.csv').to_string(index=False))
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(repo / 'common'))
    from targetfit import style, c_run, c_best, c_muted
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)
    looks = dict(ebd52c1a=dict(color=c_best, label='Starsim ebd52c1a'), rc3=dict(color=c_run, label='Starsim 4.0.0rc3'))

    # The cascade over time, each variant on each version
    s = pd.read_csv(here / 'outputs' / 'series.csv.gz').assign(version='rc3')
    base = pd.concat([pd.read_csv(A06 / 'outputs' / 'series.csv.gz'), pd.read_csv(A04 / 'outputs' / 'series.csv.gz').query('model == "v4"').assign(model='interpreted')]).assign(version='ebd52c1a')
    s = pd.concat([base, s])
    s['art_cov'] = s.on_art/s.plhiv
    panels = [('plhiv', 'People with HIV (agents)'), ('Newly Infected', 'New infections a year'), ('Died_from_HIV', 'HIV deaths a year'), ('art_cov', 'On ART, share of PLHIV')]
    fig, axs = plt.subplots(len(VARIANTS), len(panels), figsize=(14, 9))
    for row, v in zip(axs, VARIANTS):
        for ax, (col, title) in zip(row, panels):
            for ver, lk in looks.items():
                g = s[(s.model == v) & (s.version == ver) & (s.year < 2025)].groupby('year')[col]
                ax.fill_between(g.median().index, g.quantile(0.05), g.quantile(0.95), color=lk['color'], alpha=0.15, lw=0)
                ax.plot(g.median().index, g.median(), lw=1.5, **lk)
            style(ax, f'{v.title()} campaign: {title[0].lower() + title[1:]}' if col == 'plhiv' else title)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.text(0.01, -0.01, 'Lines: median over 30 runs (10 sets × 3 seeds); bands: 5-95%. Rows: compiled campaign and condensed cascade (analysis 06), interpreted campaign (analysis 04). The versions agree run by run, so rc3 (blue) covers ebd52c1a (orange).', fontsize=7, color=c_muted, va='top')
    fig.savefig(figdir / 'cascade_over_time.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # z of each target point, rc3 against ebd52c1a
    d = pd.read_csv(here / 'outputs' / 'targets_vs_ebd52c1a.csv')
    fig, axs = plt.subplots(1, len(VARIANTS), figsize=(14, 3.5), sharey=True)
    for ax, v in zip(axs, VARIANTS):
        x = d[d.variant == v].reset_index(drop=True)
        ax.axhspan(-2, 2, color=c_muted, alpha=0.12, lw=0)
        ax.plot(x.index, x.z, 'o', color=c_run, ms=2.5)
        style(ax, f'{v.title()}: {int((x.z.abs() <= 2).sum())} of {len(x)} within ±2', 'z (rc3 − ebd52c1a median)/SE' if v == VARIANTS[0] else None)
        ax.set_xlabel('Target point (common/targets.csv order)', fontsize=7)
    fig.text(0.01, -0.03, 'Shaded: ±2 standard errors. SE from both sets of 30 runs; between-set spread makes it conservative, as in analysis 06.', fontsize=7, color=c_muted, va='top')
    fig.savefig(figdir / 'targets_z.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('part', nargs='?', choices=['emod', 'optima', 'arrays', 'time', 'compare'])
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=30)
    ap.add_argument('--tag', default='rc3')
    ap.add_argument('--repo', default=str(repo))
    a = ap.parse_args()
    (here / 'outputs').mkdir(exist_ok=True)
    if a.plot:
        plot()
    else:
        dict(emod=lambda: emod(a.cores), optima=lambda: optima(3), arrays=lambda: arrays(a.tag, a.repo), time=lambda: timing(a.tag, a.repo), compare=compare)[a.part]()
