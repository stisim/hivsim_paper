"""
Analysis 05: EMOD′ with its campaign compiled into Starsim deliveries against the literal interpreter. See README.md.
Needs Starsim with PR #1538 on PYTHONPATH.

    python run.py          # 10 sets × 3 seeds × (literal, compiled), to outputs/
    python run.py --plot   # figures/
"""
import sys
import json
import time
import argparse
from collections import Counter
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
for p in ('common', 'analyses/04_emod_full', 'models/emod/reference', 'models/emod/v4'): # v4 first: common/cascade would hide its cascade.py
    sys.path.insert(0, str(repo / p))
import targetfit # noqa: E402

SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
SEEDS = [1, 2, 3]
VARIANTS = ['literal', 'compiled']


def run_one(job):
    import harness, inputs, emod_targets
    import importlib.util
    spec = importlib.util.spec_from_file_location('run04', repo / 'analyses' / '04_emod_full' / 'run.py')
    run04 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run04)
    series = run04.series # Analysis 04's yearly series
    from pmtct import tables # And its InterventionStatus and infant tables
    pid, seed, variant = job
    sim = harness.make_sim(pid, seed=inputs.load(pid).run_number*10 + seed, cascade=True, compiled=variant == 'compiled', record=False)
    t = time.time()
    sim.run()
    secs = time.time() - t
    d = sim.analyzers.emodreport.df.assign(sex=lambda x: x.Gender.map({0: 'm', 1: 'f'}))
    status, infants = tables(d)
    first = d[np.isclose(d.Year, d.Year.min())].drop(columns='sex')
    tag = dict(model=variant, pid=pid, seed=seed)
    return dict(tag=tag, secs=secs, values=emod_targets.measure(d, targetfit.load_targets()), series=series(d).assign(**tag),
                status=status.assign(**tag), infants=infants.assign(**tag), first=first.assign(**tag))


def reduce(cores):
    with Pool(cores) as pool:
        res = pool.map(run_one, [(p, s, v) for p in SETS for s in SEEDS for v in VARIANTS], chunksize=1)
    out = here / 'outputs'
    out.mkdir(exist_ok=True)
    for v in VARIANTS:
        rows = [r for r in res if r['tag']['model'] == v]
        (out / f'{v}_results.jsonl').write_text(''.join(json.dumps(dict(i=i, kind='sample', pars=dict(set=r['tag']['pid'], seed=r['tag']['seed']), values=r['values'])) + '\n' for i, r in enumerate(rows)))
    for k in ('series', 'status', 'infants', 'first'):
        pd.concat([r[k] for r in res], ignore_index=True).to_csv(out / f'{k}.csv.gz', index=False)
    pd.DataFrame([dict(**r['tag'], secs=r['secs']) for r in res]).to_csv(out / 'runtime.csv', index=False)
    import inputs
    from compile_campaign import compile_campaign
    holder, ds = compile_campaign(inputs.load(0))
    kinds = Counter(type(d.product).__name__ if d.product is not None else ('deliver=' if d.deliver is not None else 'none') for d in ds)
    delays = sum(d.delay is not None for d in ds)
    pd.DataFrame([dict(what=k, n=n) for k, n in kinds.items()] + [dict(what='deliveries', n=len(ds)), dict(what='with delay=', n=delays)]).to_csv(out / 'objects.csv', index=False)
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_muted
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)
    looks = dict(literal=dict(color=c_best, label='Literal (cascade.py)'), compiled=dict(color=c_run, label='Compiled (deliveries)'))
    note = 'Lines: median over 30 runs (10 sets × 3 seeds); bands: 5-95%.'

    def save(fig, name, text=note):
        fig.text(0.01, -0.03, text, fontsize=7, color=c_muted, va='top')
        fig.savefig(figdir / name, dpi=150, bbox_inches='tight')
        plt.close(fig)

    s = pd.read_csv(here / 'outputs' / 'series.csv.gz')
    s['off_share'] = s.art_OffART/s.plhiv
    s['art_cov'] = s.on_art/s.plhiv
    s['dx_share'] = s.diagnosed/s.plhiv
    s['supp_share'] = s.art_OnEffectiveART/s.on_art.where(s.on_art > 0)
    panels = [('plhiv', 'People with HIV (agents)'), ('on_art', 'On ART (agents)'), ('Newly Infected', 'New infections a year'),
              ('Died_from_HIV', 'HIV deaths a year'), ('dx_share', 'Diagnosed, share of PLHIV'), ('art_cov', 'On ART, share of PLHIV'),
              ('off_share', 'OffART, share of PLHIV'), ('supp_share', 'OnEffectiveART, share of those on ART')]
    fig, axs = plt.subplots(2, 4, figsize=(14, 6))
    for ax, (col, title) in zip(axs.flat, panels):
        for m, lk in looks.items():
            g = s[(s.model == m) & (s.year < 2025)].groupby('year')[col]
            ax.fill_between(g.median().index, g.quantile(0.05), g.quantile(0.95), color=lk['color'], alpha=0.15, lw=0)
            ax.plot(g.median().index, g.median(), lw=1.5, **lk)
        style(ax, title)
    axs[0, 0].legend(fontsize=7, frameon=False)
    save(fig, 'cascade_over_time.png')

    t = targetfit.load_targets()
    vals = {m: targetfit.load_runs(here, m)[1] for m in looks}
    rows = []
    for tid in t.index:
        if tid in vals['literal'] and tid in vals['compiled']:
            a, b = vals['literal'][tid].dropna(), vals['compiled'][tid].dropna()
            rows.append(dict(tid=tid, target=t.loc[tid, 'target'], literal=a.median(), compiled=b.median(), lo=a.quantile(0.05), hi=a.quantile(0.95),
                             se=np.sqrt(a.var()/len(a) + b.var()/len(b))))
    cmp = pd.DataFrame(rows).assign(ratio=lambda d: d.compiled/d.literal, z=lambda d: (d.compiled - d.literal)/d.se)
    cmp.to_csv(here / 'outputs' / 'targets_compared.csv', index=False)
    names = cmp.target.unique()
    fig, axs = plt.subplots(1, len(names), figsize=(2.2*len(names), 3.5), sharey=True)
    for ax, name in zip(axs, names):
        d = cmp[cmp.target == name]
        x = np.arange(len(d))
        ax.axhline(1, color=c_muted, lw=0.6)
        ax.vlines(x, d.lo/d.literal, d.hi/d.literal, color=c_best, alpha=0.4, lw=3)
        ax.plot(x, d.ratio, 'o', color=c_run, ms=3)
        style(ax, name.replace('_', ' '), 'Compiled / literal median' if name == names[0] else None)
        ax.set_xticks([])
        ax.set_ylim(0, 2)
    save(fig, 'targets_compared.png', 'Each target point: compiled median / literal median (blue), and the literal runs\' 5-95% over their median (orange).')

    st = pd.read_csv(here / 'outputs' / 'status.csv.gz', keep_default_na=False)
    pos = st[st.HasHIV == 1]
    tot = pos.groupby(['model', 'pid', 'seed', 'year']).n.sum().rename('tot').reset_index()
    sh = pos.groupby(['model', 'pid', 'seed', 'year', 'InterventionStatus']).n.sum().reset_index().merge(tot, on=['model', 'pid', 'seed', 'year'])
    sh['s'] = sh.n/sh.tot
    states = sorted(x for x in sh.InterventionStatus.unique() if sh[sh.InterventionStatus == x].s.max() > 0.005)
    fig, axs = plt.subplots(3, 4, figsize=(14, 8))
    for ax, x in zip(axs.flat, states):
        for m, lk in looks.items():
            g = sh[(sh.model == m) & (sh.InterventionStatus == x)].groupby('year').s.mean()
            ax.plot(g.index, g, lw=1.5, **lk)
        style(ax, x, 'Share of PLHIV')
    for ax in axs.flat[len(states):]:
        ax.set_visible(False)
    axs.flat[0].legend(fontsize=7, frameon=False)
    save(fig, 'intervention_status.png', 'Mean over 30 runs.')

    inf = pd.read_csv(here / 'outputs' / 'infants.csv.gz')
    fig, ax = plt.subplots(figsize=(5, 3.2))
    for m, lk in looks.items():
        g = inf[inf.model == m].groupby('year').new_0_4.mean()
        ax.plot(g.index, g, lw=1.5, **lk)
    style(ax, 'New infections at ages 0-4 a year (mean of runs)')
    ax.legend(fontsize=7, frameon=False)
    save(fig, 'infants.png', 'Mean over 30 runs.')
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=30)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
