"""
Analysis 04: EMOD′ with the whole campaign against the published EMOD (raw_results/emod/sub10). See README.md.

    python run.py          # run EMOD′ (10 sets × 3 seeds); measure it and EMOD on common/targets.csv; series to outputs/
    python run.py --plot   # figures/: the shared target figures for each model, and the comparisons
"""
import sys
import json
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
for p in ('models/emod/v4', 'models/emod/reference', 'common'):
    sys.path.insert(0, str(repo / p))
import targetfit # noqa: E402
import emod_targets # noqa: E402

SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
SEEDS = [1, 2, 3]
EMOD = repo / 'raw_results' / 'emod' / 'sub10'


def series(d):
    """ Yearly stocks (at x.92) and flows of a report in EMOD's format """
    d = d[d.Year >= 1990]
    yr = np.floor(d.Year).astype(int)
    flows = d.groupby(yr)[['Newly Infected', 'Died_from_HIV']].sum()
    s = d[np.isclose(d.Year % 1, 0.92)]
    s = s.assign(year=np.floor(s.Year).astype(int))
    out = s.groupby('year').agg(pop=('Population', 'sum'), plhiv=('Infected', 'sum'), on_art=('On_ART', 'sum'), diagnosed=('Diagnosed', 'sum'))
    pos = s[s.HasHIV == 1]
    for v in ['ARTNaive', 'OffART', 'OnEffectiveART', 'OnNonSuppressiveART']:
        out[f'art_{v}'] = pos[pos.ARTStatus == v].groupby('year').Infected.sum()
    if 'InterventionStatus' in s:
        for v, g in pos.groupby('InterventionStatus'):
            out[f'is_{v}'] = g.groupby('year').Infected.sum()
    return out.join(flows).fillna(0).reset_index()


def run_v4(job):
    import harness, inputs
    pid, seed = job
    sim = harness.make_sim(pid, seed=inputs.load(pid).run_number*10 + seed, cascade=True, record=False) # Some sets repeat others' parameters; EMOD tells them apart by Run_Number
    sim.run()
    d = sim.analyzers.emodreport.df.assign(sex=lambda x: x.Gender.map({0: 'm', 1: 'f'}))
    return pid, seed, emod_targets.measure(d, targetfit.load_targets()), series(d)


def run_emod(pid):
    d = emod_targets.load_report(EMOD / f'set{pid:03d}' / 'output' / 'ReportHIVByAgeAndGender.csv')
    return pid, 0, emod_targets.measure(d, targetfit.load_targets()), series(d)


def reduce(cores):
    with Pool(cores) as pool:
        res = {'v4': pool.map(run_v4, [(p, s) for p in SETS for s in SEEDS], chunksize=1), 'emod': pool.map(run_emod, SETS)}
    (here / 'outputs').mkdir(exist_ok=True)
    psets = pd.read_csv(emod_targets.PSETS).set_index('parameterization_id')
    ser = []
    for model, rows in res.items():
        lines = [json.dumps(dict(i=i, kind='sample', pars=dict(set=pid, seed=seed, likelihood=float(psets.loc[pid, 'likelihood'])), values=vals)) + '\n'
                 for i, (pid, seed, vals, _) in enumerate(rows)]
        (here / 'outputs' / f'{model}_results.jsonl').write_text(''.join(lines))
        ser += [s.assign(model=model, pid=pid, seed=seed) for pid, seed, _, s in rows]
    pd.concat(ser).to_csv(here / 'outputs' / 'series.csv.gz', index=False)
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_muted
    for model in ('emod', 'v4'):
        targetfit.plot(here, model)
    figdir = here / 'figures'
    looks = dict(emod=dict(color=c_best, label='EMOD, 10 sets'), v4=dict(color=c_run, label='EMOD′, 10 sets × 3 seeds'))
    note = 'Lines: median over runs; bands: 5-95% of runs (EMOD: 10 parameter sets, one run each; EMOD′: the same sets, 3 seeds each).'

    def save(fig, name, text=note):
        fig.text(0.01, -0.03, text, fontsize=7, color=c_muted, va='top')
        fig.savefig(figdir / name, dpi=150, bbox_inches='tight')
        plt.close(fig)

    s = pd.read_csv(here / 'outputs' / 'series.csv.gz')
    s['off_share'] = s.art_OffART/s.plhiv
    s['art_cov'] = s.on_art/s.plhiv
    s['dx_share'] = s.diagnosed/s.plhiv
    s['supp_share'] = s.art_OnEffectiveART/s.on_art.where(s.on_art > 0)
    panels = [('plhiv', 'People with HIV (agents)'), ('on_art', 'On ART (agents)'), ('Newly Infected', 'New infections a year (agents)'),
              ('Died_from_HIV', 'HIV deaths a year (agents)'), ('dx_share', 'Diagnosed, share of PLHIV'), ('art_cov', 'On ART, share of PLHIV'),
              ('off_share', 'OffART (dropped out), share of PLHIV'), ('supp_share', 'ARTStatus OnEffectiveART, share of those on ART')]
    fig, axs = plt.subplots(2, 4, figsize=(14, 6))
    for ax, (col, title) in zip(axs.flat, panels):
        for m, lk in looks.items():
            g = s[(s.model == m) & (s.year < 2025)].groupby('year')[col]
            ax.fill_between(g.median().index, g.quantile(0.05), g.quantile(0.95), color=lk['color'], alpha=0.15, lw=0)
            ax.plot(g.median().index, g.median(), lw=1.5, **lk)
        style(ax, title)
    axs[0, 0].legend(fontsize=7, frameon=False)
    save(fig, 'cascade_over_time.png')

    # Targets: EMOD′ against EMOD, each target
    t = targetfit.load_targets()
    vals = {m: targetfit.load_runs(here, m)[1] for m in looks}
    rows = []
    for tid in t.index:
        if tid in vals['emod'] and tid in vals['v4']:
            e, v = vals['emod'][tid].dropna(), vals['v4'][tid].dropna()
            rows.append(dict(tid=tid, target=t.loc[tid, 'target'], sex=t.loc[tid, 'sex'], group=t.loc[tid, 'group'], year=t.loc[tid, 'year'], age=f'{t.loc[tid, "age_lo"]}-{t.loc[tid, "age_hi"]}',
                             data=t.loc[tid, 'value'], emod=e.median(), emod_lo=e.quantile(0.05), emod_hi=e.quantile(0.95), v4=v.median(), v4_lo=v.quantile(0.05), v4_hi=v.quantile(0.95)))
    cmp = pd.DataFrame(rows)
    cmp['ratio'] = cmp.v4/cmp.emod
    cmp.to_csv(here / 'outputs' / 'targets_compared.csv', index=False)
    names = cmp.target.unique()
    fig, axs = plt.subplots(1, len(names), figsize=(2.2*len(names), 4), sharey=True)
    for ax, name in zip(np.atleast_1d(axs), names):
        d = cmp[cmp.target == name]
        x = np.arange(len(d))
        ax.axhline(1, color=c_muted, lw=0.6)
        ax.vlines(x, d.v4_lo/d.emod, d.v4_hi/d.emod, color=c_run, alpha=0.4, lw=3)
        ax.vlines(x + 0.25, d.emod_lo/d.emod, d.emod_hi/d.emod, color=c_best, alpha=0.4, lw=3)
        ax.plot(x, d.ratio, 'o', color=c_run, ms=3)
        style(ax, name.replace('_', ' '), 'EMOD′ / EMOD median' if name == names[0] else None)
        ax.set_xticks([])
        ax.set_ylim(0, 2)
    save(fig, 'targets_compared.png', 'For each target (points left to right, as in common/targets.csv): EMOD′ median / EMOD median (blue point), EMOD′ 5-95% (blue bar) and EMOD 5-95% (orange bar), both over EMOD\'s median.')

    # EMOD′ InterventionStatus
    v = s[s.model == 'v4']
    cols = [c for c in s.columns if c.startswith('is_')]
    g = v.groupby('year')[cols + ['plhiv']].median()
    fig, ax = plt.subplots(figsize=(8, 4))
    for c in cols:
        if g[c].max() > 0.01*g.plhiv.max():
            ax.plot(g.index, g[c]/g.plhiv, lw=1.3, label=c[3:])
    style(ax, 'EMOD′: InterventionStatus of people with HIV', 'Share of PLHIV')
    ax.legend(fontsize=7, frameon=False, ncol=2)
    save(fig, 'v4_intervention_status.png', 'Median over 30 runs. EMOD\'s report here has no InterventionStatus.')
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=16)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
