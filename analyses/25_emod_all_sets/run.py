"""
Analysis 25: the 129 distinct calibrated EMOD sets with both PR #4 fixes. See README.md.

    python run.py          # per run: totals, Figure 3A shares, fit to Daniel's calibration targets -> outputs/runs.csv
    python run.py --plot   # outputs/sets.csv, summary tables, figures/total_vs_fit.{png,svg}
"""
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))
sys.path.insert(0, str(here))
from common import transmission as tx # noqa: E402
import fit # noqa: E402

RUNS = repo / 'raw_results' / 'emod' / 'a25_both'
PAPER = dict(undiagnosed='undiagnosed', diagnosed_naive='diagnosed_never_treated', art_lt6m='on_art', art_supp='on_art',
             art_nonsupp='on_art', interrupted='interrupted_art')
CATS = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted_art']
CENSUS, CENSUS_YEAR = 17_563_749, 2018.42 # Daniel's scaling (mwi_transmission.Rmd), as analysis 18
PUBLISHED = 30_582
out = here / 'outputs'
OBS = fit.load_obs()


def one(folder):
    pid, seed = int(folder.name[3:6]), int(folder.name.split('_s')[1])
    o = folder / 'output'
    rep = pd.read_csv(o / 'ReportHIVByAgeAndGender.csv')
    rep.columns = rep.columns.str.strip()
    row = dict(pid=pid, seed=seed, **{f'fit_{k}': v for k, v in fit.score(rep, OBS).items()})
    factor = CENSUS / rep.loc[np.isclose(rep.Year, CENSUS_YEAR), 'Population'].sum()
    adult = rep[(rep.Age >= 15) & (rep.Age < 65)]
    yr = adult.Year.astype(int)
    ni, nhalf = adult.groupby(yr)['Newly Infected'].sum(), adult.groupby(yr).Year.nunique()
    sh = tx.shares(tx.from_emod(o), coarse=True)
    sh.pop('n')
    sh = sh.T.groupby(PAPER).sum().T.reindex(columns=CATS)
    for y in (2020, 2024):
        row[f'ni_{y}'] = ni[y] * factor * 2 / nhalf[y] # annualised over the half-years present (2024 misses the last)
        row.update({f'{c}_{y}': 100 * sh.loc[y, c] for c in CATS})
    row['factor'] = factor
    return row


def reduce(cores):
    folders = sorted(f for f in RUNS.glob('set*_s*') if (f / 'output' / 'InsetChart.json').exists()) # Written at the end of a run
    with Pool(cores) as pool:
        res = pd.DataFrame(pool.map(one, folders))
    out.mkdir(exist_ok=True)
    res.to_csv(out / 'runs.csv', index=False)
    print(f'{len(res)} runs, {res.pid.nunique()} sets')


def wq(x, w, q):
    """Weighted quantile (weights = multiplicity among the 250)"""
    i = np.argsort(x.to_numpy())
    x, cw = x.to_numpy()[i], np.cumsum(w.to_numpy()[i]) / w.sum()
    return np.interp(q, cw - w.to_numpy()[i] / w.sum() / 2, x)


def tables():
    runs = pd.read_csv(out / 'runs.csv')
    smap = pd.read_csv(out / 'set_map.csv')
    mult = smap.groupby('rep_pid').agg(n_copies=('parameterization_id', 'size'), daniel_lik=('likelihood', 'first'))
    sets = runs.groupby('pid').mean(numeric_only=True).drop(columns='seed').join(runs.groupby('pid').size().rename('n_seeds')).join(mult)
    sets = sets.sort_values('fit_total', ascending=False)
    sets['fit_rank'] = np.arange(1, len(sets) + 1)
    sets['daniel_rank'] = sets.daniel_lik.rank(ascending=False, method='min').astype(int)
    sets.round(4).to_csv(out / 'sets.csv')
    w = sets.n_copies
    lines = [f'{len(runs)} runs, {len(sets)} distinct sets ({int(w.sum())} of the 250), seeds per set: {sets.n_seeds.value_counts().to_dict()}']
    for y in (2020, 2024):
        x = sets[f'ni_{y}']
        lines.append(f'{y} new infections 15-64, all 250 (weighted by copies): median {wq(x, w, .5):,.0f}, mean {np.average(x, weights=w):,.0f}, '
                     f'5-95% {wq(x, w, .05):,.0f}-{wq(x, w, .95):,.0f}; distinct sets unweighted median {x.median():,.0f}; '
                     f'per-run median {runs[f"ni_{y}"].median():,.0f} (5-95% {runs[f"ni_{y}"].quantile(.05):,.0f}-{runs[f"ni_{y}"].quantile(.95):,.0f})')
    lines.append(f'Share of sets (weighted) with 2024 total <= {PUBLISHED:,}: {w[sets.ni_2024 <= PUBLISHED].sum() / w.sum():.2f}')
    from scipy.stats import spearmanr
    for k in ['fit_total', 'fit_Prevalence', 'fit_ARTCoverage', 'fit_OnART', 'fit_Population', 'daniel_lik']:
        r, p = spearmanr(sets[k], sets.ni_2024)
        lines.append(f'Spearman(2024 total, {k}) = {r:.2f} (p = {p:.2g})')
    r, p = spearmanr(sets.fit_total, sets.daniel_lik)
    lines.append(f'Spearman(fixed-model fit, Daniel likelihood) = {r:.2f} (p = {p:.2g})')
    fig = pd.read_csv(repo / 'common' / 'fig3a_malawi.csv', comment='#')
    fig = fig[fig.model == 'EMOD'].set_index(['year', 'category']).share
    groups = {'all 250 (weighted)': (sets, w)}
    for k in (10, 25):
        groups[f'top {k} by fixed-model fit'] = (sets.nsmallest(k, 'fit_rank'), None)
        groups[f'top {k} by Daniel likelihood (distinct)'] = (sets.nsmallest(k, 'daniel_rank'), None)
    groups['10 most-copied sets'] = (sets.nlargest(10, 'n_copies'), None)
    groups['sets with 2024 total within 10% of 30,582'] = (sets[(sets.ni_2024 / PUBLISHED - 1).abs() < 0.1], None)
    rows = []
    for g, (d, ww) in groups.items():
        ww = ww if ww is not None else pd.Series(1, index=d.index)
        r = dict(group=g, n_sets=len(d), ni_2024_mean=np.average(d.ni_2024, weights=ww), ni_2024_median=wq(d.ni_2024, ww, .5),
                 ni_2020_mean=np.average(d.ni_2020, weights=ww))
        for y in (2020, 2024):
            for c in CATS:
                r[f'{c}_{y}'] = np.average(d[f'{c}_{y}'], weights=ww)
        r['max_abs_share_diff'] = max(abs(r[f'{c}_{y}'] - fig[(y, c)]) for y in (2020, 2024) for c in CATS)
        rows.append(r)
    rows.append(dict(group='Figure 3A / caption', ni_2024_mean=PUBLISHED, **{f'{c}_{y}': fig[(y, c)] for y in (2020, 2024) for c in CATS}))
    gt = pd.DataFrame(rows).round(1)
    gt.to_csv(out / 'groups.csv', index=False)
    (out / 'summary.txt').write_text('\n'.join(lines) + '\n')
    return sets, lines, gt


def plot(sets):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': '#888888', 'xtick.color': '#555555', 'ytick.color': '#555555'})
    f, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.axhline(PUBLISHED, color='#3a3a3a', lw=1.5, zorder=1)
    ax.text(sets.fit_total.min(), PUBLISHED * 0.97, 'Published EMOD, 2024: 30,582 (Figure 3 caption)', color='#3a3a3a', va='top', fontsize=8)
    ax.scatter(sets.fit_total, sets.ni_2024 / 1, s=12 + 6 * sets.n_copies, color='#2a78d6', alpha=0.75, edgecolor='white', linewidth=0.8, zorder=2)
    ax.set_xlabel('Fit to EMOD calibration targets, fixed model (approx. log-likelihood, higher = better)')
    ax.set_ylabel('New infections 15–64, 2024 (mean of 2 seeds)')
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v / 1000:.0f}k'))
    ax.set_ylim(0, None)
    ax.grid(axis='y', color='#e6e6e6', lw=0.6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title('EMOD, both PR #4 fixes: 129 distinct calibrated sets\n(dot size = copies among the 250)', loc='left', fontsize=9)
    f.tight_layout()
    (here / 'figures').mkdir(exist_ok=True)
    f.savefig(here / 'figures' / 'total_vs_fit.png', dpi=200)
    f.savefig(here / 'figures' / 'total_vs_fit.svg')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=16)
    a = ap.parse_args()
    if a.plot:
        sets, lines, gt = tables()
        print('\n'.join(lines))
        pd.set_option('display.width', 250)
        print(gt.to_string(index=False))
        plot(sets)
    else:
        reduce(a.cores)
