"""
Analysis 18: our EMOD runs (the PR #1 files, and with PR #4's fixes) against Figure 3A. See README.md.
Needs the EMOD runs: raw_results/emod/tr10x5 (analysis 09) and fig3a_{mort,drop,both} (README.md).

    python run.py          # classify each run, to outputs/
    python run.py --plot   # figures/ and the tables
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
from common import transmission as tx # noqa: E402

SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
SEEDS = [1, 2, 3, 4, 5]
EMOD = repo / 'raw_results' / 'emod'
VARIANTS = {'bug': 'tr10x5', 'mort': 'fig3a_mort', 'drop': 'fig3a_drop', 'both': 'fig3a_both'}
LABELS = {'bug': 'PR #1 files (as merged)', 'mort': 'mortality fix', 'drop': 'dropout fix', 'both': 'both fixes'}
PAPER = dict(undiagnosed='undiagnosed', diagnosed_naive='diagnosed_never_treated', art_lt6m='on_art', art_supp='on_art',
             art_nonsupp='on_art', interrupted='interrupted_art')
CATS = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted_art']
CENSUS, CENSUS_YEAR = 17_563_749, 2018.42 # Daniel's scaling (mwi_transmission.Rmd: 2018.5; our report steps are at .42 and .92)
SCALE = 0.01
YEARS = [2010, 2015, 2020, 2024]
out = here / 'outputs'


def one(job):
    v, pid, seed = job
    folder = EMOD / VARIANTS[v] / f'set{pid:03d}_s{seed}' / 'output'
    log = tx.from_emod(folder)
    sh = tx.shares(log, coarse=True)
    n = sh.pop('n')
    sh = sh.T.groupby(PAPER).sum().T.reindex(columns=CATS)
    # Steps present per calendar year (the runs stop at 2025.0, so the last step(s) of 2024 are missing)
    tr = pd.read_csv(folder / 'TransmissionReport.csv', usecols=['SIM_TIME', 'YEAR'])
    yr = np.floor(tr.YEAR).astype(int)
    nsteps = tr.groupby(yr).SIM_TIME.nunique()
    rep = pd.read_csv(folder / 'ReportHIVByAgeAndGender.csv')
    rep.columns = rep.columns.str.strip()
    pop = rep.loc[np.isclose(rep.Year, CENSUS_YEAR), 'Population'].sum()
    adult = rep[(rep.Age >= 15) & (rep.Age < 65)]
    ni = adult.groupby(adult.Year.astype(int))['Newly Infected'].sum()
    nhalf = adult.groupby(adult.Year.astype(int)).Year.nunique()
    tot = pd.DataFrame(dict(sexual=n, steps=nsteps, ni_15_64=ni, halves=nhalf)).assign(pop2018=pop)
    tag = dict(variant=v, pid=pid, seed=seed)
    return sh.reset_index().assign(**tag), tot.rename_axis('year').reset_index().assign(**tag)


def reduce(cores):
    jobs = [(v, p, s) for v in VARIANTS for p in SETS for s in SEEDS if (EMOD / VARIANTS[v] / f'set{p:03d}_s{s}' / 'output' / 'InsetChart.json').exists()] # Written at the end of a run
    with Pool(cores) as pool:
        res = pool.map(one, jobs)
    out.mkdir(exist_ok=True)
    pd.concat([r[0] for r in res]).to_csv(out / 'shares.csv.gz', index=False)
    pd.concat([r[1] for r in res]).to_csv(out / 'totals.csv.gz', index=False)
    print(f'{len(jobs)} runs classified')


def tables():
    sh = pd.read_csv(out / 'shares.csv.gz')
    fig = pd.read_csv(repo / 'common' / 'fig3a_malawi.csv', comment='#')
    fig = fig[fig.model == 'EMOD'].pivot(index='year', columns='category', values='share')
    rows = []
    for (v, y), d in sh[sh.year.isin(YEARS)].groupby(['variant', 'year']):
        for c in CATS:
            x = 100 * d[c]
            rows.append(dict(variant=v, year=y, category=c, n_runs=len(d), median=x.median(), lo=x.quantile(0.05), hi=x.quantile(0.95),
                             mean=x.mean(), figure=fig.loc[y, c]))
    t = pd.DataFrame(rows)
    t['diff'] = t['mean'] - t.figure
    t.round(2).to_csv(out / 'shares_vs_figure.csv', index=False)

    tot = pd.read_csv(out / 'totals.csv.gz')
    tot['factor'] = CENSUS / tot.pop2018
    # Annualised: sexual transmissions by the share of the year's monthly steps present; report by half-years present
    tot['sexual_x100'] = tot.sexual / SCALE * 12 / tot.steps
    tot['sexual_census'] = tot.sexual * tot.factor * 12 / tot.steps
    tot['ni_15_64_x100'] = tot.ni_15_64 / SCALE * 2 / tot.halves
    tot['ni_15_64_census'] = tot.ni_15_64 * tot.factor * 2 / tot.halves
    cols = ['sexual_x100', 'sexual_census', 'ni_15_64_x100', 'ni_15_64_census', 'factor']
    tt = tot[tot.year.isin(YEARS)].groupby(['variant', 'year'])[cols].median().round(1)
    tt.to_csv(out / 'totals_by_variant.csv')
    return t, tt


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sh = pd.read_csv(out / 'shares.csv.gz')
    fig3 = pd.read_csv(repo / 'common' / 'fig3a_malawi.csv', comment='#')
    fig3 = fig3[fig3.model == 'EMOD'].pivot(index='year', columns='category', values='share')
    colors = dict(undiagnosed='#8fd0e6', diagnosed_never_treated='#f5b800', on_art='#9ccc65', interrupted_art='#c9585c')
    panels = [('figure', 'Published EMOD (Figure 3A, digitised)')] + [(v, f'Our EMOD: {LABELS[v]}') for v in VARIANTS]
    f, axs = plt.subplots(len(panels), 1, figsize=(8, 2.2 * len(panels)), sharex=True)
    for ax, (v, title) in zip(axs, panels):
        if v == 'figure':
            m = fig3.loc[2000:2024, CATS]
        else:
            d = sh[(sh.variant == v) & sh.year.between(2000, 2024)]
            m = 100 * d.groupby('year')[CATS].mean()
        bottom = np.zeros(len(m))
        for c in CATS:
            ax.bar(m.index, m[c], bottom=bottom, color=colors[c], width=0.7, label=c.replace('_', ' '))
            bottom += m[c].to_numpy()
        ax.set_title(title, fontsize=10, loc='left')
        ax.set_ylim(0, 100)
        ax.set_ylabel('% of new infections')
        ax.spines[['top', 'right']].set_visible(False)
    axs[0].legend(ncol=4, fontsize=8, loc='lower left', bbox_to_anchor=(0, 1.15), frameon=False)
    f.tight_layout()
    (here / 'figures').mkdir(exist_ok=True)
    f.savefig(here / 'figures' / 'shares_vs_figure.png', dpi=130)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=40)
    a = ap.parse_args()
    if a.plot:
        t, tt = tables()
        pd.set_option('display.width', 200)
        print(t[t.year.isin([2020, 2024])].round(1).to_string(index=False))
        print(tt)
        plot()
    else:
        reduce(a.cores)
