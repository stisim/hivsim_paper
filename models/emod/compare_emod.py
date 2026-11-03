"""
Compare rerun EMOD outputs with Daniel's calibration targets (original/calib_targets_data).

    python compare_emod.py --out test5
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

here = Path(__file__).resolve().parent
tdir = here / 'original' / 'calib_targets_data'
outroot = here.parents[1] / 'raw_results' / 'emod'
SCALE = 100  # 1% population


def load(outdir):
    dfs = []
    for f in sorted(outdir.glob('set*/output/ReportHIVByAgeAndGender.csv')):
        d = pd.read_csv(f, skipinitialspace=True)
        d.columns = d.columns.str.strip()
        d['set'] = int(f.parts[-3][3:])
        dfs.append(d)
    return pd.concat(dfs)


def national(d):
    """Yearly national series per set: one row per (set, Year)."""
    d = d.assign(Sex=d.Gender.map({0: 'M', 1: 'F'}), adult=(d.Age >= 15) & (d.Age < 50), a15=d.Age >= 15)
    out = {}
    for sex in 'MF':
        s = d[d.Sex == sex]
        g = s[s.adult].groupby(['set', 'Year'])
        out[f'prev_{sex}'] = g.Infected.sum() / g.Population.sum()
        out[f'artcov_{sex}'] = g.On_ART.sum() / g.Infected.sum()
        out[f'onart15_{sex}'] = s[s.a15].groupby(['set', 'Year']).On_ART.sum() * SCALE
    g = d.groupby(['set', 'Year'])
    out['pop'] = g.Population.sum() * SCALE
    out['plhiv'] = g.Infected.sum() * SCALE
    out = pd.DataFrame(out).reset_index()
    # Flows are reported per half year: sum both halves into the calendar year, shown on the mid-year row
    fy = d.assign(cy=d.Year.astype(int)).groupby(['set', 'cy'])[['Newly Infected', 'Died_from_HIV']].sum() * SCALE
    fy = fy.rename(columns={'Newly Infected': 'new_inf', 'Died_from_HIV': 'hiv_deaths'}).reset_index()
    out = out.assign(cy=out.Year.astype(int)).merge(fy, on=['set', 'cy'], how='left')
    out.loc[out.Year % 1 > 0.5, ['new_inf', 'hiv_deaths']] = None
    return out.drop(columns='cy')


def targets():
    prev = pd.read_csv(tdir / 'MalawiNationalPrevalence.csv')
    prev = prev[(prev.Province.str.strip() == 'All') & (prev.AgeBin == '[15:50)')]
    onart = pd.read_csv(tdir / 'Malawi_OnART.csv')
    cov = pd.read_csv(tdir / 'ART_Coverage_PHIA.csv')
    pop = pd.read_csv(tdir / 'Population.csv')
    pop = pop[pop.Province.str.strip() == 'All'].groupby('Year').Population.sum()
    return prev, onart, cov, pop


def plot(nat, outdir):
    prev, onart, cov, pop = targets()
    fig, axs = plt.subplots(3, 3, figsize=(14, 11))
    panels = [('prev_F', 'Prevalence 15-49, women'), ('prev_M', 'Prevalence 15-49, men'), ('pop', 'Population (all ages)'),
              ('artcov_F', 'ART coverage 15-49, women'), ('artcov_M', 'ART coverage 15-49, men'), ('plhiv', 'PLHIV'),
              ('onart15_F', 'On ART 15+, women'), ('onart15_M', 'On ART 15+, men'), ('new_inf', 'New infections')]
    for ax, (col, title) in zip(axs.flat, panels):
        for s, d in nat.groupby('set'):
            d = d.dropna(subset=[col])
            ax.plot(d.Year, d[col], color='C0', alpha=0.4, lw=0.8)
        sex = {'F': 'Female', 'M': 'Male'}.get(col[-1])
        if col.startswith('prev'):
            t = prev[prev.Gender == sex]
            ax.errorbar(t.Year, t.NationalPrevalence, yerr=[t.NationalPrevalence - t.lower, t.upper - t.NationalPrevalence], fmt='ko', ms=4)
        elif col.startswith('artcov'):
            t = cov[(cov.gender == sex) & (cov.stratification_category.str.contains('15-49', na=False) | (cov.AgeBin == '[15:50)'))]
            ax.errorbar([2015.5] * len(t), t.ARTCoverage / 100, yerr=[(t.ARTCoverage - t.confidence_lower) / 100, (t.confidence_upper - t.ARTCoverage) / 100], fmt='ko', ms=4)
        elif col.startswith('onart'):
            t = onart[onart.Gender == sex]
            ax.errorbar(t.Year, t.On_ART, yerr=t.two_sigma, fmt='ko', ms=3, label='target (15-64)')
        elif col == 'pop':
            ax.plot(pop.index, pop.values, 'ko', ms=4, label='census, ages 0-69 only')
        ax.set_title(title)
        ax.set_xlim(2000, 2025)
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(outdir / 'compare_targets.png', dpi=110)
    print('wrote', outdir / 'compare_targets.png')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='baseline')
    a = ap.parse_args()
    outdir = outroot / a.out
    nat = national(load(outdir))
    nat.to_csv(outdir / 'national.csv', index=False)
    plot(nat, outdir)
    yrs = nat[nat.Year.round(2).isin([2004.42, 2010.42, 2015.42, 2020.42])]  # report rows nearest mid-year
    print(yrs.groupby('Year')[['prev_F', 'prev_M', 'artcov_F', 'artcov_M', 'plhiv', 'new_inf', 'hiv_deaths', 'pop']].median().round(3).to_string())


if __name__ == '__main__':
    main()
