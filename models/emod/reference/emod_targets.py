"""
Measure EMOD runs against common/targets.csv and draw the shared target figures (common/targetfit.py).

    python emod_targets.py --out sub10   # raw_results/emod/sub10 -> analyses/01_reference_fits/{outputs,figures}/emod_*

Stocks are read from the ReportHIVByAgeAndGender row nearest the target year (rows are every half year); flows are
summed over the calendar year. Female sex workers are women with Risk=HIGH (set by the commercial-sex campaign
events); suppressed means ARTStatus=OnEffectiveART.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
root = here.parents[2]
sys.path.insert(0, str(root / 'common'))
import targetfit  # noqa: E402

SCALE = 100  # 1% population
PSETS = here / 'source' / 'resampled_parameter_sets_20230424.csv'


def load_report(f):
    d = pd.read_csv(f, skipinitialspace=True)
    d.columns = d.columns.str.strip()
    d = d.rename(columns={'IP_Key:Risk': 'Risk', 'IP_Key:ARTStatus': 'ARTStatus'})
    return d.assign(sex=d.Gender.map({0: 'm', 1: 'f'}))


def measure(d, targets):
    years = np.unique(d.Year)
    out = {}
    for tid, r in targets.iterrows():
        m = (d.Age >= r.age_lo) & (d.Age < r.age_hi)
        if r.group == 'fsw':
            m &= d.Risk == 'HIGH'
        elif r.group != 'all':
            raise ValueError(r.group)
        if r.sex in ('f', 'm'):
            m &= d.sex == r.sex
        if r.target in ('new_infections', 'aids_deaths'):
            col = 'Newly Infected' if r.target == 'new_infections' else 'Died_from_HIV'
            out[tid] = d.loc[m & (np.floor(d.Year) == np.floor(r.year)), col].sum() * SCALE
            continue
        s = d[m & (d.Year == years[np.abs(years - r.year).argmin()])]
        inf, pos = s.Infected.sum(), s[s.HasHIV == 1]
        supp = pos.loc[pos.ARTStatus == 'OnEffectiveART', 'On_ART'].sum()
        if r.target == 'incidence':  # new infections in the year to the survey, per susceptible person at the survey
            new = d.loc[m & (d.Year > r.year - 1) & (d.Year <= r.year), 'Newly Infected'].sum()
            out[tid] = new / (s.Population.sum() - inf)
        else:
            out[tid] = {
                'prevalence': lambda: inf / s.Population.sum(),
                'plhiv': lambda: inf * SCALE,
                'on_art': lambda: s.On_ART.sum() * SCALE,
                'population': lambda: s.Population.sum() * SCALE,
                'art_coverage': lambda: s.On_ART.sum() / inf,
                'diagnosed': lambda: pos.Diagnosed.sum() / inf,
                'on_ART_among_diagnosed': lambda: pos.On_ART.sum() / pos.Diagnosed.sum(),
                'suppressed_among_on_ART': lambda: supp / s.On_ART.sum(),
                'vls_among_plhiv': lambda: supp / inf,
            }[r.target]()
    return {int(k): float(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='sub10', help='folder under raw_results/emod/')
    a = ap.parse_args()
    targets = targetfit.load_targets()
    psets = pd.read_csv(PSETS).set_index('parameterization_id')
    expdir = root / 'analyses' / '01_reference_fits'
    (expdir / 'outputs').mkdir(parents=True, exist_ok=True)
    rows = []
    for i, f in enumerate(sorted((root / 'raw_results' / 'emod' / a.out).glob('set*/output/ReportHIVByAgeAndGender.csv'))):
        pid = int(f.parts[-3][3:])
        rows.append(dict(i=i, kind='sample', pars=dict(set=pid, likelihood=float(psets.loc[pid, 'likelihood'])),
                         values=measure(load_report(f), targets)))
    (expdir / 'outputs' / 'emod_results.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in rows))
    targetfit.plot(expdir, 'emod')


if __name__ == '__main__':
    main()
