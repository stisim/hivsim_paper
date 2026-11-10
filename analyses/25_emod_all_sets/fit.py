"""
Approximate score of an EMOD run against Daniel's calibration targets (the ingest form's Analyzers and Obs-* sheets).
dtk-tools' HIVCalibSite isn't available; see README.md for the approximation.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats

FORM = Path(__file__).resolve().parents[2] / 'models' / 'emod' / 'reference' / 'source' / 'Data' / 'calibration_ingest_form_Malawi_3_node_DTC.xlsm'
WEIGHTS = dict(Prevalence=0.8, ARTCoverage=0.15, Population=0.025, OnART=0.025) # Analyzers sheet (Incidence 0, provincial prevalence 0)
SITE_POP, SITE_YEAR, SITE_AGES = 8_258_624, 2018, (15, 50) # Site sheet: scaling for the 'Scaling' channels


def load_obs():
    wb = openpyxl.load_workbook(FORM, data_only=True)
    rows = []
    for ch in WEIGHTS:
        for r in wb[f'Obs-{ch}'].iter_rows(min_row=8, max_col=7, values_only=True):
            if isinstance(r[0], (int, float)) and r[3] and isinstance(r[4], (int, float)):
                rows.append(dict(channel=ch, year=r[0], province=r[1], gender=r[2], agebin=r[3], value=r[4],
                                 weight=r[5] if isinstance(r[5], (int, float)) else 1.0, extra=r[6]))
    obs = pd.DataFrame(rows)
    obs = obs[obs.province.astype(str).str.strip() == 'All']
    lohi = obs.agebin.str.strip('[)').str.split(':', expand=True).astype(float)
    obs['lo'], obs['hi'] = lohi[0], lohi[1]
    return obs.reset_index(drop=True)


def score(rep, obs):
    """rep: ReportHIVByAgeAndGender (stripped columns). Returns {channel: weighted mean log density, 'total': weighted sum}"""
    rep = rep.assign(yr=np.floor(rep.Year).astype(int), sex=np.where(rep.Gender == 0, 'Male', 'Female'))
    g = rep.groupby(['yr', 'sex', 'Age'])[['Population', 'Infected', 'On_ART']].sum()
    nhalf = rep.groupby('yr').Year.nunique()
    s2018 = g.loc[SITE_YEAR]
    p1549 = s2018[(s2018.index.get_level_values('Age') >= SITE_AGES[0]) & (s2018.index.get_level_values('Age') < SITE_AGES[1])].Population.sum() / nhalf[SITE_YEAR]
    factor = SITE_POP / p1549
    ages = np.array(sorted(rep.Age.unique()))
    out = {}
    for ch, d in obs.groupby('channel'):
        ll, w = [], []
        for _, o in d.iterrows():
            y = int(np.floor(o.year))
            sel = ages[(ages >= o.lo) & (ages < o.hi)]
            if y not in nhalf.index or not len(sel) or sel.max() + 5 > o.hi: # bin not a union of report bins (e.g. 65+)
                continue
            x = g.loc[y].loc[o.gender].loc[sel].sum()
            if ch == 'Prevalence':
                p, eff = x.Infected / x.Population, o.extra
            elif ch == 'ARTCoverage':
                p, eff = x.On_ART / x.Infected, o.extra
            if ch in ('Prevalence', 'ARTCoverage'):
                ll.append(stats.beta.logpdf(p, eff * o.value + 1, eff * (1 - o.value) + 1))
            else:
                v = (x.Population if ch == 'Population' else x.On_ART) / nhalf[y] * factor
                ll.append(stats.norm.logpdf(v, o.value, o.extra / 2))
            w.append(o.weight)
        out[ch] = np.average(ll, weights=w)
    out['total'] = sum(WEIGHTS[c] * out[c] for c in WEIGHTS)
    return out
