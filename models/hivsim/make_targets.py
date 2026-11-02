"""
Make results/targets.csv, the P3 calibration targets for HIVsim Malawi (see docs/targets.md), from data/sources/

    python make_targets.py

One row per target value: what it is (`target`), when (`year`, the survey's fieldwork midpoint), who (`sex`,
`age_lo`, `age_hi`, `group`), the value with its uncertainty (`lower`, `upper` as a 95% interval, and `n`, the
effective sample size for a survey proportion, where known), where it is from (`source`), and its `role`:
primary or secondary target, or a check that isn't fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

here = Path(__file__).parent
srcdir = here / 'data' / 'sources'
outfile = here.parents[1] / 'results' / 'targets.csv'

# Year labels in EMOD's files -> (survey, fieldwork midpoint)
emod_surveys = {2004: ('DHS 2004', 2004.9), 2010: ('DHS 2010', 2010.7), 2015: ('MPHIA 2015-16', 2016.2), 2015.5: ('DHS 2015-16', 2015.9)}
mphia2 = ('MPHIA 2020-21', 2020.5) # Fieldwork Jan-Mar 2020 and Mar-Apr 2021

columns = ['target', 'year', 'sex', 'age_lo', 'age_hi', 'group', 'value', 'lower', 'upper', 'n', 'source', 'role']


def agebin(s):
    """ '[15:20)' -> (15, 20) """
    lo, hi = s.strip('[)').split(':')
    return int(lo), int(hi)


def jeffreys(p, n):
    """ 95% Jeffreys interval for a proportion p from n people """
    x = p * n
    return stats.beta.ppf(0.025, x + 0.5, n - x + 0.5), stats.beta.ppf(0.975, x + 0.5, n - x + 0.5)


def emod_rows(filename, target, valcol, lowcol='lower', highcol='upper', role='primary', surveys=None):
    """ Rows from one of EMOD's target files; surveys maps EMOD's year labels to (survey, fieldwork midpoint) """
    df = pd.read_csv(srcdir/filename)
    rows = []
    for _, r in df.iterrows():
        lo, hi = agebin(r.AgeBin)
        survey, year = (surveys or {}).get(r.Year, (None, r.Year))
        rows.append(dict(target=target, year=year, sex=r.Gender[0].lower(), age_lo=lo, age_hi=hi, group='all', value=r[valcol],
                         lower=r[lowcol], upper=r[highcol], n=r.get('effective_count', np.nan), source=survey, role=role))
    return pd.DataFrame(rows)


def prevalence():
    """ DHS 2004, 2010, 2015-16 and MPHIA 2015-16 (EMOD's targets), and MPHIA 2020-21, by sex and 5-year age group """
    emod = emod_rows('emod_MalawiNationalPrevalence.csv', 'prevalence', 'NationalPrevalence', surveys=emod_surveys)
    emod = emod[emod.age_hi - emod.age_lo == 5] # The 5-year groups, not the 15-49 totals
    m = pd.read_csv(srcdir/'mphia_2020_21.csv').query('table == "6.3 prevalence"')
    lower, upper = jeffreys(m.value.values, m.n.values) # From the number tested: no design effect, so narrower than the survey's own intervals
    mphia = pd.DataFrame(dict(target='prevalence', year=mphia2[1], sex=m.sex.values, age_lo=m.age_lo.values, age_hi=m.age_hi.values, group='all',
                              value=m.value.values, lower=lower, upper=upper, n=m.n.values, source=mphia2[0], role='primary'))
    return pd.concat([emod, mphia])


def incidence():
    """ MPHIA 2015-16 and 2020-21, ages 15-49, by sex """
    df = emod_rows('emod_Malawi_Incidence_new.csv', 'incidence', 'Incidence', 'lb', 'ub', role='secondary')
    df = df[(df.age_lo == 15) & (df.age_hi == 50)].copy()
    df['year'] = df.year.map({2015.5: 2016.2, 2020: mphia2[1]}) # EMOD labels MPHIA 2015-16 incidence 2015.5
    df['source'] = df.year.map({2016.2: 'MPHIA 2015-16', mphia2[1]: mphia2[0]})
    df['n'] = np.nan # EMOD's effective counts for incidence aren't sample sizes
    return df


def unaids():
    """ UNAIDS estimates (2025 round), all ages, 1990-2021 """
    df = pd.read_csv(here/'data'/'malawi_hiv_calib.csv', index_col='year').loc[1990:2021]
    names = {'hiv.infected': 'plhiv', 'hiv.new_infections': 'new_infections', 'hiv.new_deaths': 'aids_deaths'}
    rows = [dict(target=name, year=year + 0.5, sex='all', age_lo=0, age_hi=100, group='all', value=v, source='UNAIDS 2025', role='secondary')
            for col, name in names.items() for year, v in df[col].items()] # Annual flows and mid-year stocks
    return pd.DataFrame(rows)


def fsw():
    """ Prevalence among female sex workers (Optima's data) """
    df = pd.read_csv(srcdir/'optima_fsw_prevalence.csv')
    return pd.DataFrame(dict(target='prevalence', year=df.year + 0.5, sex='f', age_lo=15, age_hi=50, group='fsw', value=df.value,
                             lower=df.lower, upper=df.upper, source='FSW surveys (Optima data)', role='secondary'))


def checks():
    """ Not fitted: numbers on ART by sex (EMOD's targets; ART is an input as a total), the MPHIA cascade, and population """
    art = emod_rows('emod_Malawi_OnART.csv', 'on_art', 'On_ART', role='check')
    art['source'] = 'Programme data via Spectrum (EMOD)'
    m = pd.read_csv(srcdir/'mphia_2020_21.csv')
    m = m[m.table.str.startswith('9.1.B')]
    cascade = pd.DataFrame(dict(target=m.table.str[6:].str.replace(' ', '_'), year=mphia2[1], sex=m.sex, age_lo=m.age_lo, age_hi=m.age_hi, group='all',
                                value=m.value, n=m.n, source=mphia2[0], role='check'))
    pop = pd.read_csv(here/'data'/'malawi_population.csv', index_col=0)
    pop = pop[pop.Time.isin([1990, 2000, 2010, 2020])].copy()
    pop['age_lo'] = (pop.AgeStart // 5 * 5).clip(upper=65)
    pop = pop.groupby(['Time', 'Sex', 'age_lo'], as_index=False).Value.sum()
    pop = pd.DataFrame(dict(target='population', year=pop.Time + 0.5, sex=pop.Sex.str[0].str.lower(), age_lo=pop.age_lo,
                            age_hi=np.where(pop.age_lo == 65, 100, pop.age_lo + 5), group='all', value=pop.Value, source='WPP', role='check'))
    return pd.concat([art, cascade, pop])


if __name__ == '__main__':
    targets = pd.concat([prevalence(), incidence(), unaids(), fsw(), checks()], ignore_index=True)[columns]
    outfile.parent.mkdir(exist_ok=True)
    targets.to_csv(outfile, index=False)
    print(targets.groupby(['role', 'target', 'source'], sort=False).size().to_string())
