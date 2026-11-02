"""
Make the Malawi data that is derived rather than downloaded:

- the calibration targets, from the UNAIDS estimates in raw_data/ (not committed);
- n_art.csv, from the UNAIDS number on ART (2010 on) and Optima's before 2010;
- malawi_deaths.csv, the WPP all-cause death rates (malawi_deaths_all_cause.csv) with AIDS deaths removed:
  at ages 15 and over by dedup_deaths, and at 0-14 by the share of child deaths that UNAIDS estimates were AIDS
  (needs malawi_population.csv, the WPP population by sex and single year of age).

    python make_data.py            # all three
    python make_data.py n_art      # without raw_data/
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from stisim.data.data_utils import dedup_deaths, deleted_fraction

rawdir = Path(__file__).parents[2] / 'raw_data'
datadir = Path(__file__).parent / 'data'
country = 'Malawi'

# UNAIDS file prefix -> HIVsim result
series = {
    'People living with HIV_People living with HIV - All ages': 'hiv.infected',
    'New HIV infections_New HIV infections - All ages': 'hiv.new_infections',
    'AIDS-related deaths_AIDS-related deaths - All ages': 'hiv.new_deaths',
    'Treatment cascade_People living with HIV receiving ART (#)': 'hiv.on_art',
    'Treatment cascade_People living with HIV receiving ART (%)': 'hiv.p_on_art',
}

# Number on ART before 2010, when the UNAIDS series starts: Optima's data (numtx in files/Malawi_20220912.prj).
# 1990-2004 is a straight line from 0 to the 2004 value
optima_n_art = {1990: 0, 1991: 769, 1992: 1537, 1993: 2306, 1994: 3075, 1995: 3843, 1996: 4612, 1997: 5380, 1998: 6149, 1999: 6918,
                2000: 7686, 2001: 8455, 2002: 9224, 2003: 9992, 2004: 10761, 2005: 28110, 2006: 57354, 2007: 96710, 2008: 152298, 2009: 198844}

# AIDS-free anchor years for removing AIDS deaths from adult all-cause mortality (see dedup_deaths)
death_anchors = dict(base_year=1985, end_year=2025)
child_aids = 'AIDS-related deaths_AIDS-related deaths - Children (0-14)' # UNAIDS file prefix
child_ages = 15 # Ages below this use the UNAIDS share


def to_number(val):
    """ '110 000 ' -> 110000; censored values ('<100', '>98') -> NaN """
    val = str(val).replace(' ', '')
    return np.nan if val[:1] in '<>' or val in ('', '...', 'nan') else float(val)


def read_unaids(prefix):
    """ The point estimates for one country, by year """
    path, = rawdir.glob(f'{prefix}*.csv')
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()
    row = df.set_index('Country').loc[country]
    row = row[[c for c in row.index if c.isdigit()]] # Drop the bounds
    return pd.Series([to_number(v) for v in row], index=row.index.astype(int))


def make_targets():
    df = pd.DataFrame({col: read_unaids(prefix) for prefix,col in series.items()})
    df['hiv.p_on_art'] /= 100
    df.index.name = 'year'
    df.to_csv(datadir/'malawi_hiv_calib.csv')
    print(df.to_string())


def make_n_art():
    unaids = pd.read_csv(datadir/'malawi_hiv_calib.csv', index_col='year')['hiv.on_art'].dropna()
    n_art = pd.concat([pd.Series(optima_n_art, dtype=float), unaids]).rename('n_art')
    n_art.index.name = 'year'
    n_art.to_csv(datadir/'n_art.csv')


def child_aids_share(all_cause, pop, aids):
    """
    The share of deaths at ages 0-14 that were AIDS, by year: UNAIDS' AIDS deaths over the deaths that the
    WPP rates give in the WPP population. It rises linearly from 0 in the first year of the rates to the first
    UNAIDS year, and holds at the last UNAIDS year's share after it
    """
    kids = all_cause[all_cause.AgeStart < child_ages].merge(pop, on=['Time', 'Sex', 'AgeStart'], suffixes=('', '_pop'))
    deaths = (kids.Value * kids.Value_pop).groupby(kids.Time).sum()
    share = (aids / deaths).dropna()
    if not share.between(0, 1).all():
        raise ValueError(f'Child AIDS share outside 0-1 (check the population units):\n{share}')
    years = np.unique(all_cause.Time)
    if years[0] < share.index[0]:
        share[years[0]] = 0.0
    return share.sort_index().reindex(years).interpolate().ffill()


def make_deaths():
    all_cause = pd.read_csv(datadir/'malawi_deaths_all_cause.csv', index_col=0)
    non_aids = dedup_deaths(all_cause, **death_anchors)
    share = child_aids_share(all_cause, pd.read_csv(datadir/'malawi_population.csv', index_col=0), read_unaids(child_aids))
    kids = all_cause.AgeStart < child_ages
    non_aids.loc[kids, 'Value'] = all_cause.Value[kids] * (1 - all_cause.Time[kids].map(share))
    non_aids.to_csv(datadir/'malawi_deaths.csv')
    print('Share of child deaths that were AIDS:\n', share.loc[::5].round(3).to_string())
    adults = deleted_fraction(all_cause, non_aids).query(f'AgeStart >= {child_ages}')
    peak = adults.loc[adults.groupby('Sex').aids_share.idxmax(), ['Sex', 'Time', 'AgeStart', 'aids_share']]
    print('Largest share of adult deaths removed as AIDS:\n', peak.round(2).to_string(index=False))


if __name__ == '__main__':
    steps = sys.argv[1:] or ['targets', 'n_art', 'deaths']
    for step in steps:
        dict(targets=make_targets, n_art=make_n_art, deaths=make_deaths)[step]()
