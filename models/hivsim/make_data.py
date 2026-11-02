"""
Make the Malawi calibration data and n_art.csv from the UNAIDS estimates in raw_data/ (not committed)

    python make_data.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

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


if __name__ == '__main__':
    df = pd.DataFrame({col: read_unaids(prefix) for prefix,col in series.items()})
    df['hiv.p_on_art'] /= 100
    df.index.name = 'year'
    df.to_csv(datadir/'malawi_hiv_calib.csv')

    n_art = df['hiv.on_art'].dropna().rename('n_art')
    n_art = pd.concat([pd.Series({2003: 0.0}), n_art]).rename('n_art') # ART began in Malawi in 2004
    n_art.index.name = 'year'
    n_art.to_csv(datadir/'n_art.csv')
    print(df.to_string())
