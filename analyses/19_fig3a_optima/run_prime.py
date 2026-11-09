"""Optima′ (all components, dt=0.02, to 2025): sexual infections by the source's Optima state, against our Optima run (env hp_rc3)"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'models' / 'optima' / 'v4'))
import optima_prime as opp, compare as cmp  # noqa: E402

COMPONENTS = cmp.COMPONENTS['E']
STATES = ['undx', 'dx', 'care', 'lost', 'usvl', 'svl']

if __name__ == '__main__':
    sim = opp.make_sim(components=COMPONENTS, dt=0.02, end=2025).run()
    by = opp.infections_by_source(sim)
    by = by.T.groupby(by.columns.str.rsplit('_', n=1).str[0]).sum().T[STATES] # Sum over CD4 stages
    years = np.arange(1991, 2025)
    t = by.index.values
    rows = by.iloc[[np.argmin(abs(t - y)) for y in years]].set_index(years) # The step starting at each whole year, Optima's convention
    rows.index.name = 'year'
    rows.to_csv(HERE / 'outputs' / 'prime_by_source.csv')
    tot = cmp.prime_inci(sim, years)
    print('check: sum over sources / all infection flows (incl. MTCT) by year 2000, 2020:', (rows.sum(1)/tot)[[9, 29]].round(4).tolist())
    o = pd.read_csv(HERE / 'outputs' / 'optima_by_source.csv')
    o = o[(o.counted == 'sexual') & (o.annual == 'step') & (o.mapping == 'raw')].pivot(index='year', columns='category', values='number')
    o = pd.DataFrame({s: o[[c for c in o.columns if c.split('_')[0] == s]].sum(1) for s in STATES})
    sel = (years >= 2000) & (years <= 2020)
    d = 100*(rows[sel]/o.loc[years[sel]] - 1)
    print('Optima′ vs Optima, % difference by source state (sexual), 2000-2020: max abs', d.abs().max().round(2).to_dict())
    print('total:', (100*(rows[sel].sum(1)/o.loc[years[sel]].sum(1) - 1)).abs().max().round(2))
    print(d.loc[[2000, 2005, 2010, 2015, 2020]].round(2))
