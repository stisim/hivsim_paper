""" Age structure of the derived populations against Optima's (shares of each sex, in %): python age_structure.py > outputs/age_structure.txt """
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.argv=['x']
import run as R, numpy as np, pandas as pd, poplabel as pl
ode, ag, n = R.optima_sizes()
for m in ('emod','hivsim'):
    d = pd.read_csv(R.OUT/f'sizes_{m}.csv', index_col='year')
    for y in (1990, 2000, 2010, 2020):
        o, x = ode.loc[y], d.loc[y]
        rows=[]
        for sex in 'FM':
            ps=[p for p in pl.POPS if p[0]==sex]
            so, sx = 100*o[ps]/o[ps].sum(), 100*x[ps]/x[ps].sum()
            gen=[p for p in ps if p[1].isdigit()]
            young=[p for p in gen if int(p[1:].split('-')[0].rstrip('+'))<15]; mid=[p for p in gen if 15<=int(p[1:].split('-')[0].rstrip('+'))<50]; old=[p for p in gen if int(p[1:].split('-')[0].rstrip('+'))>=50]
            rows.append(f'{sex}: 0-14 {so[young].sum():.1f}/{sx[young].sum():.1f}, 15-49 general {so[mid].sum():.1f}/{sx[mid].sum():.1f}, 50+ {so[old].sum():.1f}/{sx[old].sum():.1f}, 75+ {so[[p for p in old if int(p[1:].split("-")[0].rstrip("+"))>=75]].sum():.1f}/{sx[[p for p in old if int(p[1:].split("-")[0].rstrip("+"))>=75]].sum():.1f}, max |band diff| {np.abs(so[gen]-sx[gen]).max():.2f} ({np.abs(so[gen]-sx[gen]).idxmax()})')
        print(m, y, 'total', f'{x.sum()/o.sum():.2f}x', ' | '.join(rows))
