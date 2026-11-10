"""
EMOD′ with Daniel's two fixes (models/emod/v4/components_fixed.py FIXES; condensed cascade, as analysis 17), sets 0/25/50 x
seeds 1..10, 1960-2025, with common/outputs.py: yearly outputs (agents x 100, and prevalence 15-49) into data/emodp_yearly.csv.
    python run_emodp.py
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[2]
sys.path[:0] = [str(repo / 'models' / 'emod' / 'v4'), str(repo)]
PIDS, SEEDS = (0, 25, 50), range(1, 11)


def one(args):
    pid, seed = args
    import harness
    from components_fixed import FIXES
    from cascade_host import EmodHost
    from common.outputs import Outputs, yearly
    sim = harness.make_sim(pid, seed=seed, stop=2025, cascade=True, condensed=True, record=False, fixes=FIXES, analyzers=[Outputs(EmodHost, scale=100)]).run()
    out = sim.analyzers.outputs
    y = yearly(out)
    for c in ['new_infections', 'hiv_deaths', 'pop', 'plhiv', 'diagnosed', 'on_art', 'suppressed']:
        y[c] = y[c]*100
    st = out.strata[(out.strata.age >= 15) & (out.strata.age < 50)].assign(yr=lambda d: np.floor(d.year + 1e-9).astype(int))
    last = st[st.year == st.groupby('yr').year.transform('max')].groupby('yr')[['infected', 'pop']].sum()
    y['prev_15_49'] = y.year.map(last.infected/last['pop'])
    return y.assign(pid=pid, seed=seed)


if __name__ == '__main__':
    with Pool(30) as p:
        res = p.map(one, [(p_, s) for p_ in PIDS for s in SEEDS])
    pd.concat(res).to_csv(here / 'data' / 'emodp_yearly.csv', index=False, float_format='%.6g')
    print('done')
