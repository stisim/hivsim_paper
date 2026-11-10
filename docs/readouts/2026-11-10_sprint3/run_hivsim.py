"""
HIVsim reference set (models/hivsim/reference/pars.json, own care), seeds 1..N, to 2025, with common/outputs.py:
yearly outputs (incl. Figure 3A source shares) per seed into data/hivsim_yearly.csv, and targets measured into
data/hivsim_targets.csv.   python run_hivsim.py [seeds] [n_agents]
"""
import sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[2]
sys.path[:0] = [str(repo / 'models' / 'hivsim'), str(repo)]
SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 20
N = float(sys.argv[2]) if len(sys.argv) > 2 else 50e3


def one(seed):
    import swap
    from stisim_host import StisimHost
    from common.outputs import Outputs, yearly, measure
    sim = swap.make_sim('own', seed=seed, n_agents=N, stop=2025, analyzers=[Outputs(StisimHost)]).run()
    out = sim.analyzers.outputs
    y = yearly(out)
    for c in ['new_infections', 'hiv_deaths', 'pop', 'plhiv', 'diagnosed', 'on_art', 'suppressed']:
        y[c] = y[c]*out.scale
    st = out.strata[(out.strata.age >= 15) & (out.strata.age < 50)].assign(yr=lambda d: np.floor(d.year + 1e-9).astype(int))
    last = st[st.year == st.groupby('yr').year.transform('max')].groupby('yr')[['infected', 'pop']].sum()
    y['prev_15_49'] = y.year.map(last.infected/last['pop'])
    targets = pd.read_csv(repo / 'common' / 'targets.csv')
    m = measure(out.strata, targets, out.scale)
    return y.assign(seed=seed), pd.DataFrame(dict(tid=list(m), value=list(m.values()))).assign(seed=seed)


if __name__ == '__main__':
    with Pool(SEEDS) as p:
        res = p.map(one, range(1, SEEDS + 1))
    (here / 'data').mkdir(exist_ok=True)
    pd.concat([r[0] for r in res]).to_csv(here / 'data' / 'hivsim_yearly.csv', index=False, float_format='%.6g')
    pd.concat([r[1] for r in res]).to_csv(here / 'data' / 'hivsim_targets.csv', index=False, float_format='%.6g')
    print('done')
