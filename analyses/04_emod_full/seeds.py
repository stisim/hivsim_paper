"""
How many people EMOD seeds in each node, for parameter set 0, over Run_Numbers 1-20 (short runs to 1982, with the
campaign reduced to what precedes the seeding), against EMOD′'s seeding over 50 seeds. Output: outputs/seeds.csv.

    python seeds.py --emod   # run EMOD (raw_results/emod/seeds_set000/runNN), about 2 minutes each
    python seeds.py          # count seeds in both
"""
import sys
import json
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo / 'models' / 'emod' / 'reference'))
sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
PID, RUNS, SEEDS = 0, range(1, 21), range(1, 51)
OUT = repo / 'raw_results' / 'emod' / 'seeds_set000'


def run_emod(run):
    import os, subprocess
    import run_emod as r
    sets = pd.read_csv(r.orig / 'resampled_parameter_sets_20230424.csv')
    row = sets[sets.parameterization_id == PID].iloc[0]
    meta = ['iteration_number', 'likelihood', 'parameterization_id', 'run_number', 'sim_id']
    files = r.build_inputs(r.map_sample(row.drop(meta).to_dict(), r.load_param_map()), run, years=1982, nat_hist=True)
    files['config.json']['parameters'].update(Report_HIV_ByAgeAndGender_Start_Year=1980, Report_Transmission=1)
    d = OUT / f'run{run:02d}'
    d.mkdir(parents=True, exist_ok=True)
    for k, v in files.items():
        json.dump(v, open(d / k, 'w'), indent=1)
    env = dict(os.environ, LD_LIBRARY_PATH=f"{r.bindir / 'lib'}:{Path.home() / 'miniconda/envs/emod/lib'}", UCX_TLS='self,sm,tcp', UCX_NET_DEVICES='lo')
    with open(d / 'stdout.txt', 'w') as log:
        return subprocess.run([str(r.bindir / 'Eradication'), '-C', 'config.json', '-I', '.', '-O', 'output'], cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT).returncode


def emod_seeds(run):
    """ Infected people in the first report after the seeding, less those infected through relationships before it """
    out = OUT / f'run{run:02d}' / 'output'
    d = pd.read_csv(out / 'ReportHIVByAgeAndGender.csv', skipinitialspace=True)
    d.columns = d.columns.str.strip()
    t = pd.read_csv(out / 'TransmissionReport.csv', usecols=['YEAR', 'NODE_ID'])
    first = d[d.Infected > 0].Year.min()
    inf = d[d.Year == first].groupby('NodeId').Infected.sum()
    tx = t[t.YEAR <= first].groupby('NODE_ID').size()
    return {n: int(inf.get(n, 0) - tx.get(n, 0)) for n in (1, 2, 3)}


def v4_seeds(seed):
    import harness
    sim = harness.make_sim(PID, seed=seed, stop=1981.5, hiv=True)
    sim.run()
    inf = sim.diseases.hiv.infections
    s = inf[inf.type == -1].groupby('node').size()
    return {n: int(s.get(n - 1, 0)) for n in (1, 2, 3)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--emod', action='store_true')
    a = ap.parse_args()
    if a.emod:
        with ProcessPoolExecutor(10) as ex:
            print(list(ex.map(run_emod, RUNS)))
    else:
        with ProcessPoolExecutor(20) as ex:
            v = list(ex.map(v4_seeds, SEEDS))
        e = [emod_seeds(r) for r in RUNS]
        df = pd.concat([pd.DataFrame(e).assign(model='emod', run=list(RUNS)), pd.DataFrame(v).assign(model='v4', run=list(SEEDS))])
        (here / 'outputs').mkdir(parents=True, exist_ok=True)
        df.to_csv(here / 'outputs' / 'seeds.csv', index=False)
        print(df.groupby('model')[[1, 2, 3]].agg(['mean', 'std']).round(2))
        print('share with no seed in node 1:', df.groupby('model')[1].apply(lambda x: (x == 0).mean()).round(2).to_dict())
