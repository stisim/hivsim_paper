"""
Analysis 25: run the 129 distinct parameter vectors among Daniel's 250 resampled sets, with both PR #4 fixes, 2 seeds each.
Exactly analysis 18's "both fixes" runs (run_emod.py --transmission --inputs <both-fixes InputFiles>, to 2025), but seed-major
job order so every set has one seed before any has two. Writes raw_results/emod/a25_both/setNNN_sK and outputs/set_map.csv.

    python launch.py INPUTS_DIR [--cores 140]
"""
import sys
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo / 'models' / 'emod' / 'reference'))
import run_emod as r # noqa: E402

META = ['iteration_number', 'likelihood', 'parameterization_id', 'run_number', 'sim_id']
SEEDS = [1, 2]
OUT = 'a25_both'


def set_map():
    """Each of the 250 sets (Run_Number order) with the first parameterization_id that has the same parameter vector"""
    s = pd.read_csv(r.orig / 'resampled_parameter_sets_20230424.csv')
    p = [c for c in s.columns if c not in META]
    key = s[p].round(12).astype(str).agg('|'.join, axis=1)
    s['rep_pid'] = s.groupby(key).parameterization_id.transform('min')
    return s


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs')
    ap.add_argument('--cores', type=int, default=140)
    a = ap.parse_args()
    s = set_map()
    (here / 'outputs').mkdir(exist_ok=True)
    s[META + ['rep_pid']].to_csv(here / 'outputs' / 'set_map.csv', index=False)
    reps = s[s.parameterization_id == s.rep_pid]
    pmap = r.load_param_map()
    outdir = r.outroot / OUT
    jobs = [(int(row.parameterization_id), r.map_sample(row.drop(META + ['rep_pid']).to_dict(), pmap), row.run_number, 2025, outdir,
             False, False, False, [], True, seed, a.inputs) for seed in SEEDS for _, row in reps.iterrows()]
    print(f'{len(reps)} distinct sets x {len(SEEDS)} seeds = {len(jobs)} runs', flush=True)
    with ProcessPoolExecutor(min(a.cores, len(jobs))) as ex:
        for pid, rc, dt in ex.map(r.run_one, jobs):
            print(f'set {pid}: {"ok" if rc == 0 else f"FAILED rc={rc}"} ({dt / 60:.1f} min)', flush=True)
