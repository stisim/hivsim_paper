"""
Rerun Daniel Citron's EMOD-HIV Malawi model (original/) for his 250 resampled
parameter sets, using the emod-hiv 2.37.0 binary instead of his cluster build.

Replaces dtk-tools: each parameter set is mapped to model inputs exactly as in
original/optim_script.py (map_sample_to_model_input), written into the
templates at the '__KP_' markers, and EMOD is run directly.

Setup (once): see README.md in this folder.

    python run_emod.py --sets 0 25 50 --out test   # a few sets
    python run_emod.py                              # all 250 sets
"""
import argparse
import copy
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import openpyxl
import pandas as pd

here = Path(__file__).resolve().parent
orig = here / 'original'
repo = here.parents[1]
bindir = repo / 'raw_results' / 'emod_bin'
outroot = repo / 'raw_results' / 'emod'

CAMPAIGN = 'campaign_Malawi_20251215.json'  # the only row of original/scenarios.csv
DEMOG = ['Static/Demographics_Malawi.json', 'Templates/PFA_Overlay_Malawi.json',
         'Templates/Accessibility_and_Risk_IP_Overlay_Malawi.json',
         'Templates/Risk_Assortivity_Overlay_Malawi.json']
PROVINCES = ['Northern', 'Central', 'Southern']  # Site sheet of the ingest form
SCALE = 0.01  # BASE_POPULATION_SCALE_FACTOR in optim_script.py


def load_param_map():
    """(name, map_to) from the 'Model Parameters' sheet, in sheet order."""
    wb = openpyxl.load_workbook(orig / 'Data' / 'calibration_ingest_form_Malawi_3_node_DTC.xlsm', data_only=True)
    rows = list(wb['Model Parameters'].iter_rows(values_only=True))[1:]
    return [(r[0], r[5].strip() if isinstance(r[5], str) and r[5].strip() else None)
            for r in rows if r[0] not in (None, '--select--')]


def map_sample(sample, param_map):
    """Port of map_sample_to_model_input in original/optim_script.py."""
    s = dict(sample)
    t = {}
    t['Base_Infectivity'] = s.pop('Base Infectivity')
    lo, hi = s.pop('PreART Link Min'), s.pop('PreART Link Max')
    lo, hi = min(lo, hi), max(lo, hi)
    t['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Min'] = lo
    t['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Max'] = hi
    young, old = s.pop('Male To Female Young'), s.pop('Male To Female Old')
    t['Male_To_Female_Relative_Infectivity_Multipliers'] = [young, young, old]
    for p in PROVINCES:
        v = s.pop(f'{p}: LOW Risk')
        t[f'Initial_Distribution__KP_Risk_{p}'] = [v, 1 - v, 0]
    v = s.pop('Risk Assortivity')
    t['Weighting_Matrix_RowMale_ColumnFemale__KP_RiskAssortivity'] = [[v, 1 - v, 0], [1 - v, v, v], [0, v, 1 - v]]
    for name, mapto in param_map:
        if mapto and name in s:
            t[mapto] = s.pop(name)  # later rows overwrite earlier ones with the same map_to, as in the original
    assert not s, f'unused parameters: {list(s)}'
    return t


def find_markers(obj, tag, found):
    if isinstance(obj, dict):
        if tag in obj:
            found.append(obj)
        for v in obj.values():
            find_markers(v, tag, found)
    elif isinstance(obj, list):
        for v in obj:
            find_markers(v, tag, found)
    return found


def set_path(obj, path, value):
    *head, last = path
    for k in head:
        obj = obj[k]
    obj[last] = value


def strip_markers(obj):
    if isinstance(obj, dict):
        for k in [k for k in obj if '__KP_' in k]:
            del obj[k]
        for v in obj.values():
            strip_markers(v)
    elif isinstance(obj, list):
        for v in obj:
            strip_markers(v)


def to_v237(obj):
    """Changes needed for the emod-hiv 2.37 binary; each one leaves the model unchanged."""
    if isinstance(obj, dict):
        if 'Days_To_Diagnosis' in obj:  # removed in 2.37; 0 (no delay) everywhere in this campaign
            assert obj.pop('Days_To_Diagnosis') == 0
        for v in obj.values():
            to_v237(v)
    elif isinstance(obj, list):
        for v in obj:
            to_v237(v)


def build_inputs(table, run_number, years=None):
    """Return {filename: json} for one simulation."""
    config = json.load(open(orig / 'InputFiles' / 'Templates' / 'config.json'))
    files = {CAMPAIGN: json.load(open(orig / 'InputFiles' / 'Templates' / CAMPAIGN))}
    for f in DEMOG:
        files[Path(f).name] = json.load(open(orig / 'InputFiles' / f))
    cp = config['parameters']
    for key, value in table.items():
        tag, *path = key.split('.')
        if '__KP_' in tag:
            realkey = tag.split('__KP_')[0]
            hits = [o for f in files.values() for o in find_markers(f, tag, [])]
            assert hits, f'marker not found: {tag}'
            for o in hits:
                set_path(o, [realkey] + path, copy.deepcopy(value))
        else:
            assert not path and key in cp, f'unknown config parameter: {key}'
            cp[key] = value
    for f in files.values():
        strip_markers(f)
    to_v237(files[CAMPAIGN])
    cp.update(Run_Number=int(run_number), Campaign_Filename=CAMPAIGN, x_Base_Population=SCALE,
              Enable_Demographics_Builtin=0, Demographics_Filenames=[Path(f).name for f in DEMOG])
    # Reporting only, to match results/targets.csv: from 1990, 5-year age bins (lower edges) to 65+
    # Risk=HIGH marks current female sex workers; ARTStatus separates suppressive from non-suppressive ART
    cp.update(Report_HIV_ByAgeAndGender_Start_Year=1990, Report_HIV_ByAgeAndGender_Collect_Age_Bins_Data=list(range(0, 70, 5)),
              Report_HIV_ByAgeAndGender_Collect_IP_Data=['Risk', 'ARTStatus'])
    if years is not None:
        cp['Simulation_Duration'] = round((years - cp['Base_Year']) * 365)
    files['config.json'] = config
    return files


def run_one(args):
    pid, table, run_number, years, outdir = args
    simdir = outdir / f'set{pid:03d}'
    simdir.mkdir(parents=True, exist_ok=True)
    for name, obj in build_inputs(table, run_number, years).items():
        json.dump(obj, open(simdir / name, 'w'), indent=1)
    env = dict(os.environ, LD_LIBRARY_PATH=f"{bindir / 'lib'}:{Path.home() / 'miniconda/envs/emod/lib'}",
               UCX_TLS='self,sm,tcp', UCX_NET_DEVICES='lo')  # Azure NIC breaks MPI init otherwise
    t0 = time.time()
    with open(simdir / 'stdout.txt', 'w') as log:
        rc = subprocess.run([str(bindir / 'Eradication'), '-C', 'config.json', '-I', '.', '-O', 'output'],
                            cwd=simdir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    return pid, rc, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sets', type=int, nargs='*', help='parameterization_ids (default: all)')
    ap.add_argument('--years', type=float, default=2025, help='stop year (the original setup runs to 2100)')
    ap.add_argument('--out', default='baseline', help='folder under raw_results/emod/')
    ap.add_argument('--cores', type=int, default=150)
    a = ap.parse_args()

    sets = pd.read_csv(orig / 'resampled_parameter_sets_20230424.csv')
    if a.sets:
        sets = sets[sets.parameterization_id.isin(a.sets)]
    meta = ['iteration_number', 'likelihood', 'parameterization_id', 'run_number', 'sim_id']
    pmap = load_param_map()
    outdir = outroot / a.out
    jobs = [(int(r.parameterization_id), map_sample(r.drop(meta).to_dict(), pmap), r.run_number, a.years, outdir)
            for _, r in sets.iterrows()]
    with ProcessPoolExecutor(min(a.cores, len(jobs))) as ex:
        for pid, rc, dt in ex.map(run_one, jobs):
            print(f'set {pid:3d}: {"ok" if rc == 0 else f"FAILED rc={rc}"} ({dt / 60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()
