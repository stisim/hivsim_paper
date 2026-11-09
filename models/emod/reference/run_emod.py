"""
Rerun Daniel Citron's EMOD-HIV Malawi model (source/) for his 250 resampled
parameter sets, using the emod-hiv 2.37.0 binary instead of his cluster build.

Replaces dtk-tools: each parameter set is mapped to model inputs exactly as in
source/optim_script.py (map_sample_to_model_input), written into the
templates at the '__KP_' markers, and EMOD is run directly.

Setup (once): see README.md in this folder.

    python run_emod.py --sets 0 25 50 --out test   # a few sets
    python run_emod.py                              # all 250 sets
    python run_emod.py --sets 0 25 --network --out net   # with relationship reports
    python run_emod.py --sets 0 100 200 --natural-history --out nh3   # only what EMOD′ has so far
    python run_emod.py --sets 0 25 --seeds 1 2 --transmission --out tr   # transmissions by the source's cascade stage
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
orig = here / 'source'
repo = here.parents[2]
bindir = repo / 'raw_results' / 'emod_bin'
outroot = repo / 'raw_results' / 'emod'

CAMPAIGN = 'campaign_Malawi_20251215.json'  # the only row of source/scenarios.csv
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
    """Port of map_sample_to_model_input in source/optim_script.py."""
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


def network_reports(cp):
    """Relationship outputs for comparing the sexual network (output only)."""
    cp.update(Report_HIV_ByAgeAndGender_Start_Year=1960, Report_HIV_ByAgeAndGender_Add_Relationships=1,
              Report_HIV_ByAgeAndGender_Add_Concordant_Relationships=1,
              Report_Relationship_Start=1, Report_Relationship_Start_Individual_Properties=['Risk'],
              Report_Relationship_Start_Include_HIV_Disease_Statistics=0,
              Report_Relationship_Start_Include_Other_Relationship_Statistics=0,
              Report_Relationship_Start_Start_Year=1900, Report_Relationship_Start_End_Year=2200,
              Report_Relationship_Start_Min_Age_Years=0, Report_Relationship_Start_Max_Age_Years=200,
              Report_Relationship_Start_Must_Have_IP_Key_Value='', Report_Relationship_Start_Must_Have_Intervention='',
              Report_Relationship_Start_Node_IDs_Of_Interest=[], Report_Relationship_End=1)


def transmission_reports(cp):
    """Daniel's report settings for transmissions by the source's cascade stage (source/InputFiles/Templates/config_transmission.json),
    from 1960 so each person's diagnosis and ART history is complete, without ReportHIVInfection (a row per infected person
    per step, about 1.5 GB a run; not needed for the categories). Output only"""
    cp.update(Report_Transmission=1, Report_Event_Recorder=1, Report_Event_Recorder_Events=['HIVNewlyDiagnosed', 'StartedART', 'StoppedART'],
              Report_Event_Recorder_Ignore_Events_In_List=0, Report_Event_Recorder_Individual_Properties=['ARTStatus'],
              Report_Event_Recorder_Min_Age_Years=0, Report_Event_Recorder_Max_Age_Years=200, Report_Event_Recorder_Must_Have_IP_Key_Value='',
              Report_Event_Recorder_Must_Have_Intervention='', Report_Event_Recorder_Node_IDs_Of_Interest=[],
              Report_Event_Recorder_PropertyChange_IP_Key_Of_Interest='', Report_Event_Recorder_Start_Year=1960, Report_Event_Recorder_End_Year=2100,
              Report_HIV_ART=1, Report_HIV_ART_Start_Year=1960, Report_HIV_ART_Stop_Year=2100,
              Report_HIV_Mortality=1, Report_HIV_Mortality_Start_Year=2000, Report_HIV_Mortality_Stop_Year=2100)


NATURAL_HISTORY = ['Commercial_', 'ModifyStiCoInfectionStatus', 'TraditionalMaleCircumcision', 'OutbreakIndividual'] # Campaign events EMOD′ has so far


def natural_history(files, cp):
    """Keep only the campaign events for sex work, STI co-infection, traditional circumcision and seeding (no testing, care, ART, VMMC, PrEP or PMTCT), and report infections, transmissions and HIV deaths"""
    camp = files[CAMPAIGN]
    camp['Events'] = [e for e in camp['Events'] if any(k in json.dumps(e) for k in NATURAL_HISTORY) and 'NChooser' not in json.dumps(e)]
    cp.update(Report_HIV_ByAgeAndGender_Start_Year=1960, Report_HIV_ByAgeAndGender_Collect_HIV_Stage_Data=1, Report_HIV_ByAgeAndGender_Collect_HIV_Data=0,
              Report_HIV_ByAgeAndGender_Collect_Circumcision_Data=1, Report_Transmission=1, Report_HIV_Mortality=1)


def build_inputs(table, run_number, years=None, network=False, no_hiv=False, nat_hist=False, report_ip=(), transmission=False, inputs=None):
    """Return {filename: json} for one simulation; inputs: a copy of source/InputFiles to read instead"""
    inputs = Path(inputs) if inputs else orig / 'InputFiles'
    config = json.load(open(inputs / 'Templates' / 'config.json'))
    files = {CAMPAIGN: json.load(open(inputs / 'Templates' / CAMPAIGN))}
    for f in DEMOG:
        files[Path(f).name] = json.load(open(inputs / f))
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
    # Reporting only, to match common/targets.csv: from 1990, 5-year age bins (lower edges) to 65+
    # Risk=HIGH marks current female sex workers; ARTStatus separates suppressive from non-suppressive ART
    cp.update(Report_HIV_ByAgeAndGender_Start_Year=1990, Report_HIV_ByAgeAndGender_Collect_Age_Bins_Data=list(range(0, 70, 5)),
              Report_HIV_ByAgeAndGender_Collect_IP_Data=['Risk', 'ARTStatus'])
    if network:
        network_reports(cp)
    if nat_hist:
        natural_history(files, cp)
    if transmission:
        transmission_reports(cp)
    if report_ip:
        cp['Report_HIV_ByAgeAndGender_Collect_IP_Data'] = cp['Report_HIV_ByAgeAndGender_Collect_IP_Data'] + list(report_ip)
    if no_hiv: # Drop the seeding, the only source of infection, for comparisons of the network without HIV
        camp = files[CAMPAIGN]
        camp['Events'] = [e for e in camp['Events'] if e['Event_Coordinator_Config'].get('Intervention_Config', {}).get('class') != 'OutbreakIndividual']
    if years is not None:
        cp['Simulation_Duration'] = round((years - cp['Base_Year']) * 365)
    files['config.json'] = config
    return files


def run_one(args):
    pid, table, run_number, years, outdir, network, no_hiv, nat_hist, report_ip, transmission, seed, inputs = args
    simdir = outdir / (f'set{pid:03d}' if seed is None else f'set{pid:03d}_s{seed}')
    if seed is not None:
        run_number = (run_number*10 + seed) % 65536 # As EMOD′'s seeds (analyses 05-06), within EMOD's range
    simdir.mkdir(parents=True, exist_ok=True)
    for name, obj in build_inputs(table, run_number, years, network, no_hiv, nat_hist, report_ip, transmission, inputs).items():
        json.dump(obj, open(simdir / name, 'w'), indent=1)
    env = dict(os.environ, LD_LIBRARY_PATH=f"{bindir / 'lib'}:{Path.home() / 'miniconda/envs/emod/lib'}",
               UCX_TLS='self,sm,tcp', UCX_NET_DEVICES='lo')  # Azure NIC breaks MPI init otherwise
    t0 = time.time()
    with open(simdir / 'stdout.txt', 'w') as log:
        rc = subprocess.run([str(bindir / 'Eradication'), '-C', 'config.json', '-I', '.', '-O', 'output'],
                            cwd=simdir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    return f'{pid:3d}' + ('' if seed is None else f' seed {seed}'), rc, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sets', type=int, nargs='*', help='parameterization_ids (default: all)')
    ap.add_argument('--years', type=float, default=2025, help='stop year (the original setup runs to 2100)')
    ap.add_argument('--out', default='baseline', help='folder under raw_results/emod/')
    ap.add_argument('--cores', type=int, default=150)
    ap.add_argument('--network', action='store_true', help='also write relationship reports (RelationshipStart/End, relationships by age and sex from 1960)')
    ap.add_argument('--no-hiv', action='store_true', help='without the epidemic seeding, so nobody is infected')
    ap.add_argument('--report-ip', nargs='*', default=[], help='more individual properties to split the report by, e.g. InterventionStatus')
    ap.add_argument('--natural-history', action='store_true', help='campaign reduced to sex work, STI co-infection, traditional circumcision and seeding; report transmissions and HIV deaths')
    ap.add_argument('--transmission', action='store_true', help="Daniel's transmission reports (TransmissionReport, diagnosis and ART events, ReportHIVART), from 1960")
    ap.add_argument('--seeds', type=int, nargs='*', help='run each set with these seeds (Run_Number = (run_number*10 + seed) mod 65536, folders setNNN_sK); default: its own run_number')
    ap.add_argument('--inputs', help='read model inputs from this copy of source/InputFiles (e.g. with fixes) instead')
    a = ap.parse_args()

    sets = pd.read_csv(orig / 'resampled_parameter_sets_20230424.csv')
    if a.sets:
        sets = sets[sets.parameterization_id.isin(a.sets)]
    meta = ['iteration_number', 'likelihood', 'parameterization_id', 'run_number', 'sim_id']
    pmap = load_param_map()
    outdir = outroot / a.out
    jobs = [(int(r.parameterization_id), map_sample(r.drop(meta).to_dict(), pmap), r.run_number, a.years, outdir, a.network, a.no_hiv, a.natural_history, a.report_ip, a.transmission, seed, a.inputs)
            for _, r in sets.iterrows() for seed in (a.seeds or [None])]
    with ProcessPoolExecutor(min(a.cores, len(jobs))) as ex:
        for pid, rc, dt in ex.map(run_one, jobs):
            print(f'set {pid}: {"ok" if rc == 0 else f"FAILED rc={rc}"} ({dt / 60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()
