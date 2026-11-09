"""
Analysis 24: the cascade step (C) between EMOD′-fixed and Optima′, both ways (see README.md). Run with Starsim and
STIsim on PYTHONPATH: PYTHONPATH=/tmp/sti-613:/tmp/ss-nov09 conda run -n hp_rc3 python run.py ...

    python run.py emod [CORES]          # 1. EMOD′-fixed, own cascade and Optima's table, sets 0/25/50 x 16 seeds, through the runner
    python run.py ode                   # 2. Optima′ ODE: own cascade; EMOD′-fixed's table (sets 0, 25, 50); copied-bug table (set 0); set 0 without flows
    python run.py abm VARIANT SEED      # 2./3. Optima′ as agents (100,000) with Outputs(OptimaHost): own, optima (its table), emodfix0/25/50
    python run.py abms CORES            # every agent run below, in a pool (one process per run)
    python run.py tables                # outputs/*.csv and outputs/report.txt
    python run.py emod_pm [CORES]       # Addendum: own, Optima's table, and both with Optima's PMTCT per person (early.py)
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
OUT, LOGS = HERE / 'outputs', HERE / 'logs'
OPT = REPO / 'models' / 'optima' / 'v4'
EMOD = REPO / 'models' / 'emod' / 'v4'
N_AGENTS, END = 100_000, 2025
PIDS, EMOD_SEEDS = [0, 25, 50], 16
ABM = dict(own=64, optima=64, emodfix0=16, emodfix25=16, emodfix50=16) # Variant: seeds 1..n
FIG3A = ['undiagnosed', 'diagnosed_naive', 'on_art', 'interrupted'] # D8; on_art = art_lt6m + art_supp + art_nonsupp
LABELS = dict(undiagnosed='Undiagnosed', diagnosed_naive='Diagnosed, never treated', on_art='On ART', interrupted='Interrupted')
PUBLISHED = {'EMOD': [38.5, 3.7, 40.6, 17.1], 'Optima': [45.2, 28.2, 17.6, 9.0], 'our Optima (analysis 19)': [42.3, 40.7, 11.4, 5.7]}
STATES = ['undx', 'dx', 'care', 'lost', 'usvl', 'svl']


def emod_config():
    """ EMOD′-fixed with its own C; with Optima's (`optima`); and two sensitivity variants (variants/): with failure and
    regain switching EMOD′'s ART arm (`optima_sw`), and with EMOD′'s own VMMC and PrEP flows as well (`optima_swp`) """
    e = dict(D='emod_fixed', RN='emod_fixed', V='emod_fixed')
    return dict(hybrids=[dict(e, C=c) for c in ('emod_fixed', 'optima', 'optima_sw', 'optima_swp')],
                pairs={'emod_fixed': PIDS}, seeds=EMOD_SEEDS, stop=END,
                models={'emod_fixed': 'models/emod/v4/components_fixed.py', 'optima': 'models/optima/v4/components.py',
                        'optima_sw': 'analyses/24_cascade_swap/variants/optima_c_switch.py', 'optima_swp': 'analyses/24_cascade_swap/variants/optima_c_switch_prev.py'})


def emod_pm_config():
    """ The addendum's: EMOD′-fixed own; Optima's table as it is (`optima`); with Optima's PMTCT number per person
    (`optima_pm`); and with arm switching, EMOD′'s VMMC and PrEP and Optima's PMTCT (`optima_swp_pm`) """
    c = emod_config()
    e = dict(D='emod_fixed', RN='emod_fixed', V='emod_fixed')
    c['hybrids'] = [dict(e, C=k) for k in ('emod_fixed', 'optima', 'optima_pm', 'optima_swp_pm')]
    c['models'] = {k: v for k, v in c['models'].items() if k in ('emod_fixed', 'optima')} | {
        'optima_pm': 'analyses/24_cascade_swap/variants/optima_c_pmtct.py', 'optima_swp_pm': 'analyses/24_cascade_swap/variants/optima_c_switch_prev_pmtct.py'}
    return c


def run_emod(cores=96, config=None):
    from common import swaps
    import starsim, stisim
    print(starsim.__file__, stisim.__file__, flush=True)
    rows = swaps.run(config or emod_config(), cores=cores)
    print(pd.DataFrame(rows).groupby('name').error.apply(lambda e: (e.fillna('') != '').sum()), flush=True)


# %% Optima′

def optima_imports():
    for p in (OPT, REPO / 'models' / 'optima' / 'reference'):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    spec = importlib.util.spec_from_file_location('a20', REPO / 'analyses' / '20_optima_agents' / 'run.py')
    a20 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a20)
    a20.YEARS = np.arange(1990, END + 1) # extract()'s years (it reads the module's YEARS)
    return a20


def emod_table(pid, fixes=('mortality', 'reenrol')):
    """ EMOD′'s condensed cascade table for parameter set `pid`, with the fixes (EMOD′-fixed's: components_fixed.FIXES) """
    sys.path.insert(0, str(EMOD))
    import inputs
    spec = importlib.util.spec_from_file_location('emod_extract_cascade', EMOD / 'extract_cascade.py') # Not flat: Optima's folder has a module of the same name
    ex = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    from common.cascade.params import Params
    params, flows = ex.extract(inputs.load(pid, tuple(fixes)))
    assert not fixes or params.query("pathway == 'dropout' and parameter == 'p_restart'").value.astype(float).item() == 1
    return Params(params, flows)


def cascade_for(variant):
    from common.cascade.params import Params
    if variant == 'own':
        return None
    if variant == 'optima':
        return Params(REPO / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv')
    if variant.startswith('emodfix'):
        return emod_table(int(variant.removeprefix('emodfix').removesuffix('_noflows')))
    if variant == 'emodbug0':
        return emod_table(0, fixes=())
    raise ValueError(variant)


def by_state(sim):
    import optima_prime as opp
    by = opp.infections_by_source(sim)
    yr = np.floor(by.index.to_numpy() + 1e-9).astype(int)
    return by.T.groupby(by.columns.str.split('_').str[0]).sum().T.groupby(yr).sum().reindex(columns=STATES, fill_value=0)*sim.pars.dt.years


def run_ode():
    a20 = optima_imports()
    import optima_prime as opp, compare as cmp
    OUT.mkdir(exist_ok=True)
    for v in ['own', 'emodfix0', 'emodfix25', 'emodfix50', 'emodbug0', 'emodfix0_noflows']:
        T = sc.timer()
        kw = dict(cascade=cascade_for(v)) if v != 'own' else {}
        if v.endswith('_noflows'):
            kw['flows'] = False
        sim = opp.make_sim(components=cmp.COMPONENTS['E'], end=END + 1, group_results=a20.group_results(), **kw).run()
        out = a20.extract(sim)
        out.by_state, out.time = by_state(sim), T.toc(output=True)
        sc.save(OUT / f'ode_{v}.obj', out)
        print(f'ode {v}: {out.time:.0f} s', flush=True)


def run_abm(variant, seed, n=N_AGENTS):
    a20 = optima_imports()
    import agents as ag, host as H
    from common.outputs import Outputs, yearly
    T = sc.timer()
    P = cascade_for(variant)
    intvs = [H.OptimaTracker()]
    if P is not None and P.has('target', 'on_art'): # A table of pathways by rate and level: the agent cascade does the moves
        from common.cascade.pathways import TableCascade
        intvs = [TableCascade(P, H.OptimaHost), H.OptimaTracker()]
    sim = ag.make_sim(n_agents=n, rand_seed=seed, end=END + 1, cascade=P, group_results=a20.group_results(), interventions=intvs, analyzers=[Outputs(H.OptimaHost)])
    sim.run()
    out = sim.analyzers.outputs
    res = sc.objdict(yearly=yearly(out), extract=a20.extract(sim), scale=out.scale, n_agents=n, seed=seed, variant=variant, time=T.toc(output=True))
    sc.save(OUT / 'abm' / f'{variant}_{seed}.obj', res)
    print(f'abm {variant} seed {seed}: {res.time:.0f} s', flush=True)


def _abm_job(job):
    import subprocess
    v, s = job
    if (OUT / 'abm' / f'{v}_{s}.obj').exists():
        return v, s, 'stored'
    with open(LOGS / f'abm_{v}_{s}.log', 'w') as f:
        r = subprocess.run([sys.executable, __file__, 'abm', v, str(s)], stdout=f, stderr=subprocess.STDOUT)
    return v, s, r.returncode


def run_abms(cores):
    from multiprocessing import Pool
    (OUT / 'abm').mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    jobs = [(v, s) for s in range(1, max(ABM.values()) + 1) for v, n in ABM.items() if s <= n]
    with Pool(cores) as pool:
        for v, s, rc in pool.imap_unordered(_abm_job, jobs):
            print(v, s, rc, flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'emod':
        run_emod(int(sys.argv[2]) if len(sys.argv) > 2 else 96)
    elif cmd == 'emod_pm':
        run_emod(int(sys.argv[2]) if len(sys.argv) > 2 else 96, emod_pm_config())
    elif cmd == 'ode':
        run_ode()
    elif cmd == 'abm':
        run_abm(sys.argv[2], int(sys.argv[3]))
    elif cmd == 'abms':
        run_abms(int(sys.argv[2]))
    elif cmd == 'tables':
        import tables
        tables.main()
