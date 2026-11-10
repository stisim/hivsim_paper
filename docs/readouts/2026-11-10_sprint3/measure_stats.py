"""
Basic stats per model for the sprint 3 readout: wall time and peak memory of a single run, alone, 3 times.

    python measure_stats.py run MODEL [REPS [FIRST]]   # REPS runs (numbered from FIRST) one after another, each under /usr/bin/time -v; logs/MODEL_rK.log
    python measure_stats.py one MODEL SEED     # one run in this process (called by `run`); prints a JSON line
    python measure_stats.py loc                # lines of code of each model's own code
    python measure_stats.py table              # logs/runs.csv from logs/
    python measure_stats.py plot               # figures/run_time.png, .svg (after `write`)
    python measure_stats.py write              # model_stats.csv (needs logs/loc.json from `loc`)

MODEL: emod (the binary, Daniel's fixes), emod_prime (burn-in time to 1985 logged within the run), hivsim, hybrid, optima, optima_prime_ode, optima_prime_abm.
Python models run with one thread (NUMBA/OMP/MKL = 1), as the runner does.
"""
import json, os, re, shutil, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LOGS = HERE / 'logs'
ENVS = Path.home() / 'miniconda' / 'envs'
PY = dict(rc3=ENVS / 'hp_rc3' / 'bin' / 'python', optima=ENVS / 'optima' / 'bin' / 'python')
PYTHONPATH = '/tmp/ss-nov09:/tmp/sti-613'
EMOD_RUNS = Path.home() / 'hivsim_paper' / 'raw_results' / 'emod' / 'fig3a_both' # analysis 18's runs with Daniel's two fixes (PR #4)
EMOD_BIN = Path.home() / 'hivsim_paper' / 'raw_results' / 'emod_bin'
STOP = 2025
ENV = {'hp_rc3': 'rc3', 'optima': 'optima'}
MODELS = dict(emod=None, emod_prime='rc3', hivsim='rc3', hybrid='rc3', optima='optima', optima_prime_ode='rc3', optima_prime_abm='rc3')


def one(model, seed):
    """ One run; returns what the model says about itself """
    t0 = time.time()
    info = {}
    if model.startswith('emod_prime'):
        sys.path[:0] = [str(REPO / 'models' / 'emod' / 'v4'), str(REPO)]
        import harness, components
        import starsim as ss
        class Clock(ss.Analyzer): # The wall time when the run reaches the burn-in year (the runner saves there and resumes)
            def step(self):
                if 'burnin_s' not in info and self.sim.t.now('year') >= components.BURNIN:
                    info['burnin_s'] = time.time() - t0
        sim = harness.make_sim(pid=0, seed=seed, stop=STOP, record=False, hiv=True, cascade=True, condensed=True, fixes=('mortality', 'reenrol'), analyzers=[Clock(name='clock')])
        sim.init()
        info |= dict(agents_start=len(sim.people.auids), scale=components.SCALE)
    elif model == 'hivsim':
        sys.path[:0] = [str(REPO / 'models' / 'hivsim'), str(REPO)]
        import swap
        sim = swap.make_sim('own', seed=seed, n_agents=10e3, stop=STOP)
        sim.init()
        info |= dict(agents_start=len(sim.people.auids), total_pop=getattr(sim.pars, 'total_pop', None))
    elif model == 'hybrid': # Analysis 16's HHEE (fixed): HIVsim's D and RN, EMOD′ fixed's V and C; D is HIVsim's, so no burn-in
        sys.path[:0] = [str(REPO)]
        from common import swaps as s
        config = {'from': 'hivsim', 'to': 'emod_fixed', 'slots': ['C', 'V', 'RN', 'D'], 'pairs': {'hivsim': [0], 'emod_fixed': [0]},
                  'seeds': 3, 'stop': STOP, 'models': {'emod_fixed': 'models/emod/v4/components_fixed.py'}}
        os.chdir(REPO)
        sim = s.build(dict(D='hivsim', RN='hivsim', V='emod_fixed', C='emod_fixed'), s.runner.cells(config)[seed - 1], config)
        sim.init()
        info |= dict(agents_start=len(sim.people.auids))
    elif model.startswith('optima_prime'):
        sys.path[:0] = [str(REPO / 'models' / 'optima' / 'v4'), str(REPO)]
        import optima_prime as opp, compare as cmp, agents as ag
        if model.endswith('ode'):
            sim = opp.make_sim(components=cmp.COMPONENTS['E'], end=STOP)
        else:
            sim = ag.make_sim(n_agents=100_000, rand_seed=seed, end=STOP)
        sim.init()
        if model.endswith('abm'):
            info |= dict(agents_start=len(sim.people.auids))
    elif model == 'optima':
        os.chdir(REPO / 'models' / 'optima' / 'reference')
        import optima as op
        P = op.loadproj('Malawi_20220912.prj', verbose=0)
        tl = time.time()
        res = P.runsim(name='20220912', end=STOP, verbose=0)
        info |= dict(load_s=tl - t0, start=float(res.tvec[0]), stop=float(res.tvec[-1]), dt=float(P.settings.dt),
                     pops=len(P.data['pops']['short']), states=len(P.settings.statelabels))
        print('STATS ' + json.dumps(info | dict(wall_s=time.time() - t0)), flush=True)
        return
    ti = time.time()
    info |= dict(init_s=ti - t0, start=float(sim.t.start.years) if hasattr(sim.t.start, 'years') else str(sim.pars.start),
                 stop=str(sim.pars.stop), dt=str(sim.pars.dt))
    sim.run()
    info |= dict(run_s=time.time() - ti, wall_s=time.time() - t0)
    try:
        info['agents_end'] = len(sim.people.auids)
    except Exception:
        pass
    print('STATS ' + json.dumps(info), flush=True)


def run(model, reps=3, first=1):
    LOGS.mkdir(exist_ok=True)
    threads = {k: '1' for k in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS')}
    for r in range(first, first + reps):
        log = LOGS / f'{model}_r{r}.log'
        if model == 'emod': # The binary on analysis 18's set 0 inputs (Daniel's fixes, with its transmission reports), seed r
            simdir = Path('/tmp/hp-ro-stats-emod') / f's{r}'
            shutil.rmtree(simdir, ignore_errors=True)
            shutil.copytree(EMOD_RUNS / f'set000_s{r}', simdir, ignore=shutil.ignore_patterns('output', 'stdout.txt', 'status.txt'))
            env = dict(os.environ, LD_LIBRARY_PATH=f"{EMOD_BIN / 'lib'}:{ENVS / 'emod' / 'lib'}", UCX_TLS='self,sm,tcp', UCX_NET_DEVICES='lo')
            cmd, cwd = ['/usr/bin/time', '-v', str(EMOD_BIN / 'Eradication'), '-C', 'config.json', '-I', '.', '-O', 'output'], simdir
        else:
            env = dict(os.environ, PYTHONPATH=PYTHONPATH, **threads)
            cmd, cwd = ['/usr/bin/time', '-v', str(PY[MODELS[model]]), '-u', __file__, 'one', model, str(r)], HERE
        with open(log, 'w') as f:
            rc = subprocess.run(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
        if model == 'emod': # Keep the log small: drop the step-by-step lines
            lines = log.read_text().splitlines()
            log.write_text('\n'.join([l for l in lines if 'steps complete' not in l and 'Update(): Time:' not in l][-80:]) + '\n')
        print(f'{model} r{r}: rc {rc}', flush=True)


def parse(log):
    s = log.read_text()
    out = {}
    m = re.search(r'Elapsed \(wall clock\) time.*?: (?:(\d+):)?(\d+):([\d.]+)', s)
    if m:
        out['elapsed_s'] = int(m[1] or 0)*3600 + int(m[2])*60 + float(m[3])
    m = re.search(r'Maximum resident set size \(kbytes\): (\d+)', s)
    if m:
        out['maxrss_mb'] = int(m[1])/1024
    m = re.search(r'^STATS (.*)$', s, re.M)
    if m:
        out |= json.loads(m[1])
    return out


def loc():
    """ Non-blank, non-comment lines of Python in each rebuild's own code (not libraries, not tests) """
    def count(paths):
        n = 0
        for p in paths:
            for line in Path(p).read_text().splitlines():
                t = line.strip()
                n += bool(t) and not t.startswith('#')
        return n
    m = REPO / 'models'
    groups = {
        'emod_prime': [p for p in (m / 'emod' / 'v4').glob('*.py') if not p.name.startswith('test_') and p.name not in ('adapters.py',)],
        'emod_prime_adapters': [m / 'emod' / 'v4' / 'adapters.py'],
        'hivsim': [m / 'hivsim' / 'components.py', m / 'hivsim' / 'swap.py', m / 'hivsim' / 'malawi' / 'malawi.py'],
        'optima_prime': [p for p in (m / 'optima' / 'v4').glob('*.py') if p.name not in ('acceptance.py', 'compare.py')],
        'common_cascade': sorted((REPO / 'common' / 'cascade').rglob('*.py')),
        'common_swaps': [p for p in sorted((REPO / 'common' / 'swaps').rglob('*.py')) if 'tests' not in p.parts],
    }
    out = {k: dict(files=len(v), loc=count(v)) for k, v in groups.items()}
    print(json.dumps(out, indent=1))
    (LOGS / 'loc.json').write_text(json.dumps(out, indent=1) + '\n')
    return out


def table():
    import pandas as pd
    rows = []
    for log in sorted(LOGS.glob('*_r*.log')):
        model, rep = log.stem.rsplit('_r', 1)
        if rep == '0': # Warm-up (Numba compiles and caches on a first run); not counted
            continue
        rows.append(dict(model=model, rep=int(rep)) | parse(log))
    df = pd.DataFrame(rows)
    df.to_csv(LOGS / 'runs.csv', index=False)
    med = df.groupby('model')[['elapsed_s', 'maxrss_mb']].agg(['median', 'min', 'max', 'count'])
    print(df.to_string())
    print(med.to_string())
    return df


# What isn't measured here: from the model's inputs and the analyses (see model_stats.md, Notes)
ROWS = [
    dict(model='EMOD (binary)', key='emod', horizon='1960.5–2025', population='36,600 agents at start (1% scale, 100 people per agent), ~195k by 2025; 3 nodes',
         time_step='1 month (30.4 days)', start_year=1960.5, param_sets="250 resampled (Daniel's calibration; 129 distinct)", stochastic='yes; 10 sets × 5 seeds (analyses 09, 18)',
         loc='n/a (C++ binary, emod-hiv 2.37.0); inputs: Daniel\'s JSON files', fidelity='reference (with PR #4 fixes: Fig 3A shares within 5 points of the published panel; 2024 new infections 41.3k vs 30.6k published, analysis 18)',
         outputs='InsetChart, ReportHIVByAgeAndGender (by node, age, sex, Risk, ARTStatus), transmission/ART event reports'),
    dict(model="EMOD′ fixed (v4)", key='emod_prime', horizon='1960.5–2025', population='36,600 agents at start (100 people per agent), ~198k by 2025; 3 nodes',
         time_step='1 month', start_year=1960.5, param_sets="EMOD's 250 (same sets)", stochastic='yes; 10 sets × 3 seeds (analyses 04–06), 3 sets × 10 seeds (17)',
         loc_key='emod_prime', loc_extra='common_cascade', fidelity="bugs copied: target medians within 0.96–1.003 of EMOD (analysis 04); condensed cascade 343/345 target points within 2 SE of the campaign (06); fixed: 2024 new infections ~42k vs 41.3k for EMOD with the same fixes (17, 18)",
         outputs='EmodReport (EMOD\'s ReportHIVByAgeAndGender table), common Outputs (Fig 3A shares, cascade)'),
    dict(model='HIVsim Malawi (STIsim v4)', key='hivsim', horizon='1985–2025', population='~10,000 agents at start (762 people per agent), ~29k by 2025',
         time_step='1 month', start_year=1985, param_sets='1 (provisional: hivsim_malawi exp 03 best, models/hivsim/reference/pars.json)', stochastic='yes; 5 seeds (fit checks), 200 cells (analysis 16)',
         loc_key='hivsim', fidelity='is the original (STIsim); slot-factory rebuild identical, 0 differing arrays (analysis 11); prevalence 15–49 in 2020 8.1% (analysis 16)',
         outputs='STIsim results (HIV, ART, testing, VMMC), common Outputs'),
    dict(model='Optima HIV 2.12.2 (original)', key='optima', horizon='1990–2025', population='compartments: 38 health states × 39 populations = 1,482',
         time_step='0.2 year', start_year=1990, param_sets='1 (parset 20220912)', stochastic='no (deterministic)', loc='n/a (Optima HIV package)', fidelity='reference',
         outputs='results.main: prevalence, PLHIV, incidence, deaths, diagnosed, on ART, suppressed, popsize (by population)'),
    dict(model="Optima′ ODE (v4)", key='optima_prime_ode', horizon='1990–2025', population='compartments: 1,482 (as Optima)',
         time_step='0.02 year', start_year=1990, param_sets="Optima's 1", stochastic='no (deterministic)', loc_key='optima_prime',
         fidelity='within 2% of Optima each year 2000–2020 overall, except new infections 2017 (2.2%) and 2020 (2.4%) (optima_prime_plan.md)',
         outputs='Optima′ results by state and population; group results'),
    dict(model="Optima′ agents (v4)", key='optima_prime_abm', horizon='1990–2025', population='100,000 agents at start (~1% of Malawi in 1990), ~227k by 2025; 39 population groups',
         time_step='0.02 year', start_year=1990, param_sets="Optima's 1", stochastic='yes; 48 seeds at 100k, 16 at 400k (analysis 20)', loc_key='optima_prime',
         fidelity='means 2000–2020 within 0.35% of the ODE at 400k agents, within 0.65% at 100k (analysis 20)',
         outputs='as Optima′ ODE, counted by group'),
    dict(model='Swap hybrid HHEE (HIVsim D, RN; EMOD′ fixed V, C)', key='hybrid', horizon='1985–2025', population='~10,000 agents at start (HIVsim D: 762 people per agent)',
         time_step='1 month', start_year=1985, param_sets='paired: HIVsim 1 × EMOD′ set 0', stochastic='yes; 200 cells (analysis 16)', loc_key='common_swaps', loc_extra='emod_prime_adapters',
         fidelity='endpoints exact (identity, analysis 16); hybrid has no original', outputs='common Outputs (Fig 3A shares, cascade, PLHIV, deaths)'),
]


def write():
    import pandas as pd
    df = table()
    locs = json.loads((LOGS / 'loc.json').read_text())
    rows = []
    for r in ROWS:
        d = df[df.model == r['key']]
        w = d.elapsed_s
        row = dict(model=r['model'], horizon=r['horizon'],
                   wall_s_median=round(w.median(), 1) if len(d) else None, wall_s_range=f'{w.min():.1f}–{w.max():.1f}' if len(d) else None, runs=len(d),
                   burnin_s_median=round(d.burnin_s.median(), 1) if 'burnin_s' in d and d.burnin_s.notna().any() else None,
                   maxrss_mb_median=round(d.maxrss_mb.median()) if len(d) else None)
        row |= {k: v for k, v in r.items() if k not in ('model', 'key', 'loc_key', 'loc_extra', 'horizon')}
        if 'loc_key' in r:
            l = locs[r['loc_key']]
            row['loc'] = f"{l['loc']} ({l['files']} files)"
            if 'loc_extra' in r:
                row['loc'] += f" + {locs[r['loc_extra']]['loc']} ({r['loc_extra']})"
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(HERE / 'model_stats.csv', index=False)
    print(out.to_string())
    return out


COLORS = dict(emod='#2a78d6', emod_prime='#2a78d6', hivsim='#1baf7a', optima='#eb6834', optima_prime_ode='#eb6834', optima_prime_abm='#eb6834', hybrid='#7f8a99')


def plot():
    """ figures/run_time.png/.svg: median wall time of one run alone, log axis, range over the 3 runs as a whisker """
    import pandas as pd, matplotlib.pyplot as plt
    df = pd.read_csv(LOGS / 'runs.csv')
    df = df[df.rep > 0]
    order = [r['key'] for r in ROWS]
    names = {r['key']: f"{r['model']}\n{r['horizon']}" for r in ROWS}
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for i, k in enumerate(order[::-1]):
        w = df[df.model == k].elapsed_s
        if not len(w):
            continue
        ax.barh(i, w.median(), height=0.55, color=COLORS[k], alpha=0.35 if k.endswith(('_prime', '_ode', '_abm')) or k == 'hybrid' else 1)
        ax.plot([w.min(), w.max()], [i, i], color='#333', lw=1)
        ax.text(w.max()*1.12, i, f'{w.median():.0f} s', va='center', fontsize=9, color='#333')
    ax.set_yticks(range(len(order)), [names[k] for k in order[::-1]], fontsize=8.5)
    ax.set_xscale('log')
    ax.set_xlim(1, 2000)
    ax.set_xlabel('Wall time of one run alone (s, log scale; median of 3, line = range)')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e5e5', lw=1)
    ax.set_axisbelow(True)
    ax.set_title('Run time per model (solid: original; pale: v4 rebuild or hybrid)', fontsize=10, loc='left')
    fig.tight_layout()
    (HERE / 'figures').mkdir(exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(HERE / 'figures' / f'run_time.{ext}', dpi=200)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'one':
        one(sys.argv[2], int(sys.argv[3]))
    elif cmd == 'run':
        run(sys.argv[2], *[int(a) for a in sys.argv[3:5]])
    elif cmd == 'loc':
        loc()
    elif cmd == 'table':
        table()
    elif cmd == 'write':
        write()
    elif cmd == 'plot':
        plot()
