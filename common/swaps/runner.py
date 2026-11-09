"""
The swap runner (docs/swap_design.md 4): hybrids built from each model's slot factories, run in cells (parameter-set
pair x seed) in a process pool, stored with a manifest, and skipped when already run. See README.md.
"""
import os
os.environ.setdefault('NUMBA_NUM_THREADS', '1') # One thread per process (before EmodNetwork imports Numba)
import sys
import json
import time
import hashlib
import itertools
import subprocess
import traceback
import importlib.util
from types import SimpleNamespace
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
ORDER = ['D', 'R', 'N', 'RN', 'V', 'T', 'H', 'S', 'C', 'M'] # The canonical order of slots (modules are added in it)
LISTS = ['demographics', 'networks', 'diseases', 'interventions', 'connectors', 'analyzers'] # Sim keywords whose lists are joined
MODELS = dict(emod='models/emod/v4/components.py', hivsim='models/hivsim/components.py') # Default factory files, from the repo
OUT = REPO / 'raw_results' / 'swaps'
_loaded = {}


def load_model(name, path=None):
    """ The slot-factory module of model `name` (its file, from the repo; its folder goes on the path for its own imports) """
    path = REPO / (path or MODELS[name])
    key = (name, str(path))
    if key not in _loaded:
        sys.path.insert(0, str(path.parent))
        spec = importlib.util.spec_from_file_location(f'components_{name}', path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _loaded[key] = mod
    return _loaded[key]


def name(choice):
    """ A hybrid's name (2.4): `D=hivsim,RN=hivsim,V=emod,C=emod`, slots in the canonical order """
    return ','.join(f'{s}={choice[s]}' for s in sorted(choice, key=ORDER.index))


def hybrids(config):
    """ The hybrids of a path: the sequential path (`slots` swapped from `from` to `to` in that order), or with
    `lattice`, every subset; or `hybrids`, a list of choices given outright """
    if 'hybrids' in config:
        return [dict(h) for h in config['hybrids']]
    slots, a, b = config['slots'], config['from'], config['to']
    subsets = [set(c) for k in range(len(slots) + 1) for c in itertools.combinations(slots, k)] if config.get('lattice') else [set(slots[:k]) for k in range(len(slots) + 1)]
    return [{s: b if s in sub else a for s in slots} for sub in subsets]


def cells(config):
    """ The cells: pair i of `pairs` (a parameter-set ID per model) x seeds 1..`seeds`; the cell's seed is 100*i + s """
    pairs = config.get('pairs', {})
    n = len(next(iter(pairs.values()))) if pairs else 1
    assert all(len(v) == n for v in pairs.values()), 'every model needs one parameter set per pair'
    return [dict(cell=f'p{i}_s{s}', pair=i, pids={m: v[i] for m, v in pairs.items()}, seed=100*i + s)
            for i in range(n) for s in range(1, config.get('seeds', 1) + 1)]


def build(choice, cell, config, outputs=True):
    """
    The sim for hybrid `choice` ({slot: model}) in `cell`: each slot's factory, in the canonical order, called as
    `make_<slot>(pars, ctx)` with its model's parameter set for the cell; their sim keywords merged (module lists
    joined, other keywords must agree); `rand_seed` the cell's; `stop` the config's; and with `outputs`, the
    `common.outputs.Outputs` analyzer last, on the `HOST` of the model of V (or H)
    """
    import starsim as ss
    paths = config.get('models', {})
    models = {m: load_model(m, paths.get(m)) for m in sorted(set(choice.values()))} # Sorted: set order varies between processes, and the load order sets sys.path
    pars = {m: models[m].pars(cell['pids'].get(m)) for m in models}
    ctx = SimpleNamespace(choice=choice, models=models, pars=pars, seed=cell['seed'], stop=config.get('stop'), cell=cell)
    kw = dict(rand_seed=cell['seed'], verbose=0)
    if ctx.stop is not None:
        kw['stop'] = ctx.stop
    for slot in sorted(choice, key=ORDER.index):
        m = choice[slot]
        for k, v in getattr(models[m], f'make_{slot}')(pars[m], ctx).items():
            if k in LISTS:
                kw[k] = kw.get(k, []) + list(v if isinstance(v, (list, tuple)) else [v])
            elif k in kw and kw[k] != v:
                raise ValueError(f'slot {slot} ({m}) sets {k}={v!r}, but it is already {kw[k]!r}')
            else:
                kw[k] = v
    host = getattr(models[choice.get('V', choice.get('H'))], 'HOST', None) if outputs else None
    if host is not None:
        from common.outputs import Outputs
        kw['analyzers'] = kw.get('analyzers', []) + [Outputs(host, scale=getattr(models[choice['D']], 'SCALE', None))]
    sim = getattr(models[choice['D']], 'SIM', ss.Sim) if 'D' in choice else ss.Sim # A model can need its own sim class for its D (HIVsim: sti.Sim, which builds the location's demography)
    return sim(**kw)


def differences(a, b):
    """ The results (and analyzers' `df` and `strata` tables) that differ between two run sims: [] if identical """
    ra, rb = a.results.flatten(), b.results.flatten()
    out = sorted(set(ra.keys()) ^ set(rb.keys()))
    out += [k for k in ra.keys() if k in rb and not np.array_equal(np.asarray(ra[k]), np.asarray(rb[k]), equal_nan=True)]
    for k, an in a.analyzers.items():
        bn = b.analyzers.get(k)
        for attr in ('df', 'strata'):
            x, y = getattr(an, attr, None), getattr(bn, attr, None)
            if isinstance(x, pd.DataFrame) and not (isinstance(y, pd.DataFrame) and x.equals(y)):
                out.append(f'{k}.{attr}')
    return out


def git(*args):
    return subprocess.run(['git', '-C', str(REPO), *args], capture_output=True, text=True).stdout.strip()


def code_version():
    """ The code the hybrids run: git tree IDs of models/ and common/, and a hash of their uncommitted changes """
    dirty = git('diff', 'HEAD', '--', 'models', 'common') + ''.join(
        f'\n{f}\n' + (REPO / f).read_text(errors='replace') for f in git('ls-files', '-o', '--exclude-standard', '--', 'models', 'common').splitlines())
    return dict(tree_models=git('rev-parse', 'HEAD:models'), tree_common=git('rev-parse', 'HEAD:common'),
                dirty=hashlib.sha1(dirty.encode()).hexdigest()[:12] if dirty else '')


def pins():
    """ Version and (for a git checkout) commit of Starsim and STIsim """
    out = {}
    for lib in ('starsim', 'stisim'):
        try:
            mod = __import__(lib)
        except Exception: # Not installed, or not importable with this Starsim
            out[lib] = None
            continue
        d = Path(mod.__file__).parent
        commit = subprocess.run(['git', '-C', str(d), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
        out[lib] = f'{mod.__version__}' + (f'@{commit}' if commit else '')
    return out


def hybrid_id(choice, config, version):
    """ A hybrid's ID: a hash of its choice, the stop year, the model files (by contents, since one can be outside models/), the code and the pins """
    files = {m: hashlib.sha1((REPO / f).read_bytes()).hexdigest() for m, f in config.get('models', {}).items()}
    key = dict(choice=name(choice), stop=config.get('stop'), models=files, **version)
    return hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]


def save(sim, folder, cell):
    """ The run's tables: with Outputs, `<cell>.csv.gz` one row per year (common.outputs.yearly) and `<cell>.targets.json`
    (common.outputs.measure); without, the sim's results. Written to temporary names, then renamed """
    from common import outputs, targetfit
    o = sim.analyzers.get('outputs')
    df = outputs.yearly(o) if o is not None else sim.to_df()
    if o is not None:
        t = folder / f'{cell}.targets.json.tmp'
        t.write_text(json.dumps(outputs.measure(o.strata, targetfit.load_targets(), o.scale)))
        os.replace(t, folder / f'{cell}.targets.json')
    t = folder / f'{cell}.csv.gz.tmp'
    df.to_csv(t, index=False, compression='gzip')
    os.replace(t, folder / f'{cell}.csv.gz') # Last: its existence marks the run done


def start(choice, cell, config, folder):
    """ The sim of a run, ready to run. If D's model has `BURNIN` (a year after the sim's start), the run is resumed from
    its burn-in to then, `<folder>/<cell>.burnin.sim`, which is run and saved the first time (swap_design D1). Resumed
    equals continuous exactly (tests/test_resume.py), so this changes no result """
    sim = build(choice, cell, config)
    year = getattr(load_model(choice['D'], config.get('models', {}).get(choice['D'])), 'BURNIN', None) if 'D' in choice else None
    if year is None:
        return sim
    import starsim as ss
    f = Path(folder) / f'{cell["cell"]}.burnin.sim'
    if f.exists():
        return ss.load(f)
    sim.init()
    years = np.asarray(sim.t.yearvec, dtype=float)
    k = int(np.argmin(np.abs(years - year)))
    if k > 0 and k < len(years) - 1:
        sim.run(until=float(years[k - 1])) # Its last step is the one before `year`'s
        t = f.parent / f'{f.name}.tmp'
        sim.save(t)
        os.replace(t, f)
    return sim


def run_one(job):
    """ One run in a worker: build (or resume from the burn-in), run, save; returns its manifest row """
    config, choice, cell, hid, out, meta = job
    t0 = time.time()
    row = dict(hybrid=hid, name=name(choice), cell=cell['cell'], pids=json.dumps(cell['pids']), seed=cell['seed'], **meta)
    try:
        sim = start(choice, cell, config, out / hid)
        sim.run()
        save(sim, out / hid, cell['cell'])
        row['error'] = ''
    except Exception:
        row['error'] = traceback.format_exc().strip().splitlines()[-1]
        print(f'{name(choice)} {cell["cell"]}:\n{traceback.format_exc()}', file=sys.stderr)
    row['seconds'] = round(time.time() - t0, 1)
    return row


def run(config, out=OUT, cores=None):
    """ Run every hybrid in every cell not already stored under `out`; each finished run adds a row to `out/manifest.csv`.
    Returns the manifest rows of this call """
    out = Path(out)
    version = dict(code_version(), **pins())
    meta = dict(commit=git('rev-parse', '--short', 'HEAD'), dirty=version['dirty'], starsim=version['starsim'], stisim=version['stisim'])
    jobs = []
    for choice in hybrids(config):
        hid = hybrid_id(choice, config, version)
        (out / hid).mkdir(parents=True, exist_ok=True)
        (out / hid / 'hybrid.json').write_text(json.dumps(dict(name=name(choice), choice=choice, stop=config.get('stop'), models=config.get('models', {}), **version), indent=1))
        jobs += [(config, choice, c, hid, out, meta) for c in cells(config) if not (out / hid / f'{c["cell"]}.csv.gz').exists()]
    rows = []
    manifest = out / 'manifest.csv'
    if not jobs:
        return rows
    with Pool(cores or min(len(jobs), os.cpu_count())) as pool:
        for row in pool.imap_unordered(run_one, jobs):
            pd.DataFrame([row]).to_csv(manifest, mode='a', header=not manifest.exists(), index=False)
            rows.append(row)
    return rows


def load(config, out=OUT):
    """ The stored yearly tables of a config's hybrids and cells, as one DataFrame with `hybrid`, `name` and the cell's columns """
    out, version = Path(out), dict(code_version(), **pins())
    dfs = []
    for choice in hybrids(config):
        hid = hybrid_id(choice, config, version)
        for c in cells(config):
            f = out / hid / f'{c["cell"]}.csv.gz'
            if f.exists():
                dfs.append(pd.read_csv(f).assign(hybrid=hid, name=name(choice), cell=c['cell'], pair=c['pair'], seed=c['seed'], **{k: choice[k] for k in choice}))
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
