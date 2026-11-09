# Swap runner

Builds hybrids from each model's slot factories, runs them in cells in a process pool, stores one table per run with a manifest, and skips runs already stored (`docs/swap_design.md` 2.4, 3 and 4).

## Slot factories (what a model provides)

One file per model, by default `models/emod/v4/components.py` and `models/hivsim/components.py` (`runner.MODELS`; a config's `models:` overrides them). Its folder goes on `sys.path`, so it can import its model's modules flat.

```python
HOST = EmodHost          # the CascadeHost (with cascade_stage()) for when this model's V is in the hybrid; None if none
SCALE = 100              # people per agent, for counts (optional; default sim.pop_scale)
SIM = sti.Sim            # the sim class when this model's D is in the hybrid (optional; default ss.Sim)

def pars(pid): ...       # parameter set `pid`: whatever this model's factories read

def make_D(pars, ctx):   # one per slot the model fills: D, RN, V, C (and later R, N, T, H, S, M)
    return dict(n_agents=..., start=1960.5, dt=1/12, demographics=[EmodDemography(pars)])
```

A factory returns `ss.Sim` keywords. The runner calls factories in the canonical slot order (D, R, N, RN, V, T, H, S, C, M) and joins the module lists (`demographics`, `networks`, `diseases`, `interventions`, `connectors`, `analyzers`) in that order. Any other keyword (`n_agents`, `start`, `dt`, `total_pop`, ...) can be set by more than one slot only if they agree. It then sets `rand_seed` to the cell's seed and `stop` to the config's, and adds `common.outputs.Outputs(HOST of V's model)` as the last analyzer.

`ctx` has `choice` ({slot: model}), `models` ({model: its factory module}), `pars` ({model: its parameter set for the cell}), `seed`, `stop` and `cell`. A factory that has to fit another model's slot reads `ctx.choice` and does so itself. For example, EMOD's C runs the condensed cascade on `ctx.models[ctx.choice['V']].HOST`. When every slot comes from one model, nothing foreign is inserted. So identity holds by construction if the model's own `make_sim` is `build({slot: model, ...})` or calls the same factories in the same order.

For STIsim's HIV, `Outputs` reads sources from Starsim's infection log, so HIVsim's V factory should add `ss.infection_log()` to `analyzers`. It doesn't change the run.

`tests/emod_standin.py` is EMOD′'s factories as `harness.make_sim(cascade=True, condensed=True)` builds it. The tests use it until `models/emod/v4/components.py` exists, and then use that file.

## Configuration

A dict (or YAML read into one):

```python
config = {
    'from': 'hivsim', 'to': 'emod',
    'slots': ['C', 'V', 'RN', 'D'],   # sequential path: swapped from `from` to `to` in this order (5 hybrids)
    'lattice': False,                 # True: every subset (2^n hybrids, for Shapley)
    # or 'hybrids': [{'D': 'emod', 'RN': 'emod', 'V': 'emod', 'C': 'emod'}, ...] given outright
    'pairs': {'hivsim': [...], 'emod': [0, 25, 50]},  # one parameter-set ID per model per pair, paired 1:1
    'seeds': 3,                       # seeds per pair; cell `p{i}_s{s}` has seed 100*i + s, shared by every hybrid
    'stop': 2025,
    'models': {},                     # optional {model: factory file}, relative to the repo
}
```

Hybrids are named as in 2.4: `D=hivsim,RN=hivsim,V=emod,C=emod` (slots in canonical order).

## Running

```python
from common import swaps
swaps.run(config, cores=150)          # into raw_results/swaps/; returns this call's manifest rows
df = swaps.load(config)               # stored yearly tables of the config's hybrids, with name, slots, cell, pair, seed
sim = swaps.build(choice, swaps.cells(config)[0], config).run()   # one hybrid by hand
swaps.differences(sim_a, sim_b)       # [] if every result (and analyzer df/strata) is identical
```

For long runs, detach: `setsid nohup python run.py > log 2>&1 < /dev/null & disown`. Each worker uses one Numba thread (`NUMBA_NUM_THREADS=1`, set when the runner is imported).

## Store

Under `raw_results/swaps/` (gitignored):
- `<hybrid id>/hybrid.json`: the choice, stop, model files, code version and pins.
- `<hybrid id>/<cell>.csv.gz`: one row per year, from `common.outputs.yearly()`. Runs without `Outputs` store `sim.to_df()` instead. This file is written last, and its existence marks the run done.
- `<hybrid id>/<cell>.targets.json`: `{tid: value}` for `common/targets.csv` (`common.outputs.measure()`); NaN outside the run.
- `manifest.csv`: one row per finished run: hybrid id, name, cell, pids, seed, commit, dirty, Starsim and STIsim pins, seconds, and error (empty if it ran). A run that failed has a row but no table, so the next call runs it again.

The hybrid ID hashes the choice, the stop year, the model files, the git tree IDs of `models/` and `common/` (plus a hash of their uncommitted changes), and the Starsim and STIsim versions and commits. So a hybrid is rerun when the code it runs changes, but not for commits elsewhere (docs, analyses). The commit is recorded in the manifest.

## Tests

```bash
PYTHONPATH=<stisim rs-v4-art> conda run -n hp_rc3 python -m pytest common/swaps/tests -v   # EMOD′ ones take a few minutes
conda run -n hp_rc3 python -m pytest common/swaps/tests -k toy                             # seconds
```

- EMOD′ endpoint: `build(all EMOD′)` equals `harness.make_sim(cascade=True, condensed=True)`. That covers every result array, both reports' tables, and the outputs' counts and transmission log, to 2005.
- Trivial swap: C and RN from `emod_b` (EMOD′'s factories under another name) still equal EMOD′.
- Toy models (`tests/toy_a.py`, `toy_b.py`: SIS on a random network): hybrid lists, a swap for an alias changes nothing while a swap for the other model does, and pool, store, manifest and resume.
