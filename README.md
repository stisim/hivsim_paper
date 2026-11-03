# hivsim_paper

This repository holds the work for decomposing disagreement between HIV models. HIVsim, EMOD-HIV and Optima HIV, each fitted to Malawi, are re-expressed as interchangeable components in Starsim v4. Swapping those components one at a time traces how much of the difference in policy-relevant outputs each modelling choice accounts for.

## Docs

| File | What it is |
|---|---|
| [docs/vision.md](docs/vision.md) | The proposal: motivation, goals and methods |
| [docs/narrative_draft.md](docs/narrative_draft.md) | Draft introduction, with references |
| [docs/plan.md](docs/plan.md) | Weekly sprint plan to 15 Dec 2026, with the sprint 2 review |
| [docs/targets.md](docs/targets.md) | Calibration target set for HIVsim Malawi (P3) |
| [docs/optima_inventory.md](docs/optima_inventory.md) | Optima HIV Malawi: components, inputs, and its dependence on the time step |
| [docs/optima_prime_plan.md](docs/optima_prime_plan.md) | Optima′ build plan and acceptance |
| [experiments/](experiments/) | HIVsim Malawi calibration experiments, each with a README (before running) and a SUMMARY (after) |
| [docs/prompts.md](docs/prompts.md) | Agent prompt plan (P0–P5), human checkpoints, decisions, notes on Starsim v4, and a record of prompts run |

## Setup

There are two conda environments:

```bash
conda env create -f envs/starsim.yml  # HIVsim (models/hivsim) and Optima′ (models/optima)
conda env create -f envs/optima.yml   # Optima HIV 2.12.2, only for regenerating the references in models/optima/original/
```

- `starsim` pins Starsim to its first v4 release candidate (4.0.0rc1, on PyPI, Python 3.12 or later) and STIsim to a commit on `rs-art-split` (PR #607, which adds the ART split and HIV care). Change a pin on purpose and rerun the checks: `models/optima/acceptance.py` (about 3 minutes) and `models/hivsim/malawi.py`.
- Optima's references (`models/optima/original/*.obj`) are committed, so the Optima′ checks run in `starsim` alone.
- The Starsim wheel on PyPI comes with its kernels built. Installing Starsim from a git clone instead means compiling them (Numba and Rust); see Starsim's docs.
- HIVsim's Malawi data is in `models/hivsim/data/` (described in its README). `models/hivsim/make_data.py` makes the calibration targets and `n_art.csv` from UNAIDS estimates in `raw_data/`, which isn't committed.

## Layout

| Folder | What it holds |
|---|---|
| `models/<model>/` | One folder per model: the Starsim v4 version, with the original model's files in `original/` |
| `models/hivsim/` | HIVsim Malawi (`malawi.py`), its data (`data/`), and the scripts that make the data and targets and score runs against them |
| `models/optima/` | Optima′ (`optima_prime.py`) and its checks against Optima (`compare.py`, `acceptance.py`) |
| `models/optima/original/` | The Optima HIV Malawi project (`Malawi_20220912.prj`), the scripts that rerun it on Optima 2.12.2, and its reference outputs |
| `models/emod/original/` | Daniel's calibrated EMOD-HIV Malawi model: inputs, parameter sets, calibration targets, run scripts (see its README) |
| `experiments/NN_slug/` | HIVsim calibration experiments: README before running, SUMMARY after, outputs and figures |
| `spikes/` | Short technical tests, kept for the record |
| `docs/` | Proposal, plan, prompts, target set, Optima notes |
| `results/` | Committed results shared across the project (now `targets.csv`) |
| `resources/` | Source documents, such as survey reports |
| `envs/` | Conda environments |

## Conventions

- `raw_results/`: large outputs, including runs of the original models. Not committed (gitignored).
- `results/`: committed. Holds everything the plotting scripts need.
