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
| [analyses/](analyses/) | Every result, numbered, each with a README (before running) and a SUMMARY (after) |
| [docs/prompts.md](docs/prompts.md) | Agent prompt plan (P0–P5), human checkpoints, decisions, notes on Starsim v4, and a record of prompts run |

## Setup

There are three conda environments:

```bash
conda env create -f envs/starsim.yml  # Optima′ (models/optima/v4), HIVsim, and the shared scoring and plots
conda env create -f envs/optima.yml   # Optima HIV 2.12.2, only for regenerating the references in models/optima/reference/
```

The third, `emod`, runs the EMOD-HIV binary; its setup is in `models/emod/README.md`.

- `starsim` pins Starsim to its first v4 release candidate (4.0.0rc1, on PyPI, Python 3.12 or later) and STIsim to a commit on `rs-art-split` (PR #607, which adds the ART split and HIV care). Change a pin on purpose and rerun the checks, such as `models/optima/v4/acceptance.py` (about 3 minutes).
- Optima's references (`models/optima/reference/*.obj`) are committed, so the Optima′ checks run in `starsim` alone.
- The Starsim wheel on PyPI comes with its kernels built. Installing Starsim from a git clone instead means compiling them (Numba and Rust); see Starsim's docs.

## Layout

Each model goes through three stages: its **reference** version (published or finalized, never tuned here), its **v4** translation into Starsim v4 components, and **swaps** of those components between models. Code and inputs live in `models/`; every figure and result lives in `analyses/`.

| Folder | What it holds |
|---|---|
| `models/<model>/reference/` | The reference model and the scripts that rerun it |
| `models/<model>/v4/` | The Starsim v4 translation (Optima′ now; EMOD′ to come). HIVsim is already v4 |
| `models/emod/reference/source/` | Daniel's calibrated EMOD-HIV Malawi model, unchanged (see `models/emod/README.md`) |
| `models/optima/reference/` | The Optima HIV Malawi project (`Malawi_20220912.prj`), the scripts that rerun it on Optima 2.12.2, and its reference outputs |
| `models/hivsim/reference/` | HIVsim's finalized parameter set, once [hivsim_malawi](https://github.com/stisim/hivsim_malawi) has one; HIVsim is calibrated there |
| `common/` | Shared by every model: the Malawi targets (`targets.csv`, described in `docs/targets.md`) and `targetfit.py`, which scores runs against them and draws the standard figures |
| `analyses/NN_slug/` | One per result: README before running, SUMMARY after, `outputs/`, `figures/` (files prefixed by model) |
| `spikes/` | Short technical tests, kept for the record |
| `docs/` | Proposal, plan, prompts, target set, Optima notes |
| `resources/` | Source documents, such as survey reports |
| `envs/` | Conda environments |

`raw_results/` holds large outputs, including runs of the reference models. Not committed (gitignored). Anything a figure needs is copied into the analysis's `outputs/`.
