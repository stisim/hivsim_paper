# hivsim_paper

This repository holds the work for decomposing disagreement between HIV models. HIVsim, EMOD-HIV and Optima HIV, each fitted to Malawi, are re-expressed as interchangeable components in Starsim v4. Swapping those components one at a time traces how much of the difference in policy-relevant outputs each modelling choice accounts for.

## Docs

| File | What it is |
|---|---|
| [docs/vision.md](docs/vision.md) | The proposal: motivation, goals and methods |
| [docs/narrative_draft.md](docs/narrative_draft.md) | Draft introduction, with references |
| [docs/plan.md](docs/plan.md) | Weekly sprint plan to 15 Dec 2026 |
| [docs/targets.md](docs/targets.md) | Calibration target set for HIVsim Malawi (P3) |
| [docs/prompts.md](docs/prompts.md) | Agent prompt plan (P0–P5), human checkpoints, decisions, notes on Starsim v4, and a record of prompts run |

## Setup

There are two conda environments:

```bash
conda env create -f envs/starsim.yml  # HIVsim (models/hivsim) and Optima′ (models/optima)
conda env create -f envs/optima.yml   # Optima HIV 2.12.2, only for regenerating the references in optima_ref/
```

- `starsim` pins Starsim to its first v4 release candidate (4.0.0rc1, on PyPI, Python 3.12 or later) and STIsim to a commit on `rs-art-split` (PR #607, which adds the ART split and HIV care). Change a pin on purpose and rerun the checks: `models/optima/acceptance.py` (about 3 minutes) and `models/hivsim/malawi.py`.
- Optima's references (`optima_ref/*.obj`) are committed, so the Optima′ checks run in `starsim` alone.
- The Starsim wheel on PyPI comes with its kernels built. Installing Starsim from a git clone instead means compiling them (Numba and Rust); see Starsim's docs.
- HIVsim's Malawi data is in `models/hivsim/data/` (described in its README). `models/hivsim/make_data.py` makes the calibration targets and `n_art.csv` from UNAIDS estimates in `raw_data/`, which isn't committed.

## Conventions

- `raw_results/`: large outputs, including runs of the original models. Not committed (gitignored).
- `results/`: committed. Holds everything the plotting scripts need.
