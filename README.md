# hivsim_paper

This repository holds the work for decomposing disagreement between HIV models. HIVsim, EMOD-HIV and Optima HIV, each fitted to Malawi, are re-expressed as interchangeable components in Starsim v4. Swapping those components one at a time traces how much of the difference in policy-relevant outputs each modelling choice accounts for.

## Docs

| File | What it is |
|---|---|
| [docs/vision.md](docs/vision.md) | The proposal: motivation, goals and methods |
| [docs/narrative_draft.md](docs/narrative_draft.md) | Draft introduction, with references |
| [docs/plan.md](docs/plan.md) | Weekly sprint plan to 15 Dec 2026 |
| [docs/prompts.md](docs/prompts.md) | Agent prompt plan (P0–P5), human checkpoints, decisions, notes on Starsim v4, and a record of prompts run |

## Setup

There are two conda environments:

```bash
conda env create -f envs/starsim.yml  # HIVsim (models/hivsim) and Optima′ (models/optima)
conda env create -f envs/optima.yml   # Optima HIV 2.12.2, only for regenerating the references in optima_ref/
```

- `starsim` pins Starsim and STIsim to commits, since both are changing: Starsim on `v4`, and STIsim on `rs-art-split` (PR #607, which adds the ART split and HIV care). Change a pin on purpose and rerun the checks: `models/optima/acceptance.py` (about 3 minutes) and `models/hivsim/malawi.py`.
- `sciris` is pinned to 3.6.2, since 3.7.0 breaks importing Starsim.
- Optima's references (`optima_ref/*.obj`) are committed, so the Optima′ checks run in `starsim` alone.
- Starsim's Rust kernels are optional. Without them Starsim uses Numba, which gives the same numbers but compiles the first time it runs. To build them, see `starsim/solvers/rust/README.md` in a source checkout.
- HIVsim's Malawi data isn't committed yet. `models/hivsim/data/README.md` lists the files it needs.

## Conventions

- `raw_results/`: large outputs, including runs of the original models. Not committed (gitignored).
- `results/`: committed. Holds everything the plotting scripts need.
