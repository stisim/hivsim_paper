# 01: Reference models against the common targets

**Question:** how well does each reference model (the published or finalized version, before any translation to Starsim v4) match the common Malawi targets in `common/targets.csv`?

**What's run:**
- **EMOD:** Daniel's model, all 250 of his resampled parameter sets, at 1% population, to 2025: `models/emod/reference/run_emod.py --cores 100 --out all250` (about an hour on 160 cores; about 12 minutes per run when the machine isn't loaded). Then `models/emod/reference/emod_targets.py --out all250` measures the targets and draws `figures/emod_*`. (A first pass used 10 sets, ids 0, 25, ..., 225.)
- **Optima:** the 2022 Malawi project (parset 20220912) on Optima HIV 2.12.2, dt = 0.2. It is deterministic, so there is one run. `models/optima/reference/export_ref.py` (env `optima`) saved its compartment-level outputs to the committed `optima_malawi_ref.obj`; `models/optima/reference/optima_targets.py` (env `hivsim_paper` or `starsim`) measures the targets from that file and draws `figures/optima_coverage_*`. `run_optima.py` draws `figures/optima_totals.png`.
- **HIVsim:** to come, once hivsim_malawi has a finalized parameter set (`models/hivsim/reference/`).

**Outputs:** `outputs/<model>_results.jsonl` (target values per run), `<model>_coverage_table.csv` and related, from `common/targetfit.py`.
