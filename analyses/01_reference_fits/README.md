# 01: Reference models against the common targets

**Question:** how well does each reference model (the published or finalized version, before any translation to Starsim v4) match the common Malawi targets in `common/targets.csv`?

**What's run:**
- **EMOD:** Daniel's model, 10 of his 250 resampled parameter sets (ids 0, 25, ..., 225, spanning his fit scores), at 1% population, to 2025, using `models/emod/reference/run_emod.py --sets 0 25 50 75 100 125 150 175 200 225 --out sub10`. Then `models/emod/reference/emod_targets.py --out sub10` measures the targets and draws `figures/emod_*`.
- **Optima:** the 2022 Malawi project rerun on Optima HIV 2.12.2 (`models/optima/reference/run_optima.py`). So far only its totals plot (`figures/optima_totals.png`); not yet measured against the common targets.
- **HIVsim:** to come, once hivsim_malawi has a finalized parameter set (`models/hivsim/reference/`).

**Outputs:** `outputs/<model>_results.jsonl` (target values per run), `<model>_coverage_table.csv` and related, from `common/targetfit.py`.
