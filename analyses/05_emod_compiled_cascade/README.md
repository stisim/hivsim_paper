# 05: EMOD′ with its campaign compiled into Starsim deliveries

**Question:** does EMOD′ with Daniel's campaign compiled into native Starsim deliveries (`models/emod/v4/compile_campaign.py`) give the same results as EMOD′ with the campaign interpreted literally (`cascade.py`), and what does it cost to run?

**What's run:**
- Starsim v4 head with PR #1538 (branch `v4-rs-nov04`: cycles of triggers through `delay=`), from a separate checkout on `PYTHONPATH`; both variants run on it.
- `harness.make_sim(pid, seed, cascade=True)` (literal) and `harness.make_sim(pid, seed, cascade=True, compiled=True)` (compiled), for the 10 parameter sets of analysis 04 (seeded from their run numbers), 3 seeds each: 60 runs.
- Measured as in analysis 04: the common targets (`emod_targets.measure()` on each run's `EmodReport`), the cascade over time, InterventionStatus of people with HIV state by state, and infections at ages 0-4. Also the run time of each run, and the number of objects the compiler makes.
- `run.py` runs and writes `outputs/`; `run.py --plot` draws `figures/`.

**How runs compare:** the two variants share their seeds, and nothing differs before the cascade's first event (1990). After that the literal cascade draws from its own random numbers and the compiled one from Starsim's distributions, so runs with the same seed can only be compared in distribution: as medians and ranges over the 30 runs of each.
