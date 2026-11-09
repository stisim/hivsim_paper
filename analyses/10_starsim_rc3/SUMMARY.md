# 10: EMOD′ and Optima′ on Starsim 4.0.0rc3 — summary

**Answer:** apart from rc3's intended changes, nothing changes. EMOD′ gives bit-for-bit the same runs as on ebd52c1a. In Optima′, stocks and new infections are identical, and HIV deaths are identical in every year except the extra final step (open item below).

**EMOD′ (agents, 90 runs: compiled, condensed and interpreted, 10 sets × 3 seeds each)**
- All 30/30 runs of each variant are identical to their ebd52c1a runs: z = 0 at all 345 target points (`z_summary.csv`, `targets_vs_ebd52c1a_by_kind.csv`). That is better than the distribution-level agreement the README expected. EMOD′'s draws come from its own code paths, not from rc3's keyed declared streams, so the seeds carry over.
- The first rc3 run was off: 8 of 345 points were beyond 2 SE, with z up to 11.4 (`first_run_*.csv`). Every report came a month late, because the first step is now `ti == 0`. After the harness's report steps (`Snapshots`, `EmodReport`) moved down by one index (1f847d5), the runs were identical.
- 06's comparison (condensed vs compiled) gives 343/345 points within 2 SE on rc3, z mean 0.05, as on ebd52c1a. 04's comparison (EMOD′ vs EMOD) gives the same ratios to EMOD as before (`emod_ratio_by_kind.csv`).
- Run time: similar or a little faster. Paired runs at the same load took 75 s on rc3 vs 87 s on ebd52c1a for condensed, and 101 s vs 107 s for compiled (`time_paired.csv`; 1 run each, so indicative only).

**Optima′ (deterministic ODE)**
- People by state, population and year, and new infections, are identical, to within 1e-7 relative (`optima_arrays_{ebd52c1a,rc3}.npz`). This required `compare.py` to read each flow at its own row (rc3: a flow at row `ti` is over the step starting there), not at the next row.
- HIV deaths are identical in every year up to 2029.
- **Open:** in the last year, 2030, deaths are 1,230 vs 6,160. In rc3's extra final step (from the stop), every `hivdeath_*` flow is exactly 0.2× (dt×) its value in the step before, while new infections in that step are normal. A minimal declared SIS with `ss.death(I)` on rc3, stratified or not, gives a normal final step, so the cause is in how Optima′ is built or in something it uses that the minimal model doesn't. No target reads 2030. Look into it before any 2030 output is used, and if it turns out to be a Starsim bug, fix it in Starsim.
- Analysis 08's three Optima′ runs, compared with 08's outputs: stock targets are identical, and flow targets (new infections, AIDS deaths) differ by under 0.7% after 1990 (4.5% at 1990.5, the first year). That is the expected one-step shift: 08's `window` averages the rows in [y, y+0.2), which now hold the steps that start in that window rather than those that end in it. The rc3 reading is the right one, so 08 needs no change.

**Migration changes (no effect on results):** `update_results()` is split into `record_stocks()`/`record_flows()`. A delivery's trigger is now in `when=`, and `deliver=` is gone. `ss.ppl` → `ss.person`, `ss.null` → `ss.unborn`, `mode=` → `interpolation=`, `ss.transfer(dict(pop=...))`, `sim.pop_scale`, `aging=`, and `on_death()`. The env pins Starsim 4.0.0rc3. STIsim stays pinned where it was until its port (branch `rs-v4-art`) is merged.

**Next:** merge into main together with the STIsim port, so that HIVsim moves to rc3 at the same time.
