# 16: The first swap path, HIVsim → EMOD′ (uncompensated)

Step 11 of the first milestone (`docs/swap_design.md` §5): run the path and decompose the endpoints' difference into its steps.

**Question:** along HIVsim → EMOD′, swapping the coarse slots in the order C, V, RN, D (`swap_design` D4), how much of the difference between the two models' §2.2 outputs comes with each slot, before any refit (tier 1, uncompensated)? And how much do these Δs depend on conditioning on the epidemic taking off (`swap_design` D11, provisional)?

**Hybrids** (D, RN, V, C; H = HIVsim, E = EMOD′; built by the runner, `common/swaps`, as in analysis 15):

| # | Hybrid | Step into it | Scale |
|---|---|---|---|
| 0 | HHHH | (HIVsim) | 10,000 agents |
| 1 | HHHE | C: EMOD′'s condensed cascade on HIVsim | 10,000 agents |
| 2 | HHEE | V: EMOD′'s HIV on HIVsim's network and demography | 10,000 agents |
| 3 | HEEE | RN: EMOD′'s network and risk on HIVsim's demography | 10,000 agents |
| 4 | EEEE | D: EMOD′'s demography (EMOD′) | 1% (resumed from its 1960.5–1985 burn-in) |

**Code and parameters:** STIsim `rs-v4-art-mortality` (PR #613, ebbc8e7; mortality options off) and Starsim `v4-rs-nov09` (#1573, f1273734), env `hp_rc3`. HIVsim: `models/hivsim/reference/pars.json` (hivsim_malawi experiment 03's best run, provisional, and not recalibrated on this STIsim). EMOD′: parameter set 0. One pair (hivsim 0, emod 0).

**Seeds:** 200 cells (`p0_s1` … `p0_s200`, seed = s), shared by every hybrid, so each step's Δ is paired by seed. Sizing: analysis 14's rerun with the engine and condom fixes (`/tmp/agents-nov09/P/compare.txt`) needs 30–40 seeds for Δ new infections at 10% relative precision, ~125 for the share of deaths on ART, and 160–360 for Δ PLHIV; the V and RN steps have low take-off (2–3 of 5 in analysis 15), which leaves fewer seeds for the conditional Δs. Runs are cheap (~10 s for the 10,000-agent hybrids, ~80 s for EEEE), so 200 covers all but the PLHIV extreme with conditioning.

**Headline** (added 9 Nov, `docs/core_question.md`: the paper's core quantity is Bansi-Matharu et al. 2025 Figure 3A): 2020 new infections (total) and the share of sexual transmissions by the source's treatment status in Figure 3A's four categories: undiagnosed, diagnosed never treated, on ART (`common/outputs.py`'s on ART < 6 months + suppressed + not suppressed, as `tx_on_art`), interrupted. The Δ table leads with these; the other outputs are kept but not examined further.

**Outputs** (`common.outputs.yearly`, every §2.2 output but the ranking), for 2010, 2015, 2020 and 2024: ART coverage; PLHIV; new infections (all ages); HIV incidence 15–49 per 1000 person-years, women and men; HIV deaths; share of HIV deaths on ART; share of sexual transmissions by the source's cascade stage (undiagnosed, diagnosed never on ART, on ART < 6 months, on ART suppressed, on ART not suppressed, interrupted). Counts (PLHIV, new infections, HIV deaths) are scaled to people by the D's people per agent (HIVsim: `sim.pop_scale`; EMOD′: 100). Stocks at the year's end; flows over the calendar year.

**Take-off:** a run takes off if prevalence (all ages) at the end of 2020 is at least 0.5%. Reported with the distribution of PLHIV in 2020, and against analysis 15's rule (PLHIV > 0), to check that the threshold sits in the gap of a bimodal distribution.

**Analysis** (`run.py`):
- `runs`: the 5 hybrids × 200 cells through `swaps.run`, stop 2025, into `raw_results/swaps/` (gitignored).
- `identity`: endpoint identity at the full horizon (to 2025), seed 1, in separate processes: `build(HHHH)` against HIVsim's `swap.make_sim('own')` with the outputs analyzer, and EEEE (resumed from its stored burn-in, as in the runs) against EMOD′'s `harness.make_sim(cascade=True, condensed=True)`; every result array and the analyzers' tables (`swaps.differences`), and the stored yearly table against the rerun.
- `analyse`: per step (each pair of consecutive hybrids) and output and year, the paired Δ (swap − base) with a 95% t-interval over cells: (a) over all cells where the output is defined in both arms (shares are undefined, 0/0, without infections or deaths), (b) over cells where both arms took off, and the difference (b) − (a). (c) Take-off probability per hybrid (Wilson 95% interval), and per step the change in it. Also each hybrid's mean (all, and taken-off only). Tables in `outputs/`.
- Figures (`figures/`): per output, a waterfall along the path for each year: HHHH's mean, each step's Δ (all runs, and conditional, side by side) cumulated, ending at EEEE's mean, with the Δs' 95% intervals. The conditional Δs don't add up to the conditional endpoint difference, since each step keeps different cells.

**Not interpreted as findings:** HIVsim's parameters are provisional and uncalibrated on this STIsim, there is no refit (the compensated, tier 2 results are the headline later, `swap_design` §3), and EMOD′ is one parameter set. The Δs describe these two builds, not the models' published fits.

**What counts as done:** every run finishes without error; both endpoint identities hold exactly at the full horizon; every output has a Δ with an interval for every step and year, both unconditional and conditional (or the number of cells is too small to give one, and that is stated); take-off probabilities with intervals; waterfalls for every output; a SUMMARY with the Δ table and the effect of conditioning.

## Addendum: the path to `emod_fixed` (added 9 Nov, before running)

**Question:** with the EMOD′ end that matches the published figure, `emod_fixed` (analysis 17: Daniel's two fixes, `mortality` and `reenrol`; `models/emod/v4/components_fixed.py`), which steps of the path change, against the copied-bug EMOD′ above?

**Method:** as above (`run.py <cmd> --fixed`): config `to: emod_fixed`, `models: {emod_fixed: models/emod/v4/components_fixed.py}`, pair (hivsim 0, emod_fixed 0), the same 200 cells, order C, V, RN, D, from the `starsim-rc3` head. The fixes travel with the slots' parameter sets: `reenrol` with C (from HHHE on), `mortality` with D (EEEE only). So V and RN should be close to analysis 16's apart from what they inherit from C. Endpoint identity: EEEE-fixed (resumed from its burn-in) against `harness.make_sim(0, cascade=True, condensed=True, fixes=('mortality', 'reenrol'))`, and HHHH as before. Outputs and figures as above, with a `_fixed` suffix. Comparison: the 2020 headline (new infections and the four Fig 3A shares; all runs and conditional; take-off) for each step, against analysis 16's.

**What counts as done:** no run errors; both identities exact; the headline table with its comparison to analysis 16.
