# 16: The first swap path, HIVsim → EMOD′ (uncompensated) — summary

**Answer:** the path runs end to end (1000 of 1000 runs, 200 cells × 5 hybrids, no errors), and both endpoints are exact. For the core quantity (2020 transmissions by the source's treatment status, Bansi-Matharu 2025 Fig 3A's categories), almost all of the HIVsim–EMOD′ difference in the shares comes with the cascade (C). V and RN move the shares by large amounts in opposite directions, and they nearly cancel. D moves no share, but it multiplies new infections. Conditioning on take-off barely changes the share Δs (at most 0.016). It changes the count Δs a lot: for the V step it flips the sign of Δ PLHIV. This is tier 1 (no refit), with HIVsim's provisional parameters (exp 03's best set, not recalibrated on this STIsim) and EMOD′ set 0 with the copied bugs. It describes these two builds, not the published models.

**Setup:** STIsim `rs-v4-art-mortality` (ebbc8e7, mortality options off) and Starsim `v4-rs-nov09` (f1273734), commit 45c2095. Cells `p0_s1`–`p0_s200`, shared by every hybrid. Counts are in people: HIVsim's D is 762 people per agent (10,000 agents), EMOD′'s is 100. Take-off means prevalence ≥ 0.5% at the end of 2020. Hybrids are named D, RN, V, C (H = HIVsim, E = EMOD′). HIVsim now reads ART's states on `sti.ART` (STIsim #612's split); STIsim pinned in `envs/starsim.yml`; HIVsim Malawi prevalence 15–49 in 2020: 8.1% (7.6–8.7%, 5 seeds).

**Endpoint identity** (`outputs/identity.csv`, seed 1, to 2025): `build(HHHH)` equals `swap.make_sim('own')` with the outputs analyzer, and EEEE resumed from its stored burn-in equals `harness.make_sim(cascade=True, condensed=True)`. In both, 0 arrays differ and the stored yearly table is reproduced.

**Take-off** (`outputs/takeoff.csv`, `takeoff_deltas.csv`), 95% Wilson intervals:

| Hybrid | Took off | p (95% CI) |
|---|---|---|
| HHHH | 200/200 | 1.00 (0.98–1) |
| HHHE | 200/200 | 1.00 (0.98–1) |
| HHEE | 144/200 | 0.72 (0.65–0.78) |
| HEEE | 164/200 | 0.82 (0.76–0.87) |
| EEEE | 200/200 | 1.00 (0.98–1) |

Paired Δp per step: C 0; V −0.28 (−0.34, −0.22); RN +0.10 (+0.04, +0.16); D +0.18 (+0.13, +0.23). Both arms take off in 200, 144, 133 and 164 cells for the four steps. The threshold sits in a gap: in runs that didn't take off, 2020 prevalence is ≤ 0.05% (HHEE) or ≤ 0.14% (HEEE); in runs that did, ≥ 5.3% (HHEE) or ≥ 0.59% (HEEE). Analysis 15's rule (PLHIV > 0) would count 3 more runs for HHEE and 8 more for HEEE as taken off.

**Headline: 2020, paired Δ (95% CI)**, (a) all cells where the output is defined in both arms; (b) both arms took off (`outputs/headline_2020.csv`)

| Output | HHHH | Step C | Step V | Step RN | Step D | EEEE |
|---|---|---|---|---|---|---|
| New infections (a) | 25,300 | +36,200 (34,800, 37,600) | +2,800 (−3,500, 9,000) | −24,800 (−30,500, −19,000) | +123,800 (119,600, 127,900) | 163,200 |
| New infections (b) | | same (n 200) | +27,900 (24,300, 31,600) n 144 | −40,000 (−44,800, −35,200) n 133 | +115,300 (111,400, 119,100) n 164 | |
| Undiagnosed (a) | 0.723 | −0.383 (−0.398, −0.369) | +0.236 (0.214, 0.257) n 145 | −0.242 (−0.270, −0.214) n 137 | −0.011 (−0.030, 0.008) n 167 | 0.319 |
| Undiagnosed (b) | | same | +0.240 (0.220, 0.259) | −0.255 (−0.276, −0.234) | −0.005 (−0.020, 0.010) | |
| Diagnosed, never treated (a) | 0.000 | +0.027 (0.024, 0.030) | +0.002 (−0.004, 0.007) | +0.009 (−0.006, 0.024) | −0.008 (−0.020, 0.005) | 0.028 |
| Diagnosed, never treated (b) | | same | +0.002 (−0.004, 0.007) | +0.002 (−0.004, 0.009) | −0.002 (−0.007, 0.002) | |
| On ART (a) | 0.264 | −0.153 (−0.167, −0.139) | +0.017 (0.007, 0.027) | +0.013 (0.001, 0.024) | +0.004 (−0.004, 0.013) | 0.144 |
| On ART (b) | | same | +0.018 (0.008, 0.028) | +0.016 (0.005, 0.027) | +0.002 (−0.007, 0.010) | |
| Interrupted (a) | 0.013 | +0.509 (0.498, 0.520) | −0.254 (−0.275, −0.233) | +0.220 (0.195, 0.245) | +0.014 (−0.003, 0.032) | 0.509 |
| Interrupted (b) | | same | −0.259 (−0.278, −0.241) | +0.237 (0.217, 0.256) | +0.005 (−0.009, 0.020) | |

"On ART" sums `common/outputs.py`'s on ART < 6 months, suppressed and not suppressed. The shares are of sexual transmissions. They are undefined (0/0) in runs with no transmissions, so (a) for the shares already leaves out most runs that didn't take off (n 145, 137, 167 against 144, 133, 164). That is why conditioning changes them so little.

**Reading (descriptive only):**
- **C carries the share difference.** EMOD′'s condensed cascade on HIVsim moves the shares from HIVsim's (72% undiagnosed, 26% on ART, 1% interrupted) almost to EEEE's (32 / 14 / 51%) in one step. Interrupted becomes 52%, which is analysis 09's dropout-trap level (published EMOD ~18%; analysis 17 shows the trap is a copied bug). C also more than doubles new infections, and it halves ART coverage, from 0.84 to 0.43.
- **V and RN cancel in the shares.** EMOD′'s HIV on HIVsim's network and demography moves undiagnosed by +0.24 and interrupted by −0.25. EMOD′'s network then reverses this (−0.24, +0.22). Their sum is within ~0.02 of zero. A decomposition that showed only the net of the two would miss two large, opposite interactions with the network.
- **D moves no share** (all four CIs include 0). It multiplies new infections by ~4 (HEEE 39,400 → EEEE 163,200). Incidence 15–49 in women rises by 5.7 per 1000, so this is more than the larger population (EMOD′'s copied mortality bug doubles the population; analysis 17).
- **Conditioning (D11):** for counts it matters a lot. Δ new infections at V is +2,800 (CI includes 0) over all runs and +27,900 conditional. At RN it is −24,800 against −40,000. Δ PLHIV 2020 at V flips sign: −151,000 (−226,000, −77,000) over all runs, +178,000 (165,000, 191,000) conditional. For shares, conditioning changes nothing that matters (|b − a| ≤ 0.016). So the count Δs on V and RN are mostly take-off, with Δp −0.28 and +0.10 as their own outputs. Conditional count Δs on those steps compare different subsets of cells along the path, so they don't add up to the endpoint difference.

**Other outputs** (every §2.2 output, 2010/2015/2020/2024: `outputs/deltas.csv`, `means.csv`, `deltas_2020_2024.csv`; waterfalls in `figures/waterfall_<output>.png`): kept, not examined further (priority of 9 Nov). In 2020, the share of HIV deaths on ART goes 0.00 → 0.22, all at V (+0.23): HIVsim's on-ART mortality (analysis 13) against EMOD's. HIV deaths rise at C (+8,700), V (+26,900) and D (+81,300).

**Caveats:**
- The EMOD′ endpoint (set 0, Daniel's PR #1 files, with the dropout trap and the mortality bug) gives 163,000 new infections and 2.2 M PLHIV in 2020, far above Fig 3's EMOD (30,582 infections in 2024). The fixed variant (`emod_fixed`, analysis 17) is the endpoint to compare with Fig 3A.
- HIVsim's parameters are provisional and uncalibrated on this STIsim. There is no refit (tier 2 is the headline later). The Δs are of these builds.
- C on HIVsim keeps HIVsim's D, so its ART stock is no longer pinned to Malawi's numbers on ART (HIVsim's own ART is a stock target). Part of C's Δ is therefore "flow cascade against stock target", not "EMOD against HIVsim".
- These runs are at 45c2095. Analysis 17's EMOD′ fix switches were merged into `starsim-rc3` after them, so rerunning from the branch head gives new hybrid IDs and reruns.

**Friction:** STIsim `stisim/interventions/hiv_interventions.py:577` warns "ART intervention added without an HIV testing intervention" in every hybrid with EMOD's C, where the condensed cascade does the testing through `CascadeHost`. Proposed: skip the warning when HIV has another diagnosing module, or add an `ART(warn_no_test=False)` option.

**Open:** a seeding step (S) or more seeding, so that take-off is near certain on HIVsim's demography (D11, still provisional). The EMOD′-fixed endpoint. Tier 2 refit.

## Addendum: the path to `emod_fixed` (Daniel's two fixes)

**Answer:** with the EMOD′ end that roughly reproduces the published Figure 3A (analysis 17), the path's story changes. C no longer carries the share difference by itself. The dropout trap was what made C look decisive: its Δ interrupted falls from +0.51 to +0.14, and Δ on ART changes sign (−0.15 → +0.10). On ART now comes from C (+0.10) and RN (+0.21), against V (−0.11). Interrupted is +0.14 at C, −0.11 at V, +0.07 at RN. The mortality fix removes almost all of D's effect on new infections (+123.8k → +17.5k). Take-off is unchanged (same cells; the fixes don't touch seeding). Conditioning still matters only for counts. Same caveats as above: tier 1, HIVsim provisional, one parameter set.

**Setup:** config `to: emod_fixed`, `models: {emod_fixed: models/emod/v4/components_fixed.py}`, the same 200 cells, from `starsim-rc3` at 5526730 (`run.py <cmd> --fixed`, outputs and figures `*_fixed`). `reenrol` travels with C (from HHHE on), `mortality` with D (EEEE only). 1000 of 1000 runs finished with no errors.

**Endpoint identity** (`outputs/identity_fixed.csv`, seed 1, to 2025): EEEE-fixed resumed from its burn-in equals `harness.make_sim(0, cascade=True, condensed=True, fixes=('mortality', 'reenrol'))`, with 0 arrays differing and the stored table reproduced. HHHH is exact as before, and it is identical to analysis 16's HHHH. EEEE-fixed in 2020: 40,300 new infections; shares undiagnosed / diagnosed never treated / on ART / interrupted 0.359 / 0.039 / 0.469 / 0.134. Analysis 17 set 0 had 38k and 0.356 / 0.034 / 0.481 / 0.129. The published EMOD has 0.385 / 0.037 / 0.406 / 0.171.

**Take-off** (`takeoff_fixed.csv`): identical to analysis 16. HHHH, HHHE and EEEE 200/200; HHEE 144/200 (0.72, 0.65–0.78); HEEE 164/200 (0.82, 0.76–0.87). Δp per step: C 0, V −0.28, RN +0.10, D +0.18.

**Headline, 2020, paired Δ (95% CI)**, (a) all runs / (b) both arms took off, with analysis 16's Δ (copied bugs) in brackets (`outputs/headline_2020_fixed.csv`)

| Output | HHHH | C | V | RN | D | EEEE-fixed |
|---|---|---|---|---|---|---|
| New infections | 25,300 | +4,500 (3,400, 5,700) [+36,200] | (a) +20,200 (14,900, 25,500) [+2,800]; (b) +39,600 (35,500, 43,800) [+27,900] | (a) −27,200 (−31,900, −22,400) [−24,800]; (b) −40,500 (−45,000, −36,100) [−40,000] | (a) +17,500 (15,100, 19,800) [+123,800]; (b) +12,300 (10,200, 14,500) [+115,300] | 40,300 |
| Undiagnosed | 0.723 | −0.286 (−0.304, −0.268) [−0.383] | +0.220 (0.196, 0.244) [+0.236] | (a) −0.281; (b) −0.288 (−0.317, −0.258) [−0.242 / −0.255] | (a) −0.023 (−0.045, −0.002); (b) −0.018 (−0.038, 0.002) [−0.011] | 0.359 |
| Diagnosed, never treated | 0.000 | +0.047 (0.041, 0.052) [+0.027] | −0.001 (−0.010, 0.009) [+0.002] | −0.005 (−0.015, 0.005) [+0.009] | −0.001 (ns) [−0.008] | 0.039 |
| On ART | 0.264 | **+0.098** (0.081, 0.114) [−0.153] | **−0.110** (−0.129, −0.090) [+0.017] | **+0.211** (0.184, 0.238); (b) +0.218 [+0.013] | +0.013 (ns) [+0.004] | 0.469 |
| Interrupted | 0.013 | **+0.142** (0.132, 0.152) [+0.509] | **−0.109** (−0.121, −0.097) [−0.254] | **+0.074** (0.061, 0.087) [+0.220] | +0.012 (0.001, 0.022) [+0.014] | 0.134 |

Conditioning changes the share Δs by at most 0.007 here. For the shares, V's (a) and (b) are the same cells (n 144).

**Which steps change (fixed against copied bugs):**
- **C:** much smaller. Δ new infections +36.2k → +4.5k, Δ interrupted +0.51 → +0.14, and Δ on ART −0.15 → +0.10 (sign change). The dropout trap made EMOD's cascade look like the whole share difference. With re-enrolment, people who drop out come back onto ART, so EMOD's C on HIVsim leaves HHHE (0.44 / 0.05 / 0.36 / 0.16) part way to EEEE-fixed.
- **V and RN:** they still largely cancel each other in undiagnosed (+0.22, −0.28). They no longer cancel in on ART (−0.11, +0.21): EMOD′'s network on top of EMOD′'s HIV raises the on-ART share to EEEE-fixed's 0.46. In the copied-bug path, interrupted soaked this up. V's Δ new infections over all runs is now clearly positive (+20.2k against +2.8k, not significant). RN's is unchanged.
- **D:** shares still essentially unchanged (|Δ| ≤ 0.023). The mortality fix removes almost all of D's effect on new infections (+123.8k → +17.5k), as analysis 17 found for the population (2× WPP with the bug).
- **Take-off:** unchanged, cell for cell.

**Outputs:** `outputs/*_fixed.csv` (deltas, means, takeoff, takeoff_deltas, headline_2020, deltas_2020_2024, prev_2020_by_cell, identity, run_times), `figures/waterfall_<output>_fixed.png`.
