# Core question (Robyn, 9 Nov 2026)

**Decompose the difference between EMOD and Optima HIV in Malawi in Figure 3A of Bansi-Matharu et al. 2025** (Lancet Glob Health 13:e1006–19, doi:10.1016/S2214-109X(25)00121-4, PMC12500161; figure in `resources/bansi-matharu2025_fig3.jpg`): the source of new infections by the source partner's treatment status, and the total number of new infections, using the swap sequence and the common modelling architecture (`swap_design.md`). If EMOD proves too confusing, Optima (simpler) can carry more of the work.

## The quantities

| | EMOD | Optima HIV | Source |
|---|---|---|---|
| New infections, 2024 | 30,582 | 8,459 | Figure 3 caption |
| Source undiagnosed, 2024 | 33.3% | ~38% | text: 33.3–75.3% across 4 models; EMOD 33.3% (falls to 23.7% by 2040); Optima "stable at about 37.8%" |
| Source undiagnosed, 2020 | ~38% | ~46% | read off Figure 3A (approximate, to digitise) |
| Source diagnosed, never treated, 2020 | ~3% | ~30% | read off Figure 3A (approximate) |
| Source on ART, 2020 | ~41% | ~16% | read off Figure 3A (approximate) |
| Source interrupted ART, 2020 | ~18% | ~8% | read off Figure 3A (approximate); 2024 range across models 8.4–20.1% |

Year: Robyn's focus is 2020 (inside the years with treatment data, where Optima′ is verified, 2000–2020); the caption's totals are for 2024.

## What this changes

- **The paper's headline outputs** are new infections (total) and transmissions by the source's cascade stage (`common/transmission.py`, `common/outputs.py`) in 2020, for EMOD vs Optima. HIVsim becomes the common architecture (host) rather than one of the two ends.
- **The core path is EMOD′ ↔ Optima′.** `swap_design.md` §5 deferred it to January because Optima′ needs: agents, agent versions of its engine-reading pieces, the population label, `CascadeHost` on declared states, and Optima's cascade as a table. Those become the priority.
- **First check: do our references reproduce Figure 3A?**
  - EMOD: our rerun of the PR #1 files gives 50–55% of transmissions from people who interrupted ART in 2020/2024 (analysis 09), against ~18% in the published EMOD. That fits Daniel's view (PR #4) that the dropout trap (and the mortality factor) came in after the calibration. The fixed variant (D10, analysis 17) should be compared with Figure 3A.
  - Optima: Optima′ is verified against our Optima rerun for 2000–2020 (within 2%), but transmissions by source stage have not yet been computed from either.
- Optima's large "diagnosed, never treated" share around 2015–2023 (the yellow band in Figure 3A) is a structural feature to explain (answered below: it is people diagnosed but not yet linked to care, not people waiting in care).

## Do our two ends reproduce Figure 3A? (9 Nov, analyses on branch `starsim-rc3`)

| 2020 shares: undiagnosed / diagnosed never treated / on ART / interrupted | New infections | Source |
|---|---|---|
| **Published EMOD:** 38.5 / 3.7 / 40.6 / 17.1 | 30,582 (2024) | `common/fig3a_malawi.csv` (digitised) |
| EMOD′ with Daniel's PR #1 files (dropout trap, mortality bug copied): 32 / 3 / 14 / 51 | 162k | analysis 17 |
| **EMOD′ with Daniel's two fixes (`emod_fixed`), uncalibrated:** 35.6 / 3.4 / 48.1 / 12.9 | ~42k (2024, census-rescaled) | analysis 17 |
| EMOD binary with Daniel's two fixes (10 sets × 5 seeds): 36.7 / 4.4 / 45.0 / 13.9 | 41,300 (2024; 5–95% 26.5–59.4k; sets 30.6–60.5k) | analysis 18 |
| **Published Optima:** 45.2 / 28.2 / 17.6 / 9.0 | 8,459 (2024) | `common/fig3a_malawi.csv` |
| Our Optima (2022 project) and Optima′ (sexual): 42.3 / 40.7 / 11.4 / 5.7 | 29,191 (2024, projection with inputs frozen at 2020) | analysis 19 |

- **EMOD:** the EMOD binary with Daniel's dropout fix reproduces the published shares (every 2020 and 2024 share within 5 points, inside our 5–95% range; analysis 18), and EMOD′ with the same fixes agrees with the binary. Daniel's two fixes (PR #4) bring the split close to the published one without refitting, and fit the targets far better (score 108 → 7.5). This supports his view that the published model didn't have those bugs. The total is still ~40% high. Pending Daniel's answer on which files were calibrated.
- **Optima:** our project is a different calibration from MIHPSA's (shares differ from 2000 on; totals 2020 ~30k vs ~22k, and 2024 projections need inputs past 2020). We need MIHPSA's Optima Malawi project from Rowan. Optima's "diagnosed, never treated" band is people diagnosed but not yet linked to care (Optima's `linktocare` rises to 3.16 years by 2020).
- **The first swap path (HIVsim → EMOD′ copied-bug, analysis 16):** the cascade step carries almost all of the share difference; the HIV and network steps move the shares a lot in opposite directions and nearly cancel; demography multiplies infections without moving shares.

## First decomposition step: the cascade, both ways (analysis 24, 9 Nov)

Swapping the cascade between EMOD′-fixed and Optima′ (our two ends) moves the Figure 3A shares a long way toward the other model in both directions: about 0.5–0.75 of the gap in "diagnosed, never treated" and "on ART", and all of it (or more) in "interrupted"; roughly symmetric except "undiagnosed". But **new infections move opposite to the gap**: EMOD's cascade lowers infections in either host, while EMOD′-fixed has more infections than Optima′. So the cascade explains the direction of the share split, and the total-infection gap (Figure 3's 30k vs 8.5k) must come from the other components (HIV natural history and transmission, network, demography): the next steps of the path. Caveats: Optima's cascade on EMOD′'s people needs ART arm switching and PMTCT decisions (analysis 23's provisional decisions 5 and 7), and our Optima end is the 2022 project.
