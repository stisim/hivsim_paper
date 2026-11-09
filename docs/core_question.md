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
- Optima's large "diagnosed, never treated" share around 2015–2023 (the yellow band in Figure 3A) is a structural feature to explain: in Optima, diagnosed people wait in a care state before ART; in EMOD, linkage is fast or people are lost.
