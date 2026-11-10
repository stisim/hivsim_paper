# Sprint 3 readout (team check-in, 10 Nov 2026)

*Robyn, Daniel, Cliff. Covers sprints 1–3 (20 Oct–9 Nov). Plan from sprint 4: `docs/plan.md` on `main`. Analyses are on branch `starsim-rc3`.*

## Headlines

1. **The question is now sharper: decompose Figure 3A of Bansi-Matharu et al. 2025.** Why do EMOD and Optima HIV disagree, in Malawi, on who new infections come from (undiagnosed / diagnosed never treated / on ART / interrupted ART), and on how many there are (30.6k vs 8.5k in 2024)? `core_question.md` on `main`.
2. **All three models are rebuilt in Starsim v4 and verified against their originals.** EMOD′ reproduces EMOD within run-to-run noise, Optima′ is within 2% of Optima (ODE) and its agent form is within 0.35% of the ODE, and HIVsim runs on STIsim v4.
3. **The swap machinery works.** A runner, slot factories, endpoint-identity tests and an outputs analyzer are in place. The first path (HIVsim → EMOD′) ran end to end with exact endpoints, and the first step of the core path (the cascade, EMOD′ ↔ Optima′, both ways) is done.
4. **First result: the cascade carries most of the split by treatment status, not the gap in totals.** Swapping the cascade moves 0.5–0.75 of the gap in "diagnosed, never treated" and "on ART", and all of the gap in "interrupted". New infections move the other way, so the total must come from natural history, transmission, network or demography.
5. **The references are the bottleneck.** The EMOD binary with Daniel's two fixes reproduces the published shares, but its total is ~35–40% high, and we don't know which files were calibrated. Our Optima project is not the one MIHPSA published, and we need Rowan's.

## 1. How each model is decomposed

![Decomposition of EMOD′, HIVsim and Optima′ into slots](figures/decomposition.png)

Every model is cut into the same slots: D demography, R risk and sex work, N sexual contact, T transmission, H natural history and ART biology, S seeding, C cascade and prevention, M measurement, plus X solver for Optima. The first swaps use four coarse slots, **C, V (= T + H + S), RN, D**, because some slots are fused in code today. Details: `swap_design.md` §1.

## 2. The swap plan, and progress against it

![Swap plan and status](figures/swap_plan.png)

| Path | C | V | RN | D | Status |
|---|---|---|---|---|---|
| HIVsim → EMOD′ (first milestone) | ✅ | ✅ | ✅ | ✅ | All 5 hybrids run, 1,000 runs, endpoints exact (analysis 16). Uncalibrated, not yet compensated |
| **EMOD′ ↔ Optima′ (core)** | ✅ both ways (24) | ⛔ needs design | ⛔ needs design | ⛔ needs design | Prerequisites done: Optima′ as agents (20), population label and host (22), Optima's cascade as a table (23) |
| EMOD′ → HIVsim | ○ | ○ | ○ | ○ | Needs HIVsim's cascade as a table (the pathway types now exist, from 23) |
| HIVsim ↔ Optima′ | ○ | ○ | ○ | ○ | Becomes cheap once the core path exists |

✅ done · ⛔ blocked on a decision · ○ not started.

| Milestone in the plan | Planned | Status |
|---|---|---|
| CRN decision | sprint 3 | ✅ Pairing helps stocks, not shares (14); engine leak fixed (starsim#1574) |
| Component boundaries approved (P4) | sprint 4 | ◐ Drafted (`swap_design.md`, figure above); **for approval today** |
| Optima′ verified both ways | sprint 5 | ✅ ODE vs Optima 2% (P2), agents vs ODE 0.35% (20). Sign-off pending |
| Runner and endpoint test on a trivial swap | sprint 6 | ✅ (11, 15) |
| First path with endpoint identity | sprint 7 | ✅ Uncalibrated (16); compensated refit not built |
| EMOD′ ↔ Optima′ | sprint 8 / January | ◐ C step done (24) |

## 3. The models in numbers

See `model_stats.md` (run time, memory, size, fidelity).

## 4. Epidemiology

**What each model is calibrated to.** Prevalence by age and sex is the only indicator all three are fitted to. EMOD never used the UNAIDS PLHIV, new infections or AIDS deaths, and Optima only compared against them.

![Which indicators each model is calibrated to](figures/targets_matrix.png)

**Core indicators, all three models against the data.** New infections in 2020: EMOD′-fixed 38k, Optima 35k, HIVsim 26k, against UNAIDS 21k. ART coverage is close for all three. HIVsim's near-zero AIDS deaths after 2015 are a known gap: its reference set runs with the on-ART mortality options off (analysis 13; STIsim #613 has them, opt-in).

![Core indicators](figures/core_indicators.png)

**Figure 3A: published (EMOD, Optima) and reconstructed (all three).**

![Figure 3A reconstruction](figures/fig3a_reconstruction.png)

2020 shares (undiagnosed / diagnosed never treated / on ART / interrupted):

| | Shares | New infections, 2024 |
|---|---|---|
| Published EMOD | 39 / 4 / 41 / 17 | 30.6k |
| EMOD binary with Daniel's fixes | 37 / 4 / 45 / 13 | 41.3k (1.35×) |
| EMOD′ with Daniel's fixes | 36 / 3 / 48 / 13 | ~39k |
| Published Optima | 45 / 28 / 18 / 9 | 8.5k |
| Our Optima (2022 project) = Optima′ | 42 / 41 / 11 / 6 | 29.2k (3.5×; projection with inputs frozen at 2020) |
| HIVsim (provisional set) | 71 / 0 / 29 / 0 | 24.0k |

- **EMOD:** with Daniel's fixes, the published split is reproduced and the total isn't.
- **Optima:** neither is, because our project is a different calibration from MIHPSA's.
- **HIVsim's split is structural:** ART is filled to a number from whoever is diagnosed, so almost nobody stays diagnosed without treatment. There is no dropout until the coverage input falls (2021). That is a third answer to the Figure 3A question, and the swaps can explain it.

![New infections per year](figures/fig3a_totals.png)

## 5. Calibration stage gates

These are the ten steps of the calibration workflow. ✅ done · ◐ partly · ○ not started · — inherited from the original model.

| Gate | HIVsim (our calibration) | EMOD′-fixed | Optima′ | Compensated refit of hybrids (tier 2) |
|---|---|---|---|---|
| 1. Data | ✅ target set agreed (`targets.md`) | — EMOD's targets | — Optima's data | ✅ the HIVsim target set |
| 2. Visualise, model setup | ✅ 10,000 agents, seeds measured | ✅ 1% scale | ✅ | ◐ run cost known (stats table) |
| 3. Coverage check | ✅ exp 01 | — | — | ○ |
| 4. Targets and likelihood | ✅ stopping rule (3 Nov) | — | — | ◐ two knobs chosen (transmission scale, ART uptake) |
| 5. Workflow | ✅ experiments, READMEs, SUMMARYs | ✅ | ✅ | ✅ runner, store, manifest |
| 6. Edge cases | ◐ take-off on other demographies (21) | ✅ | ✅ | ◐ seeding (D11) |
| 7. Re-identification | ○ | — | — | ○ |
| 8. Parameter engineering | ◐ knobs provisional | — | — | ✅ two knobs |
| 9. Small real-data calibration | ✅ exp 03 | — 250 sets (Daniel) | — one fit | ○ |
| 10. Iterate | ◐ **exp 04: stopping rule not met** (2/600 runs); PLHIV before 2010 too low | ◐ **the two fixes invalidate the original fit:** score 108 → 7.5 uncalibrated, total ~40% high | ◐ **not MIHPSA's calibration** (19) | ○ |

Structural fixes found through calibration: declared STIsim ignored condoms (prevalence 19% → 8%); HIVsim's mortality on ART was 5–10× too low (13); EMOD's dropout trap and background mortality bug (Daniel, PR #4).

## 6. What we need to clarify

**References (critical path)**
1. **EMOD (Daniel):** which campaign and config files were calibrated for the paper? Is EMOD′-fixed, uncalibrated, the EMOD end, or do we refit its total, which is ~40% high?
2. **Optima (Rowan, via Robyn):** MIHPSA's Malawi project. Until it arrives we decompose against our 2022 project, stated as such.
3. **Year and outputs:** 2020 (inside the data, where Optima′ is verified) or 2024 (the caption's totals, which need Optima's post-2020 inputs)?

**Swap design**

4. **V, RN and D between EMOD′ and Optima′:** Optima has no network. Options: Optima's mean-field force of infection running on agents, or an agent network fitted to Optima's mixing. Or put transmission and network in one step.
5. **Slot boundaries (P4 sign-off):** PMTCT is in C (decided 10 Nov). Are EMOD's VMMC and PrEP in C, or prevention outside it? Analysis 23 has 11 more provisional decisions, plus ART arm switching on `EmodHost`.
6. **Seeding (D11):** condition on take-off, make seeding its own step, or seed more? The shares are robust to the choice; the counts aren't (21).
7. **Compensated refit:** the knobs and targets are agreed. Does the refit use the full stopping rule, or a cheaper subset, per hybrid?

**Calibration**

8. **HIVsim:** relax the stopping rule, or keep calibrating (early PLHIV)? Which runs go forward?

**Engineering and process**

9. **Starsim:** review of #1573 and #1574 (Cliff). How do we land the five local friction fixes?
10. **Sprint cycle:** keep weekly Tuesday check-ins with a readout each time?

## Sources

Analyses 01–24 (`analyses/NN_*/SUMMARY.md` on `starsim-rc3`); `docs/status_2026-11-09.md`; hivsim_malawi experiments 01–04.
