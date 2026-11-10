# Sprint plan

*Weekly sprints, Tuesday to Monday, each closed by a team check-in on the Tuesday with a readout in `docs/readouts/` (branch `starsim-rc3`). Read with `vision.md`, `prompts.md` (the P numbers), `core_question.md` and `swap_design.md` (the slots, paths and decisions D1–D11).*

*Updated 10 Nov 2026 (end of sprint 3).*

**Team:** Robyn (lead, Optima), Daniel (EMOD, prompts), Cliff (Starsim v4; Robyn reviews and contributes).

## Sprint 3 review (3–9 Nov)

Planned: HIVsim calibrated to the stopping rule; P4 drafted; a CRN test; EMOD rerun with reference outputs (Daniel); STIsim #607 merged and a pin agreed (Cliff). What happened:

- **The question sharpened (Robyn, 9 Nov):** the core question is now Figure 3A of Bansi-Matharu et al. 2025: why EMOD and Optima HIV disagree, in Malawi, on who new infections come from (undiagnosed, diagnosed never treated, on ART, interrupted) and on how many there are. HIVsim becomes the common architecture rather than one of the two ends, and **EMOD′ ↔ Optima′ becomes the core path** (`core_question.md`).
- **Engineering is well ahead of plan.** Items planned for sprints 5–8 and January are done (analyses 10–24 on `starsim-rc3`):
  - EMOD′ reproduces EMOD within run-to-run noise, piece by piece and whole (analyses 02–06, 09). The condensed cascade is the common cascade interface.
  - Swap runner, slot factories, endpoint tests and the outputs analyzer (`common/swaps/`, `common/outputs.py`; analyses 11, 15). This was the sprint 6 deliverable.
  - **First path HIVsim → EMOD′ run end to end with exact endpoints** (analysis 16, 1,000 runs). This was the sprint 7 checkpoint, though uncalibrated rather than compensated.
  - Optima′ as agents agrees with its ODE (≤0.35% at 400k agents; analysis 20; sprint 5). Population label, `OptimaHost`, Optima's cascade as a table (analyses 22–23; January items).
  - **First step of the core path: the cascade swapped both ways between EMOD′ and Optima′** (analysis 24). The cascade carries most of the share split, not the gap in totals. Optima's PMTCT, skipped on agent hosts as a shortcut, is now on by default (Robyn, 10 Nov).
- **The references moved.** Daniel's two fixes (PR #4: the dropout trap and background mortality) make the EMOD binary reproduce the published Figure 3A shares within 5 points (analysis 18; EMOD′-fixed agrees, analysis 17), but its total is ~35–40% high. Our Optima (2022 project) is not the calibration MIHPSA published (analysis 19).
- **HIVsim calibration: not done.** Exp 04 (testing for advanced HIV, EMOD's mortality on ART): 2 of 600 runs meet the stopping rule, PLHIV before 2010 too low. A provisional reference set is in `models/hivsim/reference/` (endpoint identity doesn't depend on fit). Along the way two STIsim bugs were found and fixed: declared STIsim ignored condoms (prevalence 19% → 8%) and HIVsim's mortality on ART was 5–10× too low (analysis 13).
- **CRN decided:** pairing helps stocks a lot, flows a little, cascade shares not at all (analysis 14); a Starsim engine leak fixed (starsim#1574). Use pairing and more seeds.
- **Starsim / STIsim:** pin is Starsim 4.0.0rc3 + #1573 (under Cliff's review) and STIsim #613 (on #612); #607 and #609 closed. Five more Starsim frictions fixed on a local branch (not pushed).
- **Not done:** P4 sign-off by all three (the slot boundaries in `swap_design.md`, figure in this sprint's readout); Optima′ sign-off; Daniel's P0 (we reran the binary ourselves, analyses 01, 09, 18); the compensated (tier 2) refit in the runner.

## Sprint 2 review (27 Oct–2 Nov)

Planned: both originals rerun with reference outputs, and draft prompts P1–P4. What happened:

- **Optima (Robyn), ahead of plan.**
  - P0 done: Malawi project rerun on Optima HIV 2.12.2, reference outputs in `models/optima/reference/`, component inventory in `optima_inventory.md`.
  - P2's first pass also done: Optima′ (`models/optima/`) is within 2% of Optima overall each year 2000–2020, except new infections in 2017 (2.2%) and 2020 (2.4%). By population, about half the 39 groups are 2–4% off somewhere, the worst 6–7% (85+ and 0–4 on ART). Not yet signed off.
  - Optima's ~6 components (demography, population forcing, natural history, transmission, cascade, treatment rules) agreed as a starting point for P4.
- **HIVsim (Robyn), ahead of plan.**
  - Re-ported onto Starsim 4.0.0rc1 (pinned; `envs/starsim.yml`), with ART split from HIV and an optional care module (STIsim PR #607, open). Its Optima-style options (care, loss and return, ART priority, failure and regain) default to HIVsim's behaviour, so the cascade swap is a change of settings.
  - P3, first half done: target set agreed (`targets.md`, `common/targets.csv`). Second half started: experiments 01–03 (`experiments/` in hivsim_malawi). ART coverage and viral suppression can be matched, but only with late testing, and AIDS deaths on ART are too low.
- **EMOD (Daniel):** the Malawi calibrated model is in `models/emod/reference/source/` (PR #1, merged 3 Nov): config, campaign, demographics and overlays, 250 resampled parameter sets, calibration targets and form, and the dtk-tools scripts that ran it on NYU's cluster. No outputs yet, so no rerun or reference outputs. Rerunning needs the EMOD binary it names (`EradicationMalaria-Ongoing_2023_01_05`), dtk-tools and a cluster.
- **Starsim (Cliff):** 4.0.0rc1 released (the pin). STIsim pinned to Starsim v3 on main (#608). Upstream fixes from the Optima′ build merged (Starsim #1492, #1498).
- **Not done:** P4 (component map and interfaces), the CRN test, and the prompt drafts (P1–P4 are still one-liners).

## Plan from sprint 4

Sprints 1–3 are in the reviews above. EMOD′ is no longer the critical path; the critical path is now **the references** (which EMOD and which Optima we decompose) and **the design of the V, RN and D steps between EMOD′ and Optima′** (Optima's transmission has no network).

| Sprint | Dates | Robyn | Daniel | Cliff | Deliverable |
|---|---|---|---|---|---|
| ✅ **1–3** | 20 Oct–9 Nov | See the reviews above | | | ✅ Core question; EMOD′ verified; Optima′ ODE and agents; runner; HIVsim → EMOD′ path; cascade swap EMOD′ ↔ Optima′ |
| **4** | 10–16 Nov | **P4 checkpoint:** slot boundaries and interfaces (D7, the decomposition figure). Decide the V/RN/D design for EMOD′ ↔ Optima′, D11 (seeding) and analysis 23's provisional decisions. Sign off Optima′. HIVsim exp 05: early PLHIV, or relax the stopping rule | Which EMOD files were calibrated; whether EMOD′-fixed needs a refit for the total (~40% high). Review analyses 17–18 and 24 (EMOD side) | Review and merge Starsim #1573 and #1574; STIsim #612/#613; take the friction batch | **Checkpoint:** boundaries approved by all three; references for both ends agreed (or a dated plan to get them) |
| **5** | 17–23 Nov | EMOD′ ↔ Optima′: the V step (natural history and transmission), both ways, on the agreed design. Compensated refit (tier 2: transmission scale, ART uptake) in the runner | Sign off EMOD′-fixed as the EMOD end (or refit it) | Unblock v4 issues; re-pin once #1573 lands, rerun the endpoint tests | V step both ways; runner refits each hybrid |
| **6** | 24–30 Nov | EMOD′ ↔ Optima′: RN and D steps; full path with endpoint identity; Optima′ on MIHPSA's project if Rowan's files arrive | Review the EMOD side of the path | Support the runner | **Checkpoint:** core path complete, uncalibrated and compensated |
| **7** | 1–7 Dec | Shapley lattice on the core path; HIVsim ↔ Optima′ and EMOD′ → HIVsim if cheap; treatment-interruption scenario set up (re-engagement rate in HIVsim) | Sign off the EMOD side | Interim review of the code | Decomposition of Figure 3A (shares and total) with intervals |
| **8** | 8–15 Dec | Results outline and interim figures; January plan | Interim review | Interim review | **Interim deliverable (15 Dec):** verified ends, the core path decomposed, interim figures, January plan |

## Still to schedule

These are left for January:
- the treatment-interruption scenarios and the intervention ranking (2026–2035);
- finer slots (V split into T, H, S) after the refactor, and the final orderings;
- the full manuscript, co-author round and submission (by 31 March 2027).

## Risks

| Risk | Fallback |
|---|---|
| **MIHPSA's Optima Malawi project doesn't arrive from Rowan** | decompose against our 2022 project (Optima′ = our Optima), stated as such; the published panel is a check, not the end |
| **The calibrated EMOD files are unknown, and EMOD′-fixed's total is ~40% high** | use the fixed files uncalibrated with the compensated refit (tier 2) as the headline; report the total gap with that caveat |
| Optima has no network, so V/RN/D steps on EMOD′ ↔ Optima′ need a design choice (mean-field force of infection on agents, or an agent network fitted to Optima's mixing) | do the cascade and natural history steps first; put transmission and network together in one step |
| Tests and runs slow the loop (the test suite takes ~35 min serially) | short horizons in tests (to 1995, `swap_design.md` §3), parallel runs; full horizon once per milestone |
| ~~Common random numbers don't reduce variance across swaps~~ Checked (analysis 14): they help stocks, not cascade shares | more seeds (cheap); engine leak fixed (starsim#1574) |
| ~~The EMOD′ rebuild can't reach tolerance~~ Done: EMOD′ reproduces EMOD within noise (analyses 02–06, 09) | |
| v4 API changes mid-rebuild | pinned to 4.0.0rc1; re-pin deliberately and rerun the endpoint test |
| HIVsim calibration runs on (a calibration can always be improved) | the stopping rule in `experiments/03_testing_art/SUMMARY.md` in hivsim_malawi; whatever meets it by 9 Nov goes forward |
| HIVsim's low mortality on ART is a model gap, not a parameter | fix it in STIsim before the swaps, since deaths on ART are a paper output |
| ~~MIHPSA Optima turns out to be legacy Optima HIV, not Atomica~~ Confirmed 28 Oct: legacy Optima HIV ([optimamodel/optima](https://github.com/optimamodel/optima)) | hand rebuild (declared, for its ODE form); Malawi files from Rowan |
