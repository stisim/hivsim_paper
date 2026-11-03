# Sprint plan

*Weekly sprints, Tuesday to Monday. Read with `vision.md` and `prompts.md`, where the P numbers are defined.*

**Team:** Robyn (lead, Optima), Daniel (EMOD, prompts), Cliff (Starsim v4; Robyn reviews and contributes).

## Sprint 2 review (27 Oct–2 Nov)

Planned: both originals rerun with reference outputs, and draft prompts P1–P4. What happened:

- **Optima (Robyn), ahead of plan.**
  - P0 done: Malawi project rerun on Optima HIV 2.12.2, reference outputs in `optima_ref/`, component inventory in `optima_inventory.md`.
  - P2's first pass also done: Optima′ (`models/optima/`) is within 2% of Optima overall each year 2000–2020, except new infections in 2017 (2.2%) and 2020 (2.4%). By population, about half the 39 groups are 2–4% off somewhere, the worst 6–7% (85+ and 0–4 on ART). Not yet signed off.
  - Optima's ~6 components (demography, population forcing, natural history, transmission, cascade, treatment rules) agreed as a starting point for P4.
- **HIVsim (Robyn), ahead of plan.**
  - Re-ported onto Starsim 4.0.0rc1 (pinned; `envs/starsim.yml`), with ART split from HIV and an optional care module (STIsim PR #607, open). Its Optima-style options (care, loss and return, ART priority, failure and regain) default to HIVsim's behaviour, so the cascade swap is a change of settings.
  - P3, first half done: target set agreed (`targets.md`, `results/targets.csv`). Second half started: experiments 01–03 (`experiments/`). ART coverage and viral suppression can be matched, but only with late testing, and AIDS deaths on ART are too low.
- **EMOD (Daniel):** the Malawi calibrated model is in `EMOD_MWI/` (PR #1, merged 3 Nov): config, campaign, demographics and overlays, 250 resampled parameter sets, calibration targets and form, and the dtk-tools scripts that ran it on NYU's cluster. No outputs yet, so no rerun or reference outputs. Rerunning needs the EMOD binary it names (`EradicationMalaria-Ongoing_2023_01_05`), dtk-tools and a cluster.
- **Starsim (Cliff):** 4.0.0rc1 released (the pin). STIsim pinned to Starsim v3 on main (#608). Upstream fixes from the Optima′ build merged (Starsim #1492, #1498).
- **Not done:** P4 (component map and interfaces), the CRN test, and the prompt drafts (P1–P4 are still one-liners).

## Plan from sprint 3

| Sprint | Dates | Robyn | Daniel | Cliff | Deliverable |
|---|---|---|---|---|---|
| **1** | 20–27 Oct | Vision; narrative draft; plan | — | — | ✅ `vision.md`, `narrative_draft.md`, `plan.md`, `prompts.md` |
| **2** | 27 Oct–2 Nov | See the review above | | | ✅ Optima rerun; Optima′ first pass; HIVsim on v4; target set; EMOD files in PR #1 |
| **3** | 3–9 Nov | **P3:** finish calibrating HIVsim (exp 04: testing for advanced HIV from 2004; exp 05: mortality on ART; stopping rule). Draft P4 from the Optima and HIVsim side. Test CRN on one real swap (HIVsim's cascade set to Optima's) | **P0 EMOD:** rerun the 250 parameter sets, reference outputs (`ReportHIVByAgeAndGender`) summarised into `results/`, component inventory | Review and merge STIsim #607. Agree the pin for the rebuilds | **Checkpoint:** HIVsim Malawi calibrated to the stopping rule; CRN decision |
| **4** | 10–16 Nov | **P4** with Daniel: component map and interfaces (mechanism, inputs, parameters). Optima′ by population; sign-off | **P4** with Robyn; then start **P1:** EMOD′ rebuild loop | Review P4 for fit with v4 | **Checkpoint:** component boundaries and interfaces approved by all three; Robyn signs off Optima′ |
| **5** | 17–23 Nov | Verify Optima′ as an ABM against its own ODE (time-step study). Add the explicit re-engagement mechanism to HIVsim. Score EMOD and Optima on the target set | P1 continues | Unblock v4 issues from the rebuilds | Optima′ verified both ways; fit scorecard for all three models |
| **6** | 24–30 Nov | Refactor HIVsim and Optima′ to the interfaces. Swap runner and endpoint unit test (with Cliff) | First fidelity pass for EMOD′ (median inside the original's 90% interval) | Swap runner and endpoint unit test (with Robyn) | Swap runner with endpoint test passing on a trivial swap |
| **7** | 1–7 Dec | First substitution path, HIVsim ↔ Optima′, uncalibrated | Finish verifying EMOD′, or escalate | Support the swap runner | **Checkpoint:** Daniel signs off EMOD′; first path with endpoint identity confirmed |
| **8** | 8–15 Dec | Results outline and interim figures; plan for the January sprints | Paths EMOD′ ↔ Optima′ (the published gap) and HIVsim ↔ EMOD′ | Interim review | **Interim deliverable (15 Dec):** three verified models, first decompositions, interim figures, January plan |

Changes from the 27 Oct plan (agreed with Robyn, 3 Nov): P2 and P3 moved earlier (done or under way); P4 moved from sprint 3 to sprint 4, so P1 starts a week later; Optima′ sign-off moved from sprint 6 to sprint 4. EMOD′ is now the critical path.

## Still to schedule

These are left for January:
- the compensated (tier-2) pathway;
- the treatment-interruption scenarios;
- final swap design and orderings;
- the full manuscript, co-author round and submission (by 31 March 2027).

## Risks

| Risk | Fallback |
|---|---|
| Common random numbers don't reduce variance across swaps (`SexualNet` isn't CRN-safe when people are born or die); probably minor, since HIV outputs should vary little between runs | more seeds or agents; request a CRN-safe partnership network in Starsim |
| The EMOD′ rebuild can't reach tolerance | report the diagnosed cause; partial path (EMOD components in HIVsim, no endpoint test) |
| v4 API changes mid-rebuild | pinned to 4.0.0rc1; re-pin deliberately and rerun the endpoint test |
| EMOD′ (the critical path) starts a week late, after P4 | start EMOD′'s components that don't depend on the interfaces (demography, natural history) during sprint 3–4; partial path if needed |
| HIVsim calibration runs on (a calibration can always be improved) | the stopping rule in `experiments/03_testing_art/SUMMARY.md`; whatever meets it by 9 Nov goes forward |
| HIVsim's low mortality on ART is a model gap, not a parameter | fix it in STIsim before the swaps, since deaths on ART are a paper output |
| ~~MIHPSA Optima turns out to be legacy Optima HIV, not Atomica~~ Confirmed 28 Oct: legacy Optima HIV ([optimamodel/optima](https://github.com/optimamodel/optima)) | hand rebuild (declared, for its ODE form); Malawi files from Rowan |
