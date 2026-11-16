# Notes for Daniel: EMOD-HIV Malawi rerun and EMOD′ (4 Nov 2026)

Daniel, thanks for the model. Here's where we've got to, the problems we think we've found in it, and what we'd like you to check. Details are in the files linked.

## Where things stand

- **Rerun:** your 250 parameter sets run here with the published `emod-hiv` 2.37.0 binary (`models/emod/README.md`), to 2025. Sets you scored higher also fit our common targets better, weakly (rank correlation 0.21), so the 2.37 build seems to behave like yours (`analyses/01_reference_fits/`).
- **EMOD′:** your model re-expressed in Starsim v4 (`models/emod/v4/`). It reproduces EMOD within run-to-run noise: demography and network (`analyses/02`), transmission and natural history (`03`), and the full model with the care cascade, state by state in InterventionStatus (`04`). The campaign is read three ways, all of which match: interpreted literally (`04`), compiled into Starsim deliveries (`05`), and condensed into about ten interventions with a parameter table (`06`). We copy the problems below rather than fix them, so EMOD′ is your published model.

## Problems we think we've found

Confirmed in the files and in EMOD's output:
1. **Background mortality is about 1000 times too low.** `MortalityDistributionMale/Female` in `Demographics_Malawi.json` hold annual deaths per person (infant mortality 0.18 in 1960), but `ResultScaleFactor` is 1/(365×1000), the factor for deaths per 1000. Non-HIV deaths come out at about 0.03 per 1000 a year, and the population is about twice the census by 2008. Your calibration rescaled counts to the 2018 population, which hides the size but not the growth or the age structure.
2. **People who drop out of ART can't restart, from 2009.** Dropout sets `ARTStatus` to OffART (event 76) but leaves `InterventionStatus` at OnEffectiveART or OnNonSuppressiveART (set at ART start, events 66 and 68). The re-enrolment branch (events 77 and 79) sends 85–95% of dropouts to `OnART1`, but the four ART-start choices (events 62–65) are disqualified for exactly those statuses, as are the testing and linkage pathways. By 2020, 42% of people with HIV are off ART, and ART coverage levels off at about 40%. This is probably the main reason coverage is far below the UNAIDS series.
3. **Combination PMTCT (Option B) is effectively blocked for women found at ANC.** A positive ANC test also sends them into ART staging, which sets `ARTStaging` before the combination decision is updated, and that status disqualifies it. EMOD confirms this: if the block is removed in EMOD′, infant infections in 2009 fall from about 240 to about 90 a year. After about 2008 almost nobody gets Option B.
4. **The ANC testing table is ignored.** The ANC uptake intervention is a sigmoid class that reads only its `Ramp_*` keys, so a default curve applies; ANC testing in 2001 is 0.73 instead of the intended 0.31.
5. **PrEP doesn't exclude people already on PrEP.** The NChooser excludes an intervention named `PrEP`, but the intervention is named `Oral_PrEP`.
6. **At 1% population, NChooser targets are truncated per node**, e.g. the 2011 VMMC target in the Northern region (54 men) becomes zero.

Found by reading the files, not yet confirmed in output:
7. **Calibrated parameters with no effect:** `Pr Ex Infmrl Fem LOW` maps to the male parameter, so it overwrites `Pr Ex Infmrl Male LOW` and the female value is never set; `Trns Condom Max` and `Infmrl Condoms Late` are overwritten by the per-region condom values.
8. **The male non-suppressive ART mortality table** has a 182–365-day row about 7 times the female row and its neighbours. It looks like a typo.
9. **The Risk overlay lists its nodes in a different order**, so each region's LOW-risk fraction goes to another region. Provincial prevalence had no weight in your fit, so this mislabels regions rather than changing the fit.
10. **Parameter sets and campaign dates:** the sets are from April 2023 and the config template names a 2022 campaign, but the campaign file is dated December 2025. Were the sets fitted with this campaign?

## What we'd like you to check

- **The list above:** are these real problems in the model as published, or have we misread something? In particular 2, 3 and 10.
- **The condensed parameter table** (`common/cascade/tables/malawi_emod_set000_params.csv`, with VMMC and PrEP in `..._flows.csv`): one row per cascade parameter, each naming the campaign event it came from. Does it describe your cascade as you intended it?
- **Duplicates:** your 250 resampled sets contain 129 distinct parameter vectors (one repeated 22 times), told apart by Run_Number. Is that as expected from your resampling?

Full component inventory: `docs/emod_inventory.md`.

## Review of your PR #4 (9 Nov 2026; not yet sent)

Thanks for the responses (`models/emod/reference/source/emod_responses.md`) and fixes. What we checked:

- **Mortality** (`Demographics_Malawi.json`, `ResultScaleFactor` 1/(365·1000) → 1/365, both sexes): agreed, this is the right factor for deaths per person.
- **Dropout trap** (new campaign event: on `OnART3`, the `ARTDropout` event, set `InterventionStatus` to `None`): this does release people who drop out. The ART-start events disqualify `OnEffectiveART`/`OnNonSuppressiveART`, `None` is a declared value of `InterventionStatus`, and the re-enrolment choice (`HCTUptakePostDebut8` or `LostForever9`) still follows on the same trigger. (`emod_responses.md` says "back to OffART"; `OffART` is a value of `ARTStatus`, not `InterventionStatus`, so `None` is right.)
- Path changes in `script.sbatch` and `simtools.ini` match the repo layout. The PR merges into `main` without conflicts.

Questions back:

1. **Which files produced the published (MIHPSA) results and the 250 calibrated sets?** You think your original calibration used the correct mortality factor, and that the dropout trap came in with the move to `ARTMortalityTable`. If both are later than the calibration, the files in PR #1 (which we treated as the published model, and which EMOD′ reproduces, bugs included) are not the model that was calibrated. Do you know which campaign/demographics files (or dates/commits) the 250 sets were fitted with, and which produced the MIHPSA submission?
2. **Does the calibration still hold with the fixes?** Rescaling to the 2018 census keeps prevalence and ART coverage, but the age structure and growth change, and the dropout fix changes ART coverage and the share of transmissions from people who interrupted ART (50–55% with the trap vs 8–20% across MIHPSA models). Would you expect to refit, or treat the fixed model as uncalibrated?
3. Minor: the second bullet in `emod_responses.md` stops mid-sentence ("it is unlikely this would have a").

On our side: EMOD′ already has the two bugs as switches; we're making a "fixed" EMOD′ variant that matches your two changes exactly (D10 in `docs/swap_design.md`), and will compare it with the copied-bug version.
