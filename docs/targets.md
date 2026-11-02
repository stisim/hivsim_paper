# P3 target set for HIVsim Malawi (draft for review)

**Status:** draft, 2 Nov 2026. Robyn approves this before calibrating (plan, P3). Once approved, the numbers go in `results/targets.csv`, and the same set is used to score EMOD and Optima.

**What it is for:** HIVsim is the third model of the comparison, so it has to fit the evidence that EMOD and Optima were fitted to. It doesn't need to fit better than they do. The decomposition is about outputs the data don't pin down (transmissions by cascade stage of the source, deaths on ART), so the targets cover what the data do pin down: population, prevalence by age and sex, the epidemic curve, and the cascade.

## Sources

| Source | What it gives | Kind | Notes |
|---|---|---|---|
| DHS 2004, 2010, 2015-16 | HIV prevalence by sex and 5-year age group, 15-49 | survey | In Optima's data (`hivprev`, years 2004, 2010, 2015), some age groups pooled |
| MPHIA 2015-16 | prevalence by sex and age, 15-64; share of PLHIV diagnosed, on ART, suppressed | survey | In Optima's data as 2016 (`hivprev`, and the cascade in `optpropdx`, `optproptx`, `optpropsupp`) |
| MPHIA 2020-21 | the same as MPHIA 2015-16 | survey | Not in Optima's data (the project's data ends in 2019-20). To add from the MPHIA 2020-21 report |
| FSW surveys (2006, 2014, 2018, 2019) | prevalence among female sex workers | survey | In Optima's data |
| UNAIDS estimates, 2025 round | PLHIV, new infections, AIDS deaths, 1990-2025; number on ART, 2010-2025 | model (Spectrum) | `models/hivsim/data/malawi_hiv_calib.csv` |
| Spectrum estimates in Optima's project | PLHIV, new infections, AIDS deaths with bounds, 1990-2019 | model (Spectrum) | An earlier round: PLHIV in 1990 is 406,000 here against 440,000 in the 2025 round |
| WPP | population by age and sex | estimate | Already an input (demography), so a check rather than a target |

## Proposed targets

Primary targets carry most of the weight. Secondary ones are model estimates, or they cover groups that are small or poorly measured.

| # | Target | Population | Years | Source | Role |
|---|---|---|---|---|---|
| 1 | HIV prevalence by sex and 5-year age group | 15-49 (DHS); 15-64 (MPHIA) | 2004, 2010, 2015-16 (DHS and MPHIA), 2020-21 | DHS, MPHIA | primary |
| 2 | Share of PLHIV who know their status | 15-64, by sex | 2015-16, 2020-21 | MPHIA | primary |
| 3 | Share of PLHIV on ART | 15-64, by sex | 2015-16, 2020-21 | MPHIA | check: ART is an input, so this mostly tests the size of the PLHIV denominator |
| 4 | Share of people on ART who are virally suppressed | 15-64 | 2015-16, 2020-21 | MPHIA | primary for suppression, if HIVsim's suppression isn't set from data instead (see question 5) |
| 5 | PLHIV | all ages | 1990-2021 | UNAIDS | secondary |
| 6 | New infections | all ages | 1990-2021 | UNAIDS | secondary, and the main check on the early epidemic (both originals peak too low before 2000) |
| 7 | AIDS deaths | all ages | 1990-2021 | UNAIDS | secondary |
| 8 | FSW prevalence | female sex workers | 2006, 2014, 2018, 2019 | FSW surveys | secondary, low weight (definitions and sampling differ between surveys) |
| 9 | Population | by age and sex | 1990-2021 | WPP | check only (an input) |

**Not targets:** the number on ART (an input: UNAIDS from 2010, Optima's data before); circumcision coverage (an input, once there is an `n_vmmc.csv`); anything after 2021, which is projection; the Optima and EMOD model outputs themselves, which are what the scorecard compares, not what HIVsim fits.

**Excluded from Optima's prevalence data until we know their source:** the 2018 points (one value shared by all ages 15-24, 25-34 and so on, and 1.09% at every age under 15, so probably a model estimate); 20% for every female group aged 50 and over in 2007; and the 2000 values for children.

## How the model is compared with them

- Prevalence by age and sex is compared in the survey's own age groups. Where Optima's data pools groups (DHS 2004 gives one value for 15-24), the model's prevalence is pooled the same way. HIVsim's `age_bins` will be set to 5-year groups from 15 to 65.
- Survey targets are proportions with known sample sizes, so their uncertainty comes from the survey's confidence intervals. UNAIDS series come with bounds; the 2025 round's bounds are in the UNAIDS downloads in `raw_data/`.
- The likelihood and weights are set with `calib:likelihood-design` once the set is approved. The weighting proposed here (surveys primary, UNAIDS secondary) is an input to that, not a final choice.

## The other two models

- **Optima** was fitted to the prevalence data above (in `hivprev`), with the Spectrum estimates and the 2016 cascade available for comparison. Its fit: PLHIV and population close; early new infections and AIDS deaths too low (see `optima_inventory.md`).
- **EMOD** was fitted to prevalence by age and sex, incidence, ART coverage and circumcision coverage (Bansi-Matharu 2025, Table 1), but we don't yet have its target files. **Needed from Daniel:** EMOD's Malawi calibration targets (sources, years, age groups) and the weights. If they include targets not listed here (incidence from MPHIA, ART by sex, ANC prevalence), they are added to this set.

## Calibration knobs (provisional, for the second half of P3)

Listed so the targets can be checked for coverage. Each target above should be moved by at least one knob, and every knob should move at least one target.

| Knob | Mechanism | Main targets |
|---|---|---|
| `beta_m2f` (with `rel_beta_f2m`, female to male relative to male to female, fixed at 0.5) | per-act transmission | prevalence level, new infections |
| `rel_init_prev` | size of the epidemic in 1985 | early PLHIV, new infections before 2000 |
| risk-group shares and concurrency (`prop_f0`, `prop_m0`, `f1_conc`, `m1_conc`) | sexual network | prevalence by age and sex, FSW prevalence |
| testing levels (general population and low CD4) | diagnosis, which ART needs | share diagnosed; the number on ART being reachable |

The ART input can only be reached if enough people are diagnosed, so testing is calibrated together with the epidemic, not separately.

## Questions for Robyn

1. **Which UNAIDS round?** The 2025 round (already in the repo) or the earlier Spectrum round in Optima's project? HIVsim would fit the 2025 round, but the models would then be scored against data from different rounds. Proposed: 2025 round for targets, and both rounds shown in the scorecard.
2. **Weights:** surveys primary, UNAIDS secondary, FSW low. Agree?
3. **MPHIA 2020-21:** I can't download it from this machine. Can you add the prevalence by age and sex and the cascade (to `raw_data/`)?
4. **The excluded Optima prevalence points** (2018, 2007 50+, 2000 children): leave out, or does Rowan know their source?
5. **Viral suppression:** set it from MPHIA as an input (`vls.csv`), or fit to it? Proposed: fit, since the share suppressed matters for transmissions on ART, one of the outputs the paper decomposes.
