# P3 target set for HIVsim Malawi

**Status:** agreed with Robyn, 2 Nov 2026. The values are in `results/targets.csv`, made by `models/hivsim/make_targets.py` from `models/hivsim/data/sources/`. The same set is used to score EMOD and Optima.

**What it is for:** HIVsim is the third model of the comparison, so it has to fit the evidence that EMOD and Optima were fitted to. It doesn't need to fit better than they do. The decomposition is about outputs the data don't pin down (transmissions by cascade stage of the source, deaths on ART), so the targets cover what the data do pin down: prevalence by age and sex, the size of the epidemic over time, and prevalence among sex workers.

## Targets

| Target | Population | Years | Source | Role |
|---|---|---|---|---|
| HIV prevalence by sex and 5-year age group | 15-49 (DHS); 15-64 (MPHIA 2015-16); 15-64 and 65+ (MPHIA 2020-21) | DHS 2004, 2010, 2015-16; MPHIA 2015-16, 2020-21 | EMOD's targets (DHS, MPHIA 2015-16); MPHIA 2020-21 report, Table 6.3 | primary |
| PLHIV | all ages | 1990-2021 | UNAIDS, 2025 round | secondary |
| New infections | all ages | 1990-2021 | UNAIDS, 2025 round | secondary, and the main check on the early epidemic (EMOD and Optima both peak too low before 2000) |
| AIDS deaths | all ages | 1990-2021 | UNAIDS, 2025 round | secondary |
| HIV incidence | 15-49, by sex | MPHIA 2015-16, 2020-21 | EMOD's targets; MPHIA 2020-21 report, Table 5.1 | secondary (wide intervals). One of EMOD's targets |
| Prevalence among female sex workers | FSW | 2006, 2014, 2018, 2019 | Optima's data | secondary |

**Checks, not fitted:**
- **Number on ART by sex,** 2004-2015 (programme data via Spectrum, one of EMOD's targets). ART is an input only as a total, so the split by sex is a check.
- **MPHIA 2020-21 cascade:** share diagnosed, share of those on ART, and share suppressed. HIVsim matches the number on ART to its input, so these shares follow from PLHIV and can't be fitted separately (Robyn). The way to get them right is to get PLHIV right.
- **Population** by sex and 5-year age group, 1990-2020 (WPP). Demography is an input.

**Inputs, not targets:** the number on ART (UNAIDS from 2010, Optima's data before); viral suppression among people on ART (`data/vls.csv`, from MPHIA, see below); circumcision coverage, once there is an `n_vmmc.csv`.

**Left out:** anything after 2021 (projection); the earlier Spectrum round in Optima's project, which is out of date; and Optima's prevalence points of unknown source (2018, one value shared by all ages 15-24, 25-34 and so on; 20% for every female group aged 50 and over in 2007; children in 2000). Rowan hasn't confirmed their source.

## How the model is compared with them

- **Prevalence: the shape matters more than the exact level** (Robyn). The age and sex pattern and the trend between surveys should fit; a survey's level can be missed by a little. The likelihood should reflect that, for example by widening the survey intervals or giving the age pattern more weight than the level.
- Prevalence is compared in the survey's own groups. HIVsim records results in 5-year age groups from 15 to 65 (`age_bins` in `malawi.py`).
- **Uncertainty:** the DHS and MPHIA 2015-16 values carry EMOD's effective sample sizes and 95% intervals. For MPHIA 2020-21 the intervals come from the number tested (Jeffreys), so they leave out the survey design effect and are narrower than the survey's own. UNAIDS bounds aren't in the file yet (they're in the UNAIDS downloads in `raw_data/`).
- **Dates:** each survey is dated at the middle of its fieldwork: DHS 2004 2004.9, DHS 2010 2010.7, DHS 2015-16 2015.9, MPHIA 2015-16 2016.2, MPHIA 2020-21 2020.5 (fieldwork January-March 2020 and March-April 2021). EMOD's files label these 2004, 2010, 2015.5 and 2015.
- **Weights:** surveys primary, UNAIDS secondary, sex workers secondary. The likelihood and the numbers are set with `calib:likelihood-design`.

## Viral suppression (input)

`data/vls.csv` sets the share of people on ART who are suppressed: 91.3% in 2015-16 and 96.9% in 2020-21 (MPHIA 2020-21 report, Table 9.1.B), interpolated between. Before 2016 it holds at 91.3% back to the start of ART, since there is no national measurement earlier.

The 2015-16 value is from Optima's data (`optpropsupp`, 2016). Its three values (94.6%, 91.3%, 88.0%) are in the order high, best, low rather than Optima's best, low, high; 91.3% is the middle one and agrees with MPHIA 2015-16. Optima's share diagnosed for 2016 has the same order (81.0%, 76.8%, 72.6%).

## The other two models

- **Optima** was fitted to prevalence by population (its `hivprev` data), with the Spectrum estimates and the 2016 cascade for comparison. Its fit: PLHIV and population close; early new infections and AIDS deaths too low (see `optima_inventory.md`).
- **EMOD**'s targets are in its calibration files (PR #1, `EMOD_MWI/calib_targets_data/`): prevalence by sex and age (DHS 2004, 2010, 2015-16; MPHIA 2015-16), incidence (MPHIA 2015-16 and 2020-21), numbers on ART by sex (2004-2015), ART coverage by sex and age (MPHIA 2015-16), and population (WPP). Prevalence and incidence are targets here, numbers on ART by sex and population are checks, and ART coverage by sex and age isn't used, for the same reason as the MPHIA cascade. EMOD also has prevalence by region (`ARTCoverage.csv`, despite its name), which this national set leaves out. Its weights are in its calibration form (`EMOD_MWI/Data/calibration_ingest_form_Malawi_3_node_DTC.xlsm`), not yet read.

## Calibration knobs (provisional, for the second half of P3)

Each target should be moved by at least one knob, and every knob should move at least one target.

| Knob | Mechanism | Main targets |
|---|---|---|
| `beta_m2f` (with `rel_beta_f2m`, female to male relative to male to female, fixed at 0.5) | per-act transmission | prevalence level, new infections |
| `rel_init_prev` | size of the epidemic in 1985 | early PLHIV, new infections before 2000 |
| risk-group shares and concurrency (`prop_f0`, `prop_m0`, `f1_conc`, `m1_conc`) | sexual network | prevalence by age and sex, FSW prevalence |
| testing levels (general population and low CD4) | diagnosis, which ART needs | the number on ART being reachable |

Testing has no target of its own now that the cascade shares are checks. It has to diagnose enough people for the ART input to be reached, so it is fitted with the epidemic, not separately.
