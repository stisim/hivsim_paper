# EMOD-HIV Malawi: component inventory

Source: Daniel Citron's calibrated Malawi model in `models/emod/reference/source/` (config template, campaign `campaign_Malawi_20251215.json`, demographics and three overlays, 250 parameter sets in `resampled_parameter_sets_20230424.csv`, calibration form `Data/calibration_ingest_form_Malawi_3_node_DTC.xlsm`). We run it on the `emod-hiv` 2.37.0 binary via `models/emod/reference/run_emod.py`; a built input set is in `raw_results/emod/sub10/set000/`. Mechanisms below are from the EMOD source (github.com/EMOD-Hub/EMOD, `main` at 01f8b6f, Nov 2026). The binary is a build of `main` from Aug 2026, so line references may be slightly off for it.

Values are given for parameter set 0 unless they were calibrated. For calibrated parameters the table shows the median and [min–max] across the 250 sets.

## Structure

- **Agent-based**, 1960.5 to 2100 (we stop at 2025). The time step is 30.42 days (monthly).
- **3 nodes** (1 Northern, 2 Central, 3 Southern), with no migration. Nothing links the nodes apart from shared parameters, so they are effectively three separate epidemics with some node-specific inputs.
- **Population scale** `x_Base_Population` = 0.01, so there are 36,600 agents in 1960 (3.66M people). Every agent has weight 1 (`FIXED_SAMPLING`, rate 1). The calibration rescaled modelled population and number on ART to 8.26M aged 15–49 in 2018 (the form's Site sheet), so the absolute size of the population was never fitted.
- **Random numbers**: `USE_AES_COUNTER`, one stream per node. The seed is `Run_Number`, a random seed stored with each parameter set. The calibration ran 3 replicates per set. We run one.
- **Individual properties (IPs)**: everyone has `Risk` (LOW/MEDIUM/HIGH) and `InterventionStatus` (17 values, the cascade state). `ARTStatus` (ARTNaive/OffART/OnEffectiveART/OnNonSuppressiveART) is used only for reporting. `Accessibility` and `OnPrEP` exist but have no effect (see "Switched off").
- **Campaign**: 89 events. Nearly all of them are `NodeLevelHealthTriggeredIV` listeners that hand out one intervention when a named event fires. The interventions are diagnostics or random choices that fire further events. The care cascade is a state machine built out of these events.

## Calibration

- **Method**: OptimTool (dtk-tools), with 83 samples per iteration, 3 replicates and up to 25 iterations. The 250 sets come from iterations 4–16, and their log-likelihoods range from 93.6 to 98.3.
- **Parameters**: 55 dynamic and 1 static (circumcision efficacy, fixed at 0.6). They are listed under each component below.
- **Targets and weights**:

| Target | Likelihood | Weight |
|---|---|---|
| Prevalence, national, by age and sex (DHS 2004 and 2010, MPHIA 2015–16) | Beta | 0.8 |
| ART coverage | Beta | 0.15 |
| Population (rescaled) | Gaussian | 0.025 |
| Number on ART (rescaled) | Gaussian | 0.025 |
| Incidence; provincial prevalence | Beta | 0 |

- **Quirks in how parameters reach the model** (reproduced in `run_emod.py`):
  - `Pr Ex Infmrl Fem LOW` maps to the male parameter, so the sampled `Pr Ex Infmrl Male LOW` is overwritten and never used. The female value stays at the template value of 0.
  - **`Trns Condom Max` and `Infmrl Condoms Late` have no effect.** They write to the `Society__KP_Defaults` marker in every node. The per-node `<Region>: Trns/Infmrl Condoms Late` parameters come later in the sheet and write the same field, so they overwrite them (checked in set 0: each node's maximum equals its regional value).
  - **The per-region `LOW Risk` fractions go to the wrong nodes.** The Risk overlay labels node 1 "Southern", node 2 "Northern" and node 3 "Central". The demographics, the PFA overlay, the campaign and the calibration form all have node 1 = Northern, 2 = Central, 3 = Southern. So `Southern: LOW Risk` sets node 1 (Northern), `Northern:` sets node 2 (Central), and `Central:` sets node 3 (Southern). Provincial prevalence had zero weight, so this changes the meaning of the labels but not the fit.
  - When `PreART Link Min` > `Max` the two are swapped (8 sets).
- **Which campaign was used for calibration isn't known.** The parameter sets are dated 2023-04-24. The config template still refers to `campaign_Malawi_20220422_recalib.json`, and the only campaign in the folder is dated 2025-12-15. The sets may have been fitted with a different campaign from the one we run. Not checked.
- **Parameters at their bounds** across the 250 sets, which suggests the fit pushed against the ranges:
  - The transitory formation rate is at its lower bound in more than half the sets.
  - `ART Link Max` reaches 1, and `Male To Female Young` reaches 5.
  - `Southern: LOW Risk` reaches 0.2, and `Trns Condom Mid` reaches 2010.
  - `Infrml Dur Het` reaches 0.9999.
  - The regional condom maxima reach 0.8 (transitory) and 0.6 (informal).
  - The maximum number of transitory partners reaches 2 (LOW) and 3.5 (MEDIUM).

## Components

### 1. Demography

- **Initial population**:
  - 475,800 / 1,573,800 / 1,610,400 people in the three nodes (1960).
  - Ages are drawn from a 0–90 age distribution (`DISTRIBUTION_COMPLEX`).
  - Nobody is infected at the start.
- **Fertility** (`INDIVIDUAL_PREGNANCIES_BY_AGE_AND_YEAR`):
  - Each woman becomes pregnant at an age- and year-specific rate (WPP 2022 births per 1000 women; 7 five-year age groups, 1950–2100). The birth happens 40 weeks later.
  - The scale factor is correct (1/365,000).
  - The demographics file notes that "fertility still ends before HIV had a demographic impact".
  - Whether HIV affects fertility was not checked; EMOD has no such parameter here.
- **Background mortality** (`NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER`):
  - Table by sex, 21 age groups (0–100) and year (1952.5–2062.5).
  - The table holds deaths per person per year, but `ResultScaleFactor` = 1/365,000. **Every rate is therefore 1000 times too low** (known; `models/emod/README.md`).
- **Migration**: none (`NO_MIGRATION`).
- **Outputs driven**: population size and age structure, and through them births, ART counts and everything scaled. Because of the mortality bug the population grows about 4% a year and is about twice the census by 2008.

### 2. Risk groups and sex work

- **`Risk` is drawn at birth** from the node's initial distribution: LOW with probability v, MEDIUM with 1 − v, HIGH never. There are no transitions apart from those made by the campaign.

| Parameter | Value |
|---|---|
| LOW fraction, node 1 (from `Southern: LOW Risk`) | 0.243 [0.200–0.382] |
| LOW fraction, node 2 (from `Northern: LOW Risk`) | 0.388 [0.306–0.465] |
| LOW fraction, node 3 (from `Central: LOW Risk`) | 0.373 [0.301–0.448] |

- **Female sex workers (FSW) and clients are people whose `Risk` is HIGH** (campaign events 0–12):
  - At sexual debut (`STIDebut`), a woman "will ever be FSW" with probability 0.062 / 0.082 / 0.098 (nodes 1/2/3). A man "will ever be a client" with probability 0.088 / 0.143 / 0.239.
  - The delay from debut to uptake is U(0, 5 y) for women and U(0, 10 y) for men. At uptake `Risk` becomes HIGH.
  - Time in sex work is Weibull (shape 3.31, scale 2215 d; mean about 5.4 y) for women and U(2, 30 y) for men. At the end `Risk` becomes **MEDIUM**, even for people who were LOW before.
- **Changing `Risk` redraws the person's concurrency limits and flags** and clears their cached assortativity group (`STIInterventionsContainer.cpp:140`). Their existing relationships continue.
- **Outputs driven**: FSW and client prevalence (2019 FSW prevalence is about half the survey value). Commercial sex is also the route by which the epidemic is seeded (section 12).

### 3. Sexual network (pair-formation algorithm, PFA)

There are four relationship types: TRANSITORY, INFORMAL, MARITAL and COMMERCIAL. Each type has its own PFA, rates and duration (`PFA_Overlay_Malawi.json`). Relationships are always between a man and a woman, and both partners must be in the same node.

**How a relationship forms** (`NodeSTI.cpp:145`):

1. **Eligibility** (`IndividualSTI.cpp:1070`). A person is available for type *t* if all of the following hold:
   - they are past their sexual debut age;
   - for every type in which they already have a relationship, they hold the extra-relational flag for that type;
   - their active plus queued relationships of type *t* number fewer than their maximum for *t*.
2. **Entry rates** (`FlowControllerImpl.cpp:17,111`):
   - The desired number of entries is `Formation_Rate_Constant` (per day) × the number of eligible people. It is split half to men and half to women, and across age bins by the marginals of the joint age table.
   - The per-person rate in an age-sex bin is desired / (eligible without flag + 10 × eligible with flag). People who hold the extra-relational flag for that type enter at 10 times the rate (`Extra_Relational_Rate_Ratio`).
   - **These per-person rates are recomputed only during the PFA burn-in** (`PFA_Burnin_Duration_In_Days` = 5475, i.e. until 1975.5; `NodeSTI.cpp:160`). After that they are frozen at their 1975 values, so formation per person no longer responds to the size or age structure of the population.
3. **Entering the queue** (`IndividualSTI.cpp:761`). Each step, each available person draws exponential waiting times at the sum of their rates over types. They join the queue of one type per event, until the end of the step or until they are no longer available.
4. **Matching** (`BehaviorPfa.cpp:65`, every `Update_Period`; with 7 or 30 days and a monthly step, this is every step). For each man in the queue:
   - Weight the female age bins that have someone queued by the row of the joint age table (20 × 20 bins of 2.5 y from 17.5 for the three main types; one bin 15–50 for commercial). A bin whose first woman has `EnterRelationshipNow` set gets 10 times the weight.
   - If the summed weight is above `PFA_Cum_Prob_Selection_Threshold` = 0.2, pick a bin.
   - Pick a woman in that bin by the `Risk` assortativity weights (`Assortivity.cpp:437`): each candidate is weighted by the matrix entry for (his Risk, her Risk) and one is drawn.
   - Men left unmatched stay in the queue.
5. **Duration** (`Relationship.cpp:157`) is Weibull(scale, heterogeneity = 1/shape), drawn at formation. A relationship ends when its timer runs out or a partner dies.

**Concurrency** (`ConcurrencyParameters.cpp:324,365`):

- **Extra-relational flags** are drawn at debut and again whenever `Risk` changes. They are correlated: types are drawn in the order COMMERCIAL → TRANSITORY → INFORMAL → MARITAL, and drawing stops at the first failure. So a LOW person holds the INFORMAL flag with probability p_T × p_I.
- **The maximum number of relationships** of each type is the parameter value randomly rounded (2.3 gives 2 with probability 0.7 and 3 with probability 0.3).

| | LOW | MEDIUM | HIGH |
|---|---|---|---|
| Prob. extra, TRANSITORY M / F | 0.375 [0.266–0.453] / 0.109 [0.046–0.180] | 0.497 [0.328–0.678] / 0.526 [0.274–0.741] | 1 |
| Prob. extra, INFORMAL M / F | 0.018 [0–0.131] (quirk) / 0 | 0.182 [0–0.508] / 0.604 [0.263–0.797] | 1 |
| Prob. extra, MARITAL | 0 | 1 | 1 |
| Prob. extra, COMMERCIAL | 1 | 1 | 1 |
| Max TRANSITORY M / F | 1.81 [1.40–2.0] / 1.70 [1.33–2.0] | 3.24 [2.24–3.5] / 3.21 [2.36–3.5] | 1 |
| Max INFORMAL M / F | 0.62 [0.60–0.86] / 0.81 [0.60–1.13] | 1.72 [1.39–2.09] / 1.26 [1.00–1.53] | 1 |
| Max MARITAL M / F | 1 | 0.96 [0.90–1.20] / 0.98 [0.90–1.25] | 1 |
| Max COMMERCIAL | 0 | 0 | 59 |

**Relationship parameters** (the same in every node except the condom maxima):

| Type | Formation rate (/day) | Duration (Weibull scale, heterogeneity) | Coital acts/day | Condom use: min → max, mid year, rate |
|---|---|---|---|---|
| TRANSITORY | 1.37e-4 [1.37e-4–4.99e-4] | 0.957 y, 0.833 | 0.33 | 0 → node max, 2005 [2000–2010], 2.60 [1.82–3.43] |
| INFORMAL | 2.73e-4 [9.1e-5–8.8e-4] | 2.03 y, 0.967 [0.802–0.9999] | 0.33 | 0 → node max, 1996 [1992–2002], 1.80 [1.07–2.42] |
| MARITAL | 1.15e-4 [7.3e-5–1.5e-4] | 22.15 y, 0.667 | 0.33 | 0 → 0.170 [0.104–0.238], 1997 [1991–2002], 2.35 [1.30–3.24] |
| COMMERCIAL | 0.15 | 7 days, 1 (exponential) | 0.0027 | 0.5 → 0.85, 1999.5, 1 |

- **Regional condom maxima**:
  - Transitory: Northern 0.767 [0.657–0.8], Central 0.616 [0.531–0.722], Southern 0.791 [0.674–0.8].
  - Informal: Northern 0.574 [0.477–0.6], Central 0.527 [0.445–0.588], Southern 0.474 [0.401–0.543].
- **Risk assortativity** (rows are men, columns women; order LOW, MEDIUM, HIGH): [[v, 1−v, 0], [1−v, v, v], [0, v, 1−v]], with v = 0.610 [0.390–0.811]. It applies to transitory, informal and marital relationships; commercial relationships are not assortative. So LOW women never partner HIGH men, and LOW men never partner HIGH women.
- **Sexual debut**: Weibull with a minimum age of 13. Women: scale 16.85 [16.08–18.04], heterogeneity 0.086 [0.01–0.187]. Men: scale 16.65 [16.00–17.52], heterogeneity 0.058 [0.01–0.113].
- **Outputs driven**: prevalence by age and sex, and incidence. This is the main calibrated block (36 of the 55 parameters).

### 4. Transmission

Each step, every relationship (handled from the man's side) has a Poisson number of coital acts (`Relationship.cpp:292`):

- The mean is dt × coital act rate × a dilution factor. Dilution is 0.75 / 0.6 / 0.45 when the partner with more relationships has 2 / 3 / 4+ relationships.
- **The first act of every relationship is forced.** With a monthly step and 7-day commercial relationships, almost every commercial relationship has exactly one act.

The infection probability per act (`IndividualSTI.cpp:379,645`) is the product of:

- **Base_Infectivity**, 0.00129 [0.00073–0.0022]. It is defined female-to-male.
- **Stage**: × 26 during acute infection (2.9 months), × 10 during AIDS (the last 9 months of the prognosis), × 1 during latent infection.
- **ART**: × 0.04 for suppressive ART or × 0.65 for non-suppressive ART, reached linearly over 183 days. Back to × 1 (and AIDS × 10) when ART is failing (section 9).
- **Condom**: × (1 − 0.9) if a condom is used. Use is drawn per act from the relationship type's sigmoid.
- **STI co-infection**: × 5.5 if either partner has the STI flag. The larger of the two multipliers is used, never both.
- **Circumcision**: × (1 − 0.6) for a circumcised man.
- **PrEP**: × (1 − 0.75).
- **Male-to-female**: × the multiplier for women, by her age: `Male To Female Young` = 4.49 [3.65–5.0] up to age 15, then linear to `Male To Female Old` = 2.27 [1.21–3.0] at 25 and above.

Other settings:

- There is no superinfection (`Max_Individual_Infections` 1) and no individual heterogeneity in infectiousness (lognormal scale 0).
- **Outputs driven**: incidence (UNAIDS falls to about 20k a year; the model stays around 150k) and prevalence.

### 5. Natural history

- **Prognosis** (time from infection to HIV death without ART) is drawn at infection (`FerrandAgeDependentDistribution.cpp`):
  - Age 15 and over: Weibull with shape 2 and scale 21.182 − 0.2717 × age (age capped at 50). At age 20 the median is 13.1 y; at age 40 it is 8.6 y.
  - Under 15: 57% are rapid progressors (exponential, rate 1.52/y); the rest follow a Weibull with median 16 y and shape 2.7.
- **Stages** (`InfectionHIV.cpp:181`): acute for 2.9 months, AIDS for the last 9 months, latent in between. These are the only stages that affect transmission.
- **CD4 is deterministic given the prognosis** (both heterogeneities are 0; `SusceptibilityHIV.cpp:132`):
  - √CD4 falls linearly from √540.55 at infection to √31.63 at death.
  - So the ART thresholds sit at fixed fractions of the prognosis: CD4 500 at 5%, 350 at 26%, 250 at 42% and 200 at 52%.
- **WHO stage** comes from Weibull-drawn fractions of the prognosis (constants built into EMOD; values not checked). It is used only in CD4-agnostic ART staging.
- **Symptoms**: `NewlySymptomatic` fires a Weibull time before death (scale 618 d, heterogeneity 0.5; mean about 550 d). It triggers symptomatic testing.
- **Seeded infections start partway through their course.** `OutbreakIndividual` with `Incubation_Period_Override` = 0 fast-forwards each seeded infection to a uniform random point in its prognosis (`InfectionHIV.cpp:271`).
- **Outputs driven**: HIV deaths, CD4-based eligibility, and the diagnosis rate through symptomatic testing.

### 6. STI co-infection

- **This is a lifelong flag, not a disease.**
  - On day 1 (1960.5), 30% of post-debut MEDIUM/HIGH people and 10% of LOW people get the flag.
  - Afterwards the same probabilities apply at sexual debut, using the person's `Risk` at debut.
  - It is never cleared, and later becoming an FSW doesn't change it.
- **Effect**: × 5.5 on the per-act probability (section 4).

### 7. Testing and diagnosis

All of these pathways are switched off for anyone whose `InterventionStatus` is in the event's `Disqualifying_Properties`, which usually means LostForever or anyone already in staging, linkage, pre-ART or on ART. All tests are perfect (sensitivity and specificity 1).

- **Voluntary testing (HCT): uptake loop**, from 1990:
  - At debut, the person enters the testing loop with probability sigmoid(−0.005 → 0.05, mid 2005, rate 1).
  - Otherwise they wait an exponential delay (mean 365 d) and then enter with probability sigmoid(−0.01 → 0.15, mid 2006, rate 1), repeating.
  - People already past debut in 1990 start in this loop.
  - The probabilities are the same for both sexes.
- **HCT testing loop**: a rapid test after an exponential delay (mean 365 d), repeated while negative.
  - **The dropout state (`HCTTestingLoop2`) is never triggered.** Once in the loop, a person tests about once a year for life.
  - A positive test goes to "Consider immediate ART" (section 8).
- **Symptomatic testing**: at `NewlySymptomatic`, the person presents with a probability that is stepwise by year: 0.14 (1990), 0.61 (2004), 0.82 (2006), 0.94 (2010), 0.99 (2015).
  - Presenters go to ART staging.
  - Non-presenters get a second chance of 0.5 from 2016.
- **ANC** (from 2001): at 12 weeks of pregnancy, a woman is tested with probability by year, then gets a rapid test.
  - **The intended uptake table is ignored.** The intervention is `HIVSigmoidByYearAndSexDiagnostic`, which has no `Time_Value_Map`, so its default sigmoid applies (0 → 1, mid 2000, rate 1). That gives 0.73 in 2001 and 0.95 in 2003, against the intended 0.31 and 0.55; the two agree from 2004.
  - A positive result goes to ART staging (probability rising to 1 by 2007) and to PMTCT (section 11).
- **Infant testing**: at 6 weeks of age (from 2004), any child is tested with probability 0 → 0.34 (2009). Positives go to ART staging.
- **Dead code**: `HCTUptakePostDebut7` (events 81–82, "increased testing uptake after 2009") and `OnART1_defaulter_tracer` are never broadcast.
- **Outputs driven**: share diagnosed, tested in the past year, and timing of ART start.

### 8. Linkage, staging and pre-ART

- **"Consider immediate ART"** (HCT positives only, events 83–84): the person goes straight to ART initiation, skipping CD4 staging, with probability:
  - women: 0 before 2007, then 1;
  - men: 0.75 from 2006, 0.9 from 2012, 0.95 from 2018.
  - Otherwise they go to ART staging.
  - This is test-and-treat for people found by HCT from 2006/2007, years ahead of Malawi's policy (2016).
- **ART staging** (symptomatic, ANC, infant and pre-2007 HCT positives):
  1. A rapid test.
  2. 5% are lost after diagnosis (back to the HCT uptake loop).
  3. Blood draw, which records CD4; staging uses the lowest CD4 recorded.
  4. CD4-agnostic eligibility check: WHO stage ≥ 4 from 2004, ≥ 3 from 2008, everyone from 2016; pregnant women from 2004; children under 1/2/5 y from 2008/2011/2014.
  5. CD4 eligibility check: ≤ 200 / 250 / 350 / 500 from 2004 / 2008.3 / 2011.6 / 2014.3; pregnant ≤ 1000 from 2011.
  6. **Nobody is eligible before 2004.**
- **Linkage to ART** for eligible people:

| Parameter | Value |
|---|---|
| Linkage probability | sigmoid(0 → `ART Link Max`, mid `ART Link Mid`, rate 1) |
| `ART Link Max` | 0.947 [0.697–1.0] |
| `ART Link Mid` | 2009 [2002–2012] |

  - With a median mid-year of 2009, linkage is near 0 before about 2006. **This is the main calibrated control on the number on ART.**
  - Failure returns the person to the HCT uptake loop.
- **Pre-ART** (not yet eligible):

| Parameter | Value |
|---|---|
| Linkage to pre-ART | sigmoid(`PreART Link Min` → `Max`, mid `PreART Link Mid`) |
| `PreART Link Min` | 0.447 [0.112–0.707] |
| `PreART Link Max` | 0.696 [0.369–0.941] |
| `PreART Link Mid` | 2003 [1999–2006] |

  - In pre-ART, every 182 days 25% are lost (back to the HCT uptake loop) and 75% are re-staged.
  - Eligible people go **directly** to ART initiation, without the ART linkage sigmoid.
- **Outputs driven**: number on ART (fitted), ART coverage (fitted), and people diagnosed but not on ART.

### 9. ART

- **Initiation**:
  - A delay of U(1, 120) days.
  - Then a random choice between suppressive ("Effective") and non-suppressive ART: P(suppressive) = 0.63 for men under 35, 0.849 for men 35+, 0.746 for women under 35, 0.891 for women 35+.
  - The choice is fixed until dropout. There is no switching and no viral-load monitoring.
  - Both arms are `ARTMortalityTable`.
- **Effects**, in both arms (`HIVInterventionsContainer.cpp:491`):
  - Natural disease progression is paused at ART start.
  - CD4 recovers (+15.6 m − 0.21 m² after m months, saturating at +287 after 3.1 y).
  - Transmission is multiplied by 0.04 or 0.65, reached linearly over 183 days.
  - Mother-to-child transmission is × 0.0333.
  - So "non-suppressive" ART differs from suppressive ART only in its transmission multiplier and its mortality table.
- **Mortality on ART** (`ARTMortalityTable.cpp:111`):
  - The time to death is drawn at initiation from annual rates, by ART-duration bin (<182 d, 182–365 d, 1–2 y, 2–3 y, 3 y+) × age (<25, 25–35, 35–45, 45+) × CD4 at initiation (interpolated over 0–624.5), with separate tables by sex and arm.
  - **The draw for each later duration bin is timed from ART start, not from the start of that bin.** So the time to death is not the piecewise-exponential distribution the table implies.
  - In the last 9 months before that death, ART "fails": the person moves to the AIDS stage with no ART effect on transmission, i.e. × 10 (`HIVInterventionsContainer.cpp:109`).
- **Dropout** (events 74–80):
  - Exponential with mean 2051 d (5.6 y), the same for everyone. Dropout removes ART; progression resumes from the point reached at ART start.
  - Then 5% are lost forever (`InterventionStatus:LostForever`, excluded from everything).
  - The other 95% go to "Consider not dropping", which sends the person to restart ART with probability 0 before 2009, 0.85 from 2009, 0.9 from 2016 and 0.95 from 2018; otherwise they go back to the HCT uptake loop.
- **The restart branch is broken.** See the first problem below. From 2009, most people who drop out stay off ART and untested.
- **Outputs driven**: number on ART, ART coverage, viral suppression, HIV deaths, and transmission from people who interrupted treatment.

### 10. Male circumcision

- **Traditional**:
  - At 1975 every living person, and from 1975 every newborn, is circumcised with probability 0.025 / 0.0605 / 0.189 (nodes 1/2/3). It happens at birth, not in adolescence.
  - The event is applied to both sexes but has an effect only on men.
- **VMMC** (`NChooserEventCoordinatorHIV`, from 2011):
  - A target number of HIV-negative, uncircumcised men aged 5–40 per node and year (2011–2024 programme data, then flat to 2040).
  - Targets are multiplied by 0.01 and **truncated** (`NChooserEventCoordinator.cpp:346`). So Northern 2011 (54 men) gets zero, and every target loses up to one agent.
  - Each year's target is spread evenly over its steps, with the remainder given to the earliest steps.
  - The exclusion uses the intervention name `Any_MC`. Traditionally circumcised men don't carry that name, so they can be chosen. `MaleCircumcision` then refuses to distribute to them; whether the NChooser still counts them toward the target was not checked.
- **Efficacy**: 60% reduction in acquisition for both.
- **Outputs driven**: incidence in men; circumcision coverage (not reported in our runs).

### 11. Vertical transmission and PMTCT

- **Transmission** happens at birth only, with probability 0.3 (there is no separate breastfeeding route). It is × 0.0333 if the mother is on ART (either arm), otherwise × (1 − PMTCT efficacy).
- **PMTCT** after a positive ANC test:
  - sdNVP (efficacy 0.66), or combination with probability sigmoid(0 → 1, mid 2008.4).
  - Combination is always Option B (0.967).
- **Efficacy persists**: a woman's PMTCT efficacy is never reset after the birth (`HIVInterventionsContainer.cpp:632`), so it seems to carry over to later pregnancies. Not tested.
- **Outputs driven**: infections in children, and prevalence in under-15s.

### 12. Epidemic seeding

- At `SeedYr HIGH` = 1983 [1978–1991], a fraction `Demographic Coverage` = 0.020 [0.001–0.044] of HIGH-risk people (current FSW and clients) are infected in every node.
- The infections are fast-forwarded to random points in their course (section 5).

### 13. PrEP

- **Oral PrEP** from 2019:
  - Fixed numbers of HIV-negative women aged 15–35 (2020: 15–65) with `Risk` MEDIUM or HIGH (`NChooser`).
  - The targets are 1060 in 2019, then 219, 6472, 31,026, 76,122 and 106,858 (2024), then 5.3M over 2025–2076, all × 0.01.
- **Effect** is `ControlledVaccine`: 75% acquisition blocking for 90 days, then it expires, with no re-uptake.
  - The NChooser excludes women with an intervention named `PrEP`, but the intervention is called `Oral_PrEP`, so women already on PrEP are not excluded (minor).
- **Outputs driven**: negligible before 2025 (about 10–1000 agents).

### 14. Reporting

- **`ReportHIVByAgeAndGender`**: yearly rows by node, sex, age bin and IP. Columns: population, infected, newly infected, on ART, died, died from HIV, tested in the past year or on ART, tested ever, diagnosed, newly tested positive or negative.
- Daniel's setup reported from 2000 with bins 0/15/20/…/50. `run_emod.py` widens this to 1990, 5-year bins, and splits by `Risk` and `ARTStatus`. That changes the output only.
- `InsetChart`, `BinnedReport` and `DemographicsSummary` are also written.

## Switched off or unused

- **Migration** and between-node relationships; behavioural super-spreaders (0); individual heterogeneity in infectiousness; TB; `Min_Days_Between_Adding_Relationships` (0).
- **Accessibility**: 100% `Yes` in both the demographics file and the overlay (the overlay comment says 80%). So the 2010–2012 `ReferenceTrackingEventCoordinator` that raises it to 100% has no effect.
- **IPs defined but never used**: `Uptake` (75/25), `OnDrug`, `PrEPWilling`, `OnPrEP_CAB`, `WillingToStartARTFromPrEPCampaign`.
- **Event names declared but never used**: many of the `Custom_Individual_Events`, e.g. CAB-LA PrEP, PEPFAR disruption, vaccines and boosters.
- **Unused campaign markers**: `__KP_` markers (90-90-90 scenario hooks, symptomatic linkage, HCT post-debut, PrEP waning) that no calibrated parameter writes to.
- **Config parameters not in the 2.37 schema**, and therefore ignored: `ART_Viral_Suppression_Multiplier` (0.08; the intervention parameters set the ART multipliers), `Infectious_Period_*`, and `Circumcision_Reduced_Acquire` at config level.

## Problems found

1. **The ART restart branch traps people who drop out from 2009 onwards.**
   - Dropout (`ARTDropout`) stops ART and sets `ARTStatus` to OffART, but leaves `InterventionStatus` at `OnEffectiveART` or `OnNonSuppressiveART`.
   - "Consider not dropping" then sends 85–95% to `OnART1`. The suppressive/non-suppressive choice at `OnART1` lists those two statuses as disqualifying, so no ART is given. `BaseIntervention::Distribute` checks `Disqualifying_Properties`: `campaign/Interventions.cpp:147,178`.
   - Only the dropout timer restarts. Every testing and linkage event is also disqualified for these statuses.
   - Each 5.6-year cycle, only 0.95 × (1 − p), i.e. 5–14%, escape to the HCT loop, and 5% are lost forever. Before 2009 (p = 0) everyone who drops out returns to the HCT loop.
   - The output fits this:
     - OffART rises from 6% of PLHIV in 2010 to 46% in 2024 (set 0).
     - In 2020 only 354 of 18,921 OffART PLHIV agents (1.9%) were tested in the past year.
   - This is probably the main reason ART coverage levels off at about 40%. How much the mortality bug adds is not known.
2. **Background mortality is 1000 times too low** (known): `ResultScaleFactor` in `Demographics_Malawi.json`.
3. **The male non-suppressive ART mortality table has one row about 7 times too high.** For ART duration 182–365 days and CD4 ≥ 25, it holds 0.56 / 0.31 / 0.155 / 0.071 / 0.047 per year. The female row is 0.085 / 0.044 / 0.022 / 0.010 / 0.007, and the neighbouring male rows are 2–3 times the suppressive values, not 11 times. This looks like a typo.
4. **Two calibrated condom parameters have no effect** (`Trns Condom Max`, `Infmrl Condoms Late`), and **one calibrated concurrency parameter is never used** (`Pr Ex Infmrl Male LOW`). See Calibration.
5. **The regional LOW-risk fractions are applied to the wrong nodes** because of the node labels in the Risk overlay.
6. **The ANC testing uptake table is ignored** (wrong intervention class for a time table), so ANC testing is higher than intended in 2001–2003.
7. **The draw of time to death on ART restarts from ART start for each duration bin** (EMOD behaviour, `ARTMortalityTable.cpp:111`).
8. **Small-number effects at 1% scale**: VMMC and PrEP targets are truncated per node, and PrEP is about 10 agents in 2019.
9. **Test-and-treat for HCT positives from 2006/2007** is a modelling choice, not Malawi policy; the event is labelled "after 2009?".

## Implications for EMOD′ (Starsim v4 / STIsim)

EMOD′ copies the problems above, so it reproduces the published model (decided 3 Nov; see `docs/prompts.md`, P1).

**Probably easy:**

- **Demography**: age- and year-specific pregnancy and death tables, with no migration. The three nodes are separate populations, so a node attribute with partnerships restricted within a node, or three sims, would do.
- **Natural history**: an age-dependent Weibull prognosis, deterministic √CD4 decline, and stage multipliers on fixed windows. This is a small coded module.
- **Transmission**: per-act probability × Poisson acts per relationship with a forced first act. The multipliers are all simple. Note that EMOD multiplies probabilities, whereas Starsim multiplies hazards (see the v4 notes in `docs/prompts.md`).
- **STI flag, circumcision, PrEP, PMTCT, seeding**: simple; STIsim's `VMMC` and `Prep` already take numbers to distribute.
- **Reporting**: a custom analyzer by node, sex, age and IP.

**Probably hard:**

- **The PFA network.** STIsim's `StructuredSexual` differs substantially. A custom network module would be needed for:
  - queue-based matching with 20 × 20 joint age tables;
  - entry rates frozen after the 1975 burn-in;
  - correlated extra-relational flags, randomly rounded concurrency limits redrawn when `Risk` changes, and Risk-weighted partner choice;
  - coital dilution, the forced first act, and 7-day commercial relationships inside a monthly step.
- **The cascade as a state machine.** About 60 events with `InterventionStatus` disqualification rules and zero-delay event chains, which EMOD resolves within one step because newly given interventions are updated in the same loop (`InterventionsContainer.cpp:212`). Reproducing the cascade means reproducing those semantics, including the dropout trap and the dead branches. A declarative transition table keyed on status would probably fit better than coded interventions.
- **ART mortality table**: the draw rule above and the 9-month failure before death.
- **Random numbers**: EMOD's AES-counter streams can't be matched, so comparisons have to be distributional across seeds. P1's criterion already is.

**Decisions needed before P1:**

- Reproduce or fix the dropout trap, the mortality bug, the male mortality-table row and the ANC table. Each changes the fitted model.
- Whether EMOD′ keeps the 1% population scale, including its truncation effects.
