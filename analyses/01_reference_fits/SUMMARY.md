# 01: Summary (EMOD so far, 3 Nov 2026)

## EMOD

Rerun with the published `emod-hiv` 2.37.0 binary instead of Daniel's 2023 cluster build (setup in `models/emod/README.md`). Runs take about 13 minutes each.

**What matches**
- Prevalence by age and sex (DHS 2004, 2010, 2015-16; MPHIA 2015-16, 2020-21). Prevalence was 80% of his fit's weight.
- The number on ART (programme data to 2015), and the share of PLHIV diagnosed in 2020-21 (about 90% for women, 85% for men, as in MPHIA).

**What doesn't**
- Population is about twice the census, and PLHIV is about double (about 2M by 2020; UNAIDS about 1M).
- New infections stay around 150k a year (UNAIDS falls to about 20k), and AIDS deaths peak too high and rise again after 2015.
- ART coverage levels off at about 40% from 2014 (UNAIDS: 86% by 2021), and 33% of PLHIV are virally suppressed in 2020-21 (MPHIA: 80-86%).
- FSW prevalence is about half the survey values (2019: 26% against 50%).

**Cause of the population problem: background mortality about 1000 times too low**
- `MortalityDistributionMale/Female` in `Demographics_Malawi.json` hold annual deaths per person (infant mortality 0.18 in 1960), but their `ResultScaleFactor` is 1/(365*1000), the factor for deaths per 1000. Every rate is divided by an extra 1000.
- In the output, non-HIV deaths are about 0.03 per 1000 per year (Malawi: about 10-15), so the population grows about 4% a year.
- His runs use the same input files, so they almost certainly have the same problem. It didn't show up because his calibration rescaled population and on-ART counts to 8.26M aged 15-49 in 2018 and gave population only 2.5% of the weight.

**ART coverage: mostly dropout**
- In 2020 (set 0), 42% of PLHIV have `ARTStatus=OffART`: diagnosed, started ART and stopped. 33% are on suppressive ART, 8% on non-suppressive ART, and 17% have never started.
- The number on ART matches the programme data, so with PLHIV about double, coverage comes out about half. How much of the dropout is his calibrated parameters and how much the mortality bug isn't known yet.

**Check of the 2.37 build:** our scores agree only weakly with his likelihoods across the 10 sets (rank correlation 0.43, in the right direction); too few sets to say much.

**Open decision:** fix the mortality scale factor (one line, but his parameter sets were fitted without background mortality, so they would probably need recalibrating) or reproduce the model as published. Raised with Daniel on 3 Nov.

## Optima, HIVsim

Not yet measured against the common targets.
