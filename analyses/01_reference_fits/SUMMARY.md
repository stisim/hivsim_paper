# 01: Summary (EMOD and Optima, 3 Nov 2026)

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
- Most of this is a campaign bug, found by the EMOD inventory (`docs/emod_inventory.md`) and confirmed in the campaign file. Dropping out of ART sets `ARTStatus` to `OffART` (event 76) but leaves `InterventionStatus` at `OnEffectiveART` or `OnNonSuppressiveART` (set at ART start, events 66 and 68). From 2009, 85-95% of dropouts are sent to restart (event 79 to `OnART1`), but the four ART-start choices (events 62-65) are blocked for exactly those statuses, as are the testing and linkage pathways. So almost nobody who drops out restarts. Not checked against Daniel's 2023 build, but nothing in the campaign resets the status.
- How much of the low coverage is this trap and how much the doubled PLHIV from the mortality bug isn't known yet.

**Check of the 2.37 build:** our scores agree only weakly with his likelihoods across the 10 sets (rank correlation 0.43, in the right direction); too few sets to say much.

**Open decision:** fix the mortality scale factor (one line, but his parameter sets were fitted without background mortality, so they would probably need recalibrating) or reproduce the model as published. Raised with Daniel on 3 Nov.

## Optima

One deterministic run of the 2022 project (parset 20220912), measured by `models/optima/reference/optima_targets.py`. Every target has a value; none is NaN. Optima has no age bands for FSW, clients and MSM (all 15-49), so for targets by age they are spread over the 15-49 bands of their sex in proportion to each band's size. With one run there is no band, and the coverage table's "in range" counts are 0 by construction; read the figures instead.

**What matches**
- Population (WPP, within 1.5% in 1990-2020; Optima forces population size to its input).
- Prevalence by age and sex for ages 25 and over in all five surveys, mostly within the survey intervals. Men aged 30-55 are a little low in most surveys.
- PLHIV from 2005 to 2015 (within 3% of UNAIDS), and new infections from 2003 to 2012.
- ART coverage to 2017 (within about 5 points of UNAIDS), and the share of PLHIV diagnosed in 2020-21 (91% for women against 90%; 89% for men against 85%).
- FSW prevalence in 2006 (72% against 71%) and 2019 (40%, at the bottom of the survey interval).

**What doesn't**
- Prevalence at 15-24 is too high: 15-19 is 4.5% for women and 2.4% for men in 2020-21 (MPHIA: 1.7% and 1.4%), and two to three times the data in every survey. Most of it is in Optima's own 15-19 populations (3.9% and 1.7% without the key populations spread in).
- The early epidemic is too small: new infections peak at 87k (UNAIDS 130k), AIDS deaths at 56k (UNAIDS 76k), and PLHIV is 13-15% low around 2000. Deaths stay 20-40% low throughout, 39% low by 2020.
- New infections fall too slowly after 2012: 34k in 2020 (UNAIDS 21k), and 2020 is 40% of the early peak (UNAIDS 16%). PLHIV keeps rising to 1.15M in 2021 (UNAIDS flat at 1.0M).
- Incidence 15-49 in both MPHIA surveys is at or above the top of the survey interval (2020-21: 0.44% for women, 0.26% for men; MPHIA 0.31% and 0.15%).
- The cascade is short at the treatment step in 2020-21: 86% of diagnosed people are on ART (MPHIA 97-98%), so ART coverage is 79% for women and 77% for men (MPHIA 89% and 83%), and 73% and 71% are suppressed (MPHIA 86% and 80%). From 2018 ART coverage falls behind UNAIDS (78% against 84% in 2020). Optima holds the share in care on ART fixed from 2020.
- FSW prevalence is 9-13 points low in 2014 and 2018 (53% against 63%, 42% against 55%).
- The number on ART 15-64 matches the programme data to 2008, then is 10-30% low (2010: 150k women against 218k). Optima fits its own all-ages on-ART series exactly, so the two data sources differ.

## HIVsim

Not yet measured against the common targets.
