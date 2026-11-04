# 02: EMOD′ sexual network against EMOD — summary

**Answer:** EMOD′'s network and demography reproduce EMOD's within its run-to-run noise. On the same parameter sets without HIV, every statistic compared agrees with EMOD to within a few percent from 1961 to 2025: population by age, sex and node; the share in each relationship type by sex and age; concurrency; formation rates; durations; coital acts; age differences; Risk mixing; and contacts per FSW and per client. The 5–7% commercial-sex excess in the first version was noise: it reversed sign with a change that leaves the mechanism unchanged. The frozen pair-formation rates match the tables EMOD logs, but because they are frozen, each run's commercial intensity varies by about ±15%. Against the published model (with HIV), the two diverge after about 1990, by as much as HIV deaths change the population, most of all for marital relationships.

## What matches

Solid lines and bands in the figures are medians and ranges over the 10 sets. Dashed lines are medians over sets 0, 100 and 200 for both models, so each pair of lines is like for like. Numbers are v4 against EMOD without HIV on those 3 sets unless stated.

- **Demography** (`demography.png`, the demography alone, 20 seeds):
  - Total population: −0.2% in 1984 against all 13 EMOD runs (107,757 ± 395 against 107,976 ± 453); +0.04% in 2023 against the 3 without HIV.
  - All 5-year bands by sex within about 1% in 2023; yearly births within about 0.5%.
  - The bands needed EMOD's own ageing: EMOD keeps age in single-precision days, so after whole steps a monthly cohort falls on the other side of some 5-year edges, and EMOD reports before the step's births. Both are copied.
  - Single network runs differ by −0.3% to +1.9% in 2023, which is seed noise.
- **In a relationship** (`in_relationship.png`, share of those past debut): transitory, informal and marital within −1.6% to +1.1% in 1980, 2000 and 2020; against all 10 EMOD sets in 1980, −0.3% to +0.2% for every type. EMOD's marital peak (28% around 1978) and decline to 25% by 2005 are reproduced.
- **By age** (`by_age.png`; 1985 over the 10 sets, 2020 over the 3): the shapes match.
- **Concurrency** (`concurrency.png`): −0.5% to +0.1%; mean number of partners −1.5% to +0.2%.
- **Formation** (`formation.png`, 1970–2024): transitory +1.8%, informal +2.3%, marital −0.3%.
- **Frozen rates**: EMOD's tables were logged at 1975.5 for all 10 sets (`raw_results/emod/rates_debug/`). v4 / EMOD: FSW commercial 1.005 ± 0.026 (s.e.), clients 1.009, transitory 0.992 ± 0.013, informal 0.983 ± 0.013, marital 1.031 ± 0.060.
- **Queues** (`queues.png`): match, including the runaway commercial queue of men (about 1,800 entries per 1000 past debut).
- **Durations and acts** (`durations.png`): drawn means 328 / 728 / 7300 / 7.0 days against 329 / 731 / 7302 / 7.0; actual durations and acts per relationship also match.
- **Age differences** (`age_difference.png`): means 3.3 / 4.6 / 4.6 / 6.8 years against 3.3 / 4.6 / 4.9 / 6.7.
- **Commercial contacts** (`commercial.png`), 1990–2024 means: per FSW 98.0, against 99.8 (10 sets with HIV) and 99.5 (3 without); per client 23.1, against 24.9 and 23.6. On the 3 matched sets alone v4 gives 92.2 per FSW; the first version gave 101.6 on the same sets, so this is the spread of frozen rates in single runs, not a bias.
- **Risk mixing**: within 0.9 points for every type.

## What doesn't, and likely causes

1. **With HIV, EMOD diverges from about 1990**: by 2020 it has marital −21% and concurrency −7% relative to v4, as expected with no HIV in this comparison (analysis 03 adds HIV). Up to the late 1980s v4 is within about 1% for every type.
2. **Commercial intensity varies a lot between runs** (frozen rates), so agreement can only be judged as an average over runs or sets.
3. **Small simplifications**: when a woman queued twice for one type is matched, the matched entry is removed rather than EMOD's most recent one; EMOD's 50,000-candidate cap isn't applied (it never binds).

## Caveats

- One v4 run per set against one EMOD run per set (3 without HIV).
- Snapshots are taken after the step's ageing and before its births, as EMOD's report is.
- v4's random numbers aren't safe for common random numbers.
- These outputs were produced before the STI flag was added to the risk module; that changes only the random-number sequence.

## Next steps

- HIV: analysis 03.
- For commercial sex, compare averages over many seeds, or over several EMOD runs per set.
