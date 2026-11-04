# 04: EMOD′ against the published EMOD — summary

**Answer:** with the whole campaign interpreted from Daniel's campaign file, EMOD′ reproduces the published EMOD within its run-to-run noise, on the common targets and on the cascade over time, including the dropout trap that leaves 42–44% of people with HIV OffART by 2020.

## What matches

EMOD′ (10 sets × 3 seeds) against EMOD (the same 10 sets, one run each, `raw_results/emod/sub10`). The 10 sets are 8 distinct parameter vectors: Daniel's 250 resampled sets contain only 129 distinct ones (0 = 125, 25 = 50), told apart in EMOD by Run_Number, so EMOD′ seeds from the run number too.

Median ratio EMOD′ / EMOD over all target points (`outputs/targets_compared.csv`, `targets_compared.png`):

| Target | Median ratio | Range over target points |
|---|---|---|
| Diagnosed | 1.012 | |
| ART coverage | 1.015 | 1.00–1.12 |
| On ART among diagnosed | 0.990 | |
| Suppressed among on ART | 1.000 | |
| VLS among PLHIV | 1.007 | |
| Prevalence | 0.99 | 0.88–1.08 |
| Population | 1.00 | |
| People with HIV | 0.97 | 0.87–1.01 |
| New infections | 0.97 | 0.94–1.04 |
| Incidence | 0.95 | |
| AIDS deaths | 0.97 | 0.90–1.05 |

Over time (`cascade_over_time.png`, medians), 2020:

| | EMOD′ | EMOD |
|---|---|---|
| People with HIV (agents) | 22.4k | 23.5k |
| On ART | 8.9k | 9.5k |
| ART coverage | 40.0% | 40.1% |
| Diagnosed | 85.4% | 84.6% |
| OffART share of PLHIV | 43.5% | 42.0% |
| OnEffectiveART share of those on ART | 80.1% | 80.1% |
| New infections a year | 1,606 | 1,682 |
| HIV deaths a year | 1,094 | 1,114 |

- The diagnosed, coverage and OffART curves overlie each other from 2004.
- People with HIV in 2020 by set (median of 3 seeds against EMOD's run): within ±2% for 8 of 10 sets; −3% for set 50, +9% for set 25.
- The shared target figures for each model (`emod_coverage_*`, `v4_coverage_*`) show the same misfit to the data as analysis 01.
- EMOD′'s InterventionStatus (`v4_intervention_status.png`), 2020: about 63% of people with HIV OnEffectiveART and 16% OnNonSuppressiveART (most of both have dropped out and aren't on ART; this is the trap), 9% HCTUptakePostDebut, 4% LostForever.

## What doesn't, and likely causes

1. **EMOD′ runs 3–5% low on people with HIV, new infections and deaths** (medians over all runs), while per-set medians are close: EMOD′'s band has a lower tail of runs where node 1's epidemic doesn't start.
2. **Seeding** (`seeds.py`, `outputs/seeds.csv`), set 0 with 20 EMOD runs (Run_Number 1–20, to 1982) against 50 EMOD′ seeds: mean seeds in node 1 are 1.40 against 0.92, and none in 20% of EMOD runs against 38% of EMOD′'s. Nodes 2 and 3 agree (5.65 / 5.34, 8.10 / 9.14). Part is real: EMOD′ has 11% fewer Risk-HIGH people in node 1 in 1980, consistently across the 10 sets (ratio 0.889; nodes 2 and 3 within 1%), cause not found. Scaled by that deficit, EMOD′ would seed about 1.03 in node 1 against 1.40 ± 0.23, so the rest is plausibly chance.
3. **PMTCT, not checked against EMOD:** women found positive at ANC also enter ART staging, whose status disqualifies the combination-PMTCT branch, so after about 2008 almost nobody gets Option B (only single-dose nevirapine before that). This follows from EMOD's update order, which EMOD′ copies. Also unchecked: PMTCT efficacy is reset after 280 days in EMOD′, following the source, where the inventory says it's never reset.
4. **PrEP exclusion:** the NChooser excludes people with an intervention named `PrEP`, but the intervention is named `Oral_PrEP`, so nobody already on PrEP is excluded. Copied, as another of Daniel's quirks.

## Caveats

- EMOD has one run per set; EMOD′ three.
- EMOD's report here has no InterventionStatus, so the cascade is compared through ARTStatus, ART and diagnosis.

## Next steps

- Find the node-1 Risk-HIGH deficit.
- Run EMOD with InterventionStatus in the report; check the PMTCT branch and infant infections.
- Then the component swaps.
