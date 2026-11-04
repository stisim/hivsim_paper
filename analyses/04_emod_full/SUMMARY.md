# 04: EMOD′ against the published EMOD — summary

**Answer:** with the whole campaign interpreted from Daniel's campaign file, EMOD′ reproduces the published EMOD within its run-to-run noise: on the common targets, on the cascade over time (including the dropout trap that leaves 42% of people with HIV OffART by 2020), and state by state in InterventionStatus. Two behaviours EMOD′ derives from the campaign were checked against EMOD and confirmed: combination PMTCT (Option B) is blocked for women found at ANC, and PMTCT efficacy is reset after 280 days.

## What matches

EMOD′ (10 sets × 3 seeds) against EMOD (the same 10 sets, one run each, `raw_results/emod/sub10`). The 10 sets are 8 distinct parameter vectors: Daniel's 250 resampled sets contain only 129 distinct ones (0 = 125, 25 = 50), told apart in EMOD by Run_Number, so EMOD′ seeds from the run number too.

- **Targets** (`outputs/targets_compared.csv`, `targets_compared.png`), median ratio EMOD′ / EMOD over target points:

  | Target | Median ratio | Range |
  |---|---|---|
  | Diagnosed | 0.999 | |
  | ART coverage | 0.998 | 0.98–1.01 |
  | On ART among diagnosed | 0.998 | |
  | Suppressed among on ART | 0.999 | |
  | VLS among PLHIV | 1.003 | |
  | Prevalence | 0.996 | 0.87–1.09 |
  | Population | 0.998 | |
  | People with HIV | 0.97 | 0.87–1.01 |
  | New infections | 0.99 | 0.93–1.04 |
  | Incidence | 0.99 | |
  | AIDS deaths | 0.98 | 0.90–1.07 |
  | On ART | 0.96 | (down to 0.66 in 2004–06, with a few hundred agents on ART) |

- **Over time** (`cascade_over_time.png`), 2020 medians:

  | | EMOD′ | EMOD |
  |---|---|---|
  | People with HIV (agents) | 22.3k | 23.5k |
  | On ART | 9.0k | 9.5k |
  | ART coverage | 39.8% | 40.1% |
  | Diagnosed | 84.5% | 84.6% |
  | OffART share of PLHIV | 41.9% | 42.0% |
  | OnEffectiveART share of those on ART | 80.1% | 80.1% |
  | New infections a year | 1,632 | 1,682 |
  | HIV deaths a year | 1,074 | 1,114 |

- **People with HIV in 2020 by set**: within ±2% for 8 of 10 sets (set 25 +10%, set 50 −5%).
- **InterventionStatus, state by state** (`intervention_status_compared.png`; EMOD sets 0, 100, 200 with InterventionStatus in the report, `raw_results/emod/ip3`, `run_emod.py --report-ip InterventionStatus`): every state overlies EMOD's from 1995 to 2025, except TestingOnSymptomatic before 2004 (below).
- **This needed EMOD's delay timing**: EMOD's countdown timers fire at the first update after their time is up (`CountdownTimer::Decrement`), so a 1-day muxer holds people in ARTStaging for a step and the 1–120-day delay to ART keeps them LinkingToART a step longer. EMOD′ now does the same.
- **PMTCT** (`pmtct.py`, `pmtct_infants.png`), new infections at ages 0–4 a year in 2009 / 2012 / 2020: EMOD 245 / 177 / 214, EMOD′ 242 / 170 / 218. Test variants of EMOD′: with Option B not blocked, 89 / 97 / 197; with PMTCT efficacy not reset, 171 / 137 / 192. So EMOD blocks Option B for women at ANC (ART staging sets a status that disqualifies the combination decision, from EMOD's update order) and resets PMTCT after 280 days. The inventory's "never reset" is wrong for this binary.
- **Seeding** (`seeds.py`, `outputs/seeds.csv`; set 0, EMOD Run_Numbers 1–20 against 50 EMOD′ seeds): mean seeds by node 1.40 / 5.65 / 8.10 in EMOD against 0.92 / 5.34 / 9.14; no seed in node 1 in 20% of EMOD runs and 38% of EMOD′'s. The number of Risk-HIGH people at seeding agrees by node and sex (node-1 men 71.3 ± 1.6 against 68.9 ± 1.1), and both node-1 seed means are within about 1.4 standard errors of the expected 1.1, so the difference is chance. (The 11% fewer node-1 HIGH people seen earlier came from particular runs.)
- **PrEP exclusion**: the NChooser excludes people with an intervention named `PrEP`, but the intervention is named `Oral_PrEP`, so nobody already on PrEP is excluded. Copied, as one of Daniel's quirks.

## What doesn't, and likely causes

1. **TestingOnSymptomatic before 2004** is 12% lower in EMOD′ (4.2% against 4.8% of PLHIV). Not found; possibly when EMOD draws the time from symptoms to death (per person, not per infection).
2. **People with HIV run 3% low over all runs** while per-set medians are close: EMOD′'s lower tail is runs in which node 1's epidemic doesn't start, which is chance in the seeding.

## Caveats

- EMOD has one run per set; EMOD′ three. The InterventionStatus check uses 3 sets.
- The `unblocked` and `persist` variants exist only in `pmtct.py`, for the test.

## Next steps

- TestingOnSymptomatic before 2004.
- The component swaps.
