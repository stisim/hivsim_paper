# 04: EMOD′ against the published EMOD

**Question:** with the whole campaign (testing, care, ART, VMMC, PrEP and PMTCT, interpreted from Daniel's campaign file by `models/emod/v4/cascade.py`), does EMOD′ reproduce the published EMOD Malawi model, on the common targets and in its cascade states?

**What's run:**
- **EMOD:** the 10 published parameter sets of analysis 01 (`raw_results/emod/sub10`: sets 0, 25, ..., 225, one run each, to 2025, reporting by Risk and ARTStatus from 1990).
- **EMOD′:** `harness.make_sim(pid, seed, cascade=True)` for the same 10 sets, 3 seeds each. Its `EmodReport` writes EMOD's `ReportHIVByAgeAndGender` table, so both are measured by the same code (`models/emod/reference/emod_targets.py`'s `measure()`).
- **Seeding:** `seeds.py --emod` runs EMOD for set 0 with Run_Numbers 1-20 to 1982 (campaign reduced to what precedes the seeding; `raw_results/emod/seeds_set000`), and `seeds.py` counts the people seeded in each node by EMOD and by EMOD′ (50 seeds). This settles whether set 0's node 1 epidemic failing to start in EMOD′ (analysis 03) is chance.
- **PMTCT and InterventionStatus:** EMOD for sets 0, 100 and 200 with InterventionStatus in the report (`run_emod.py --sets 0 100 200 --out ip3 --report-ip InterventionStatus`), against EMOD′ as it is and in two variants made in `pmtct.py` for the test only: the combination-PMTCT decision not disqualified by ART staging (`unblocked`), and PMTCT efficacy not reset after 280 days (`persist`). Compared: InterventionStatus of people with HIV over time, state by state, and infections at ages 0-4. If EMOD gives Option B to women found at ANC, its infant infections fall below EMOD′'s from about 2009; if it doesn't, they follow EMOD′ as it is.
- `run.py` runs EMOD′, measures both, and writes `outputs/`; `run.py --plot` draws `figures/`.

**Compared:** the common targets (`common/targets.csv`: prevalence by age and sex, people with HIV, new infections, AIDS deaths, ART coverage, diagnosed, on ART among diagnosed, suppression, FSW prevalence, population), as distributions over sets and seeds; and over time, people with HIV by ARTStatus (including the share OffART, the dropout trap), on ART, diagnosed, new infections and deaths.

**Caveat:** EMOD's report here has no InterventionStatus, so the cascade is compared through ARTStatus, ART and diagnosis. EMOD′ also reports InterventionStatus, which is shown for EMOD′ alone.
