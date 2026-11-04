# 06: The condensed cascade against the compiled campaign — summary

**Answer:** replacing Daniel's campaign (95 compiled deliveries) by the condensed cascade in `common/cascade/` (43 deliveries driven by a 49-row parameter table read from the same campaign, talking to HIV only through `CascadeHost`) changes EMOD′ by less than run-to-run noise on every measure compared: 343 of 345 target points within 2 standard errors (worst z = 2.5, one on-ART point inside the compiled 5–95% range), and target medians within 0.2% for people with HIV, infections, deaths and ART (prevalence +0.6%, incidence +2%). It runs in 74 s against 97 s. Starsim v4 at ebd52c1a (includes PR #1538).

## By measure (condensed vs compiled, 10 sets × 3 seeds, the same seeds)

- **2020 medians** (`cascade_over_time.png`): people with HIV 22.9k vs 22.8k; on ART 9,201 vs 9,125; diagnosed 84.5% vs 84.6%; ART coverage 40.3% vs 39.9%; OffART 42.3% vs 42.0%; suppressed among on ART 80.1% vs 80.1%.
- **Places in care** (`care.png`; compiled InterventionStatus grouped, with OnEffective/OnNonSuppressiveART split by ARTStatus into on ART and out): staging, pre-ART, linking, on ART, out and lost overlie within 0.1–0.3 points; the testing loop is 0.2 points higher after 2012 (6.0% vs 5.8% in 2020). None and uptake differ by up to 7 points in 2001–2004 because EMOD keeps people tested at ANC, on symptoms or at 6 weeks in transient TestingOn* statuses, which the condensed cascade doesn't have; a labelling difference, the dynamics agree.
- **Infant infections a year** (`infants.png`), 2005 / 2009 / 2012 / 2020: 171.7 / 234.7 / 163.1 / 187.8 vs 180.2 / 235.2 / 165.0 / 188.7. The PMTCT Option-B block, now a parameter (`pmtct.blocked_if_staged`), reproduces.

## The component

- `common/cascade/interface.py`: `CascadeHost`, what the cascade needs from an HIV model (people: age, sex, pregnancy, debut, dying; events: debut, first symptoms, 12 weeks pregnant, 6 weeks old; HIV: infected, test, CD4, WHO stage, on ART, `start_art(uids, arm)`, `stop_art`, PMTCT, circumcision, PrEP). `models/emod/v4/cascade_host.py` implements it for EMOD′.
- `common/cascade/params.py` and `tables/malawi_emod_set000_{params,flows}.csv`: one row per pathway parameter (forms: constant, sigmoid, step, linear, delays), each naming the campaign event it came from; VMMC and PrEP as flows per period. `models/emod/v4/extract_cascade.py` builds the tables for any parameter set (only the two linkage sigmoids differ between sets).
- `common/cascade/cascade.py`: pathways for routine testing, symptomatic, ANC and infant testing, staging and eligibility, pre-ART, ART, dropout and return, PMTCT, VMMC and PrEP; one place in care per person (`people.care`). A new testing pathway (e.g. FSW) is a table row plus an entry in `testing=`.
- EMOD's quirks are explicit parameters: the dropout trap (`dropout.p_restart = 0`), the Option-B block (`pmtct.blocked_if_staged = 1`), the ignored ANC table (`anc.p_test` as EMOD's default sigmoid).

## What was condensed

One eligibility rule for staging and the pre-ART recheck (EMOD's two tables differ slightly); EMOD's disqualifying properties reduced to one place in care per person; waits that end without effect when the person has moved on, rather than being cancelled; ordering interactions as parameters. Traditional circumcision stays in EmodHIV. None of these shows at the precision of 30 runs.

## How the interface would fit the other models

- **STIsim `sti.HIV` / `sti.ART`:** diagnosis, CD4, ART start (with the arm set by `p_effective_art`) and stop, PrEP and circumcision map directly. Missing: WHO stage, a first-symptoms event, PMTCT efficacy per mother, and choosing the ART arm outside `ART`.
- **Optima′:** compartmental, so a per-agent interface can't drive it; the same table would have to become transition rates between undiagnosed, diagnosed, in care, treated and lost.

## Next steps

- Send `common/cascade/tables/malawi_emod_set000_params.csv` to Daniel to check.
- Implement `CascadeHost` for STIsim's `sti.HIV` (WHO stage, a symptom event, per-mother PMTCT).
- Decide how the cascade enters Optima′ (rates from the table).
- Use the condensed cascade in the component swaps.
