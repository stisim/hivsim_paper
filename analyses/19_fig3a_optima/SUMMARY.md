# 19: Optima and Optima′ against Figure 3A (Malawi) — summary

**Answer:** neither our Optima rerun nor Optima′ reproduces the published Optima HIV panel of Figure 3A (Bansi-Matharu 2025), in shares or totals. Optima′ matches our Optima, so the gap is between our project (saved by Optima 2.11.4, September 2022) and the one the paper used, not in how transmissions are counted. Until MIHPSA's Optima Malawi project arrives (asked of Rowan), "Optima" in the core question means our 2022 project and Optima′. Against the published panel in 2020 that stand-in is within ~3 points on undiagnosed, has ~12 points more "diagnosed, never treated", and ~6 points less "on ART". Optima's large "diagnosed, never treated" band comes from people diagnosed but not yet linked to care, not people waiting in care for ART.

## 1. Method

- **What Optima records.** Its force of infection is a sum over source compartments (cascade state × CD4 × population), and `model.py` already stores infections by source state and population (`raw_incibypop`). Mother-to-child transmissions are spread over the mother's states in proportion to the number of people in them.
- **Ever-treated tracer.** `optima_attr.py` patches a copy of `models/optima/reference/model_dt.py` at import (`~/src/optima` untouched). It counts, in each in-care and lost compartment, the people who have ever been on ART, moved by the same operations that move Optima's people (step transitions, ageing and risk moves, population forcing, the treatment and care rules). Results are bit-for-bit identical to the unpatched model, and the total reproduces Optima's `numinci`. In 2020, 24,101 of the 25,310 people in "lost" had been on ART, and 508 of the 4,014 in care.
- **Mapping to the paper's categories (D8):** undiagnosed → undiagnosed; diagnosed not in care, plus never-treated people in care or lost → diagnosed never treated; on ART, suppressed or not → on ART; ever-treated people in care or lost → interrupted. The simple version (all of lost = interrupted) differs by at most 0.5 points.
- **Counted:** sexual transmissions are the primary measure (the paper's network models count only sexual); all transmissions including mother-to-child are reported separately. Annual figures use Optima's convention: the rate in the step starting each whole year. The suppressed transmit nothing in this project (`efftxsupp` = 1, U=U).
- **Optima′:** an output only. `Transmission` records each step's infections by the source's state, and `infections_by_source(sim)` scales them to the solver's actual infection flows; dynamics unchanged. Starsim copies the disease but not the rate functions, so `Transmission` is shared rather than copied (`__deepcopy__` returns itself).

## 2. Optima′ against our Optima (sexual), 2000–2020

- Undiagnosed within 1.7%, diagnosed-not-in-care within 2.0%; lost up to 2.7%, unsuppressed on ART up to 2.6%, total up to 2.55% (2020).
- In care differs by up to 85% at single whole-year steps (2018, 2020); it carries under 2% of transmissions, and the treatment rule empties it each step, so a one-step snapshot is unstable.
- Shares in 2020 within 1.2 points (below).

## 3. Against the published panel (`common/fig3a_malawi.csv`)

| 2020, % of new infections | Undiagnosed | Diagnosed, never treated | On ART | Interrupted |
|---|---|---|---|---|
| Published Optima HIV | 45.2 | 28.2 | 17.6 | 9.0 |
| Our Optima, sexual | 42.3 | 40.7 | 11.4 | 5.7 |
| Our Optima, incl. mother-to-child | 38.0 | 36.9 | 20.0 | 5.1 |
| Optima′, sexual | 42.7 | 39.4 | 11.9 | 5.9 |

In 2024 (a projection for us, inputs held at 2020): 39.4 / 37.2 / 17.4 / 6.0, against the published 50.5 / 5.3 / 35.6 / 8.5.

**Totals** (Figure 1A saved as `resources/bansi-matharu2025_fig1.jpg`; published values read by eye):

| New infections, ages 15–64 | Ours | Published Optima (Figure 1A) |
|---|---|---|
| 2000 | 61k | ~71k |
| 2010 | 52k | ~45k |
| 2015 | 48k | ~45k |
| 2020 | 30k | ~22k |
| 2024 | 25k | ~7k |

At all ages our 2024 total is 29,191, against the caption's 8,459.

## 4. Why they differ

- **The published run is a different calibration.** The model structure looks the same (Figure 1A shows 65% of people on ART suppressed in 2008, the same step-size artefact as in our project), but outputs differ from 2000 on: in 2000 the published run has 70% undiagnosed (ours 55%), 3.7% interrupted when almost nobody had started ART (ours 0.2–0.4%), and ~5% on ART (ours 1%); new infections differ throughout. MIHPSA must have used a later project with treatment inputs extended past 2020.
- **The published 2023–24 changes need inputs we don't have.** The published "diagnosed, never treated" share falls from ~30% to 5% while "on ART" rises. Our project holds the treatment and suppression proportions and the linkage time fixed from 2020, so our shares barely move after 2020.
- **Where Optima's "diagnosed, never treated" band comes from:** not people waiting in care (in 2020, in care holds 4k people causing 0.6k infections) but people diagnosed and not yet linked to care (123k people causing 11.7k infections). Optima's `linktocare` is a mean time to linkage, and the project's value rises from 0.1 to 3.16 years by 2020: this is how the calibration squares 90% diagnosed with a forced number on ART. The published run has the same band, smaller (28% in 2020).
- **Why our "on ART" share is low in 2020 (11% sexual):** only unsuppressed people on ART transmit, and their number falls from 2015 to 2020 as viral-load monitoring rises.

## 5. What this means

- Until we have MIHPSA's Optima Malawi project, the Optima end of the core question is our 2022 project and Optima′, with the differences above.
- `run_optima.py` should rerun on that project file once it arrives. The tracer only covers projects older than Optima 2.12.0 run at the calibrated step (dt = 0.2), and stops with an error otherwise.

**Files:** `README.md`, `optima_attr.py`, `run_optima.py` (env `optima`), `run_prime.py` (env `hp_rc3`), `compare.py`; `outputs/optima_by_source.csv`, `optima_new_infections.csv`, `prime_by_source.csv`, `shares_vs_published.csv`, `compare.txt`, `run_prime.log`; `figures/fig3a_optima.png`. Changed: `models/optima/v4/optima_prime.py` (output only).
