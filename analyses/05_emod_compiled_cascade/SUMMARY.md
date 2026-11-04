# 05: EMOD′ with its campaign compiled into Starsim deliveries — summary

**Answer:** compiled into 95 native Starsim deliveries, Daniel's campaign gives the same results as the literal interpreter, in distribution, on every measure compared: all 345 target points are within 2 standard errors of each other, the cascade over time and every InterventionStatus state overlie, and infant infections match, so the PMTCT Option-B blocking and the 280-day reset reproduce. The compiled run takes about 1.6 times as long (96 s against 60 s). It needs Starsim PR #1538 (cycles through `delay=`).

## What compiles natively, and what doesn't

From `outputs/objects.csv` and `models/emod/v4/compile_campaign.py`:
- **95 deliveries**, 9 with `delay=`, from about 60 campaign events (dead branches pruned).
- **Native:** listeners as event triggers with conditions on people (node, sex, age, Disqualifying_Properties on `istatus`); probabilistic decisions as `coverage=` with the complement on `trigger & ~yes`; random choices as conditional chains; muxers as `delay=` (Starsim's wait, max(1, ceil(d/dt)) steps from the trigger, equals EMOD's countdown timing, and every muxer here has Max_Entries = 1); the 1990 initialization as `at=`; PMTCT's reset as `delay=`; cycles through delays.
- **Coded products:** `SetIP` (54, status and properties), `Decide` (10, decisions on the person's state), `Act` (12: ART, dropout, PMTCT, circumcision, PrEP), `MuxEntry` (6), `ResetPMTCT` (3).
- **`deliver=` adapters:** 4 for events from coded modules (STIDebut, NewlySymptomatic, SixWeeksOld, TwelveWeeksPregnant), and the NChoosers for VMMC and PrEP.
- **A state holder** that marks who is dying this step and ends muxer waits on disqualification.

## What matches

Literal against compiled, the same 10 sets × 3 seeds:
- **Targets** (`targets_compared.png`, `outputs/targets_compared.csv`), compiled median over literal: people with HIV 1.001 (0.99–1.01), new infections 1.000 (0.97–1.04), AIDS deaths 1.000 (0.97–1.03), ART coverage 0.999, diagnosed 0.998, suppressed among on ART 0.999, prevalence 0.999 (0.93–1.06). All 345 points within 2 standard errors.
- **Over time** (`cascade_over_time.png`), 2020 medians:

  | | Compiled | Literal |
  |---|---|---|
  | People with HIV (agents) | 22.8k | 22.7k |
  | On ART | 9,125 | 9,128 |
  | OffART share of PLHIV | 42.0% | 42.0% |
  | ART coverage | 39.9% | 40.1% |
  | Diagnosed | 84.6% | 84.7% |

- **InterventionStatus** (`intervention_status.png`): all 12 states overlie each other.
- **Infant infections** (`infants.png`), new infections at ages 0–4 a year, 2005 / 2009 / 2012 / 2020: 180 / 235 / 165 / 189 compiled against 176 / 238 / 165 / 191 literal.

## Comparing runs with the same seed

The two variants are identical until the cascade's first event (the 1990 report matches exactly in every run). After that the literal cascade draws from its own random numbers and the compiled one from Starsim's distributions, so they can be compared only in distribution.

## What didn't compile cleanly, and why

1. **Order of updates.** EMOD updates a person's interventions breadth first along their chain, so a status set earlier can disqualify a later intervention; Starsim runs each delivery once per step for everyone, after all its possible triggers. Here it matters for ART staging (triggered from many chains) against the combination-PMTCT decision. The compiler makes deliveries that interact through InterventionStatus wait for each other along each chain using empty event terms (`event & ~event`; 55 terms on 14 deliveries). Without them, infant infections fell to 81 a year in 2009 (as if Option B weren't blocked); with a cruder global order they rose to 399 in 2005. This holds for this campaign but isn't guaranteed in general.
2. **Leaving a wait.** `delay=` checks eligibility only when the wait is due; EMOD drops a muxer as soon as the person is disqualified. The state holder ends those waits directly.
3. **Targets that vary by step** (NChooser) have no native form (`capacity=` is fixed), so they use `deliver=`.
4. **Smaller items:** forward references use `Cond('event', name, outcome)` from Starsim's internals; a delivery whose `start=` is after the sim's end raises; deliveries need `disease=` even when they touch no disease; deliveries reach people whose death is requested this step.

## Run time

Per run, one Numba thread per process, 30 at once: compiled 96 s (91–103), literal 60 s (57–63).

## Next steps

- In Starsim: an ordering between deliveries (`after=`), leaving a wait when ineligible, public references to deliveries by name, and capacity by step. These would remove most of the coded glue.
- Profile the compiled run.
- The condensed cascade (analysis 06), checked against this one.
