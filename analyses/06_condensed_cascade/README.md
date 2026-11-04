# 06: The condensed cascade against the compiled campaign

**Question:** how much does EMOD′ change when Daniel's campaign (about 60 events, compiled into 95 deliveries in analysis 05) is replaced by the condensed cascade in `common/cascade/` (one table of parameters by pathway, a `CascadeHost` interface to the disease), with its parameters read from the same campaign? The condensed cascade is the one the component swaps will use, so the differences here are the cost of condensing, measure by measure.

**What's run:**
- Starsim v4 from `~/starsim` (commit ebd52c1a, 4.0.0rc2, with PR #1538) on `PYTHONPATH`, for both variants.
- `harness.make_sim(pid, seed, cascade=True, compiled=True)` (compiled) and `harness.make_sim(pid, seed, cascade=True, condensed=True)` (condensed, parameters from `models/emod/v4/extract_cascade.py` for each set), for the 10 parameter sets of analyses 04 and 05, 3 seeds each (seeded from their run numbers): 60 runs.
- Measured as in analyses 04 and 05: the common targets (`emod_targets.measure()`), the cascade over time, places in care of people with HIV, and infections at ages 0-4; and run time. The compiled run's InterventionStatus is grouped into the condensed cascade's places in care (`run.py`, `GROUPS`).
- `run.py` runs and writes `outputs/`; `run.py --plot` draws `figures/`.

**What the condensed cascade leaves out or changes, so where differences are expected:**
- One eligibility rule, staging's, for the pre-ART recheck too (EMOD's two tables differ a little).
- Traditional circumcision stays in EmodHIV (it is a property of the population, not the programme).
- No InterventionStatus properties: a person has one place in care, and testing pathways reach those not in care (none, uptake, testing). EMOD's per-intervention disqualifying lists are condensed to that.
- Waits (uptake, testing loop, pre-ART, dropout) belong to a place in care: one ends without effect if the person has moved on, rather than being removed when they move.
- Ordering interactions are parameters: PMTCT's block for women sent to staging (`pmtct.blocked_if_staged`), the dropout trap (`dropout.p_restart` 0).

**How runs compare:** as in analysis 05, in distribution: medians and ranges over the 30 runs of each.
