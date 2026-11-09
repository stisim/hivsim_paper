# 09: Transmissions by the source's cascade stage, EMOD against EMOD′ — summary

**Answer:** EMOD′ with the condensed cascade reproduces EMOD's split of transmissions by the source's cascade stage, within run-to-run noise. Over 1995–2024, 133 of 144 coarse cells (category × year) are within 2 standard errors. The exceptions are 2001–2003 (diagnosed, never on ART: EMOD′ 0.7–1.0 points higher, z 3.5–5.2) and a few single years, each under 0.7 points. In both models, transmission from people who interrupted ART rises from 7% in 2010 to 50% in 2020 and 54–55% in 2024. That is 2.5–3 times the top of MIHPSA Malawi's range across models (8.4–20.1%). The cause is EMOD's dropout trap (`docs/emod_inventory.md`, problem 1), which EMOD′ copies, so this endpoint is dominated by a known bug in the reference model.

Versions: emod-hiv 2.37.0 binary; Starsim v4 at ebd52c1a (Numba kernels, since the checkout's Rust kernels are stale, as in analyses 05–08); branch `emod-transmission`. Runs: parameter sets 0, 25, …, 225 × seeds 1–5 = 50 runs per model, to 2025.

## Shares of sexual transmissions (median, 5–95% over 50 runs each; `outputs/compare_coarse_years.csv`)

| Category | Year | EMOD | EMOD′ | EMOD′ − EMOD (mean ± SE, points) |
|---|---|---|---|---|
| Undiagnosed | 2010 | 51.6% (47.8–58.9) | 51.7% (48.0–57.8) | +0.2 ± 0.7 |
| | 2015 | 39.4% (35.1–42.9) | 39.1% (34.4–44.3) | −0.2 ± 0.6 |
| | 2020 | 32.4% (26.7–36.7) | 32.1% (28.4–36.1) | 0.0 ± 0.6 |
| | 2024 | 30.9% (25.2–34.7) | 29.9% (25.0–34.5) | −0.6 ± 0.6 |
| Diagnosed, never on ART | 2010 | 32.3% (13.4–36.8) | 31.9% (13.6–36.7) | +0.2 ± 1.5 |
| | 2015 | 9.6% (4.9–12.4) | 9.1% (4.9–12.5) | +0.1 ± 0.5 |
| | 2020 | 3.1% (2.2–4.3) | 3.5% (1.9–5.2) | +0.2 ± 0.2 |
| | 2024 | 2.2% (1.2–3.3) | 2.1% (1.3–3.2) | 0.0 ± 0.1 |
| On ART < 6 months | 2010 | 5.6% (3.8–8.4) | 5.8% (3.9–8.4) | 0.0 ± 0.3 |
| | 2015 | 6.2% (5.1–7.8) | 6.5% (5.0–7.7) | +0.1 ± 0.2 |
| | 2020 | 4.6% (3.6–5.3) | 4.6% (3.8–5.3) | +0.1 ± 0.1 |
| | 2024 | 4.0% (3.1–4.9) | 4.1% (3.5–5.1) | +0.2 ± 0.1 |
| On ART ≥ 6 months, suppressive arm | 2010 | 1.9% (1.3–3.7) | 1.8% (1.1–3.1) | −0.2 ± 0.1 |
| | 2015 | 3.8% (2.9–5.0) | 3.8% (2.8–4.9) | 0.0 ± 0.1 |
| | 2020 | 3.5% (2.7–4.4) | 3.3% (2.5–4.0) | −0.2 ± 0.1 |
| | 2024 | 3.0% (2.4–4.1) | 2.9% (2.0–4.0) | −0.2 ± 0.1 |
| On ART ≥ 6 months, non-suppressive arm | 2010 | 3.2% (2.0–4.7) | 3.0% (1.8–5.1) | −0.1 ± 0.2 |
| | 2015 | 6.8% (5.5–8.0) | 6.5% (5.3–7.9) | −0.2 ± 0.2 |
| | 2020 | 6.1% (5.2–7.4) | 6.0% (5.2–7.4) | 0.0 ± 0.1 |
| | 2024 | 5.5% (4.3–6.5) | 5.6% (4.6–6.5) | +0.2 ± 0.1 |
| Interrupted (off ART after starting) | 2010 | 7.4% (5.6–12.3) | 7.3% (5.5–12.0) | −0.1 ± 0.4 |
| | 2015 | 34.1% (31.6–39.8) | 34.7% (31.3–40.6) | +0.3 ± 0.5 |
| | 2020 | 50.2% (46.8–56.1) | 50.2% (47.3–55.6) | 0.0 ± 0.6 |
| | 2024 | 54.2% (49.3–60.8) | 55.4% (50.4–61.0) | +0.4 ± 0.7 |

Within "interrupted" (`compare_fine_years.csv`), almost all are first interruptions: 47.7% (EMOD) vs 47.3% (EMOD′) in 2020. Later interruptions add 2.9% (≥ 6 months since stopping) and 0.4% (< 6 months). Restarted and on ART < 6 months is 0.4%. The non-suppressive arm transmits more than the suppressive one despite being chosen for only 11–37% of starts, because it multiplies transmission by 0.65 against 0.04. Sexual transmissions a year (agents, median, EMOD vs EMOD′): 1299 vs 1308 (2010), 1434 vs 1455 (2020), 1560 vs 1660 (2024).

## Daniel's categories (`mwi_transmission.Rmd`)

Inputs:
- `TransmissionReport`: source ID, sex, age, `YEAR`.
- `ReportEventRecorder`: `HIVNewlyDiagnosed` only, from 2016, ages 15 and over, first per person.
- `ReportHIVART`: `StartingART` 1 = start, 0 = dropout, from 2016.

Transmissions from 2000 are counted by `floor(YEAR)` × source sex × source 5-year age bin 15–80; other ages are dropped by the join. Counts are scaled by 17,563,749 over the sim's 2018.5 population and averaged over runs. All comparisons are strict (event year < transmission year), and the windows are ≤ 0.5 years.

- **Undiagnosed:** all transmissions minus those after a diagnosis event.
- **Diagnosed, not on ART:** transmissions after diagnosis minus transmissions after the first ART start, as two separate counts.
- **On ART < 6 months after first initiation, CD4 at initiation < 200 or ≥ 200:** after the first start, before the first dropout (or with none), within 0.5 years.
- **On ART < 6 months after re-initiation:** the latest start k ≥ 2 before the transmission, with dropout k after the transmission or none, within 0.5 years; split by CD4 at that start.
- **On ART, VL suppressed or not:** currently on ART (the latest start, with its dropout not yet), split by fixed shares, which are not model output. `fun.onart.vl` assigns 8% suppressed and 92% not, while the under- and over-6-month versions assign 92% suppressed (a bug).
- **On ART < 6 months or > 6 months, by VL:** as above, split at 0.5 years since the latest start.
- **Dropped out (all):** the latest dropout k before the transmission, with start k+1 after it or none.
- **Dropped out, first interruption:** the same with k = 1.
- **Dropped out, 2nd+ interruption, < 6 or > 6 months:** k ≥ 2, split at 0.5 years since that dropout.
- A second section counts people, not transmissions, by `ARTStatus` from ReportHIVByAgeAndGender: diagnosed and ART-naive, on ART, OffART, and undiagnosed as the remainder.

Here these are made into one partition (`common/transmission.py`): `undiagnosed`, `diagnosed_naive`, `art_first_lt6m`, `art_reinit_lt6m`, `art_supp`, `art_nonsupp` (≥ 6 months, by arm), `off_first`, `off_later_lt6m`, `off_later_ge6m`. Coarse groups: undiagnosed, diagnosed_naive, art_lt6m, art_supp, art_nonsupp, interrupted.

Rules used:
- Starting ART counts as diagnosis.
- Events in the same monthly step count as before the transmission.
- 6 months means six monthly steps (182.5 days).
- The ART arm stands in for viral load, since neither model has one; EMOD's `SRC_VIRAL_LOAD` is a constant 10000.
- Daniel's CD4-at-initiation split, age bins and population scaling are not reproduced. Shares here are of all sexual transmissions, at all ages.

## Why Daniel's "diagnosed, not on ART" went negative, and the fix

It is the difference of two counts whose sets aren't nested: the "after first ART start" set includes sources with no diagnosis event before the transmission. With his settings this happens because the event recorder and the ART report both start in 2016, and diagnoses before age 15 are dropped, so someone diagnosed before 2016 who starts (or restarts) ART after 2016 is in the second count but not the first. Recomputed his way on our runs (`outputs/emod_checks.csv`, transmissions 2016–2025):
- **His settings:** 315–587 such transmissions per run, and the yearly difference is negative in all 50 runs (minimum −25 to −76).
- **Reports from 1960 at all ages:** 0–5 such transmissions per run (ART starts with no `HIVNewlyDiagnosed` event), and the difference is positive (8 to 59).

The fix has two parts: start the reports before anyone is diagnosed (here 1960), and classify each transmission once, as a set difference by event rather than by subtracting counts. Taking only the first ART start, as his comment suggests, doesn't fix it on its own: with the 2016 start, the "first" start in the report isn't the person's first.

## Checks

- **ARTStatus against the start and stop events:** they disagree for 0.34–0.56% of transmissions. Every disagreement has an ART start or stop in the same step, with EMOD's ordering going either way. None is on ART for 6 months or more, where the arm is read from ARTStatus.
- **EMOD event semantics** (set 0): `StartedART` and `StoppedART` alternate for every person, with no double starts; `ReportHIVART` has the same counts, so it adds nothing; `HIVNewlyDiagnosed` fires once per person; `StoppedART` means dropout only; about 6% of ART starts are restarts.
- The categories partition all transmissions by construction.

## Things found

- **Interrupted share 50–55% by 2020–2024** in EMOD as calibrated, against MIHPSA's 8.4–20.1%. This is the dropout trap: from 2009, 85–95% of dropouts are sent to a restart that never happens. This rerun of Daniel's model can't produce a value in MIHPSA's range. Fixing the trap (`dropout.p_restart` in the condensed cascade) would move this endpoint more than anything else.
- **EMOD′ is 0.7–1.0 points higher for "diagnosed, never on ART" in 2001–2003**, the years when ANC testing starts. Analysis 06 found differences in labels only for those years. Not investigated.
- **EMOD′'s suppressive-arm share is about 0.2 points low** in most years (z −1 to −2). Small, not investigated.

## Left to do

- HIVsim: `cascade_stage()` on its `CascadeHost` (diagnosed, on ART, and the arm), and the source of each infection from STIsim's infection log for `TransmissionLog`.
- Optima′: in agents, the same two hooks; in the ODE, attribution by each stage's share of the force of infection (`swap_design.md` 2.2), with `lost` split into before and after ART.
- `TransmissionLog` builds ART history by watching from the start of the run; a sim resumed from a 1985 snapshot needs that history carried in the snapshot.
- Only the condensed cascade is supported (`EmodHost.cascade_stage` reads the arm from the cascade's `people.arm`).
- Ask Daniel which EMOD value went into MIHPSA's interrupted range.
- Move to Starsim 4.0.0rc3 with the rest of the repo.
