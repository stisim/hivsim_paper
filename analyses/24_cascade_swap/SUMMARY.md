# 24: The cascade step (C) between EMOD′-fixed and Optima′, both ways — summary

**Answer:** in both directions, swapping the cascade moves the Figure 3A shares a long way towards the other model.
- It moves diagnosed never treated and on ART by about 0.5–0.75 of the gap between our two ends, and interrupted by the whole gap or more.
- On those shares the effect is roughly symmetric. Undiagnosed is not (0.5–1.1 of the gap into EMOD′, about 2.0 into Optima′).
- New infections don't follow the gap. EMOD′-fixed has more infections than Optima′ (39.7k against 33.6k in 2020), yet EMOD's cascade lowers infections in either host: −13 to −19% in Optima′, and Optima's cascade raises EMOD′'s by +53% as it stands, +9% sexual once host artefacts are removed. So the cascade explains the direction of the share split but not the total: that comes from the other components.

**Ends** are ours, measured the same way at both: EMOD′-fixed with its own cascade (48 cells), and Optima′ ODE (the 2022 project). Published 2020: EMOD 38.5 / 3.7 / 40.6 / 17.1, Optima 45.2 / 28.2 / 17.6 / 9.0. Our Optima (analysis 19) is 42.3 / 40.7 / 11.4 / 5.7; Optima′ gives 42.1 / 39.3 / 12.7 / 5.8. Shares are undiagnosed / diagnosed never treated / on ART / interrupted, % of sexual transmissions.

**Runs** (Starsim /tmp/ss-nov09 f1273734, STIsim /tmp/sti-613 ebbc8e7, env hp_rc3; `run.py`, `tables.py`; all numbers in `outputs/report.txt`, `deltas.csv`, `gap.csv`, `open_item.csv`):
- **Direction 1:** EMOD′-fixed own cascade against Optima's table on EMOD′'s people (`models/optima/v4/components.py`), EMOD sets 0/25/50 × 16 seeds = 48 cells, paired through the runner. Sensitivity arms (`variants/`): 1s, the ART arm switches on treatment failure and regain; 1sp, as 1s plus EMOD′'s own VMMC and PrEP. No errors; a rerun of the main arms was byte-identical.
- **Direction 2:** Optima′ own against EMOD′-fixed's condensed table (`p_restart` = 1, checked), ODE for sets 0/25/50 (25 and 50 have identical tables, so two distinct EMOD tables); agents at 100k with `Outputs(OptimaHost)`, 16 seeds per set against 64 of the own cascade (unpaired).
- **Analysis 23's open item:** Optima′ agents, Optima's table against own cascade, 64 seeds each, same code and environment.

## Direction 1 (Optima's cascade into EMOD′-fixed), 2020, paired Δ (95% CI)

| 2020 | EMOD′-fixed | + Optima's cascade | Δ main | Δ 1s | Δ 1sp |
|---|---|---|---|---|---|
| New infections | 39.7k | 60.6k | +20.9k (17.9, 23.9) | +11.8k (8.9, 14.6) | +8.6k (5.8, 11.4) |
| Sexual transmissions | 37.8k | 53.2k | +15.5k | +6.5k | +3.5k (1.2, 5.7) |
| Mother-to-child | 1.9k | 7.3k | +5.4k | +5.3k | +5.1k |
| Shares, % | 35.3 / 4.9 / 45.9 / 13.9 | 39.0 / 24.2 / 33.8 / 3.1 | +3.7 / +19.3 / −12.2 / −10.8 (each CI about ±1–2) | +7.4 / +24.6 / −22.4 / −9.6 | +6.3 / +24.8 / −21.2 / −9.9 |

2024, main arm: +27.2k (23.7, 30.7) new infections; shares +7.8 / +20.2 / −15.3 / −12.8; 1sp +8.6k.

## Direction 2 (EMOD′-fixed's cascade into Optima′), 2020

| 2020 | Optima′ (ODE) | + EMOD's cascade, ODE (set 0 / sets 25 and 50) | Δ agents (95% CI) |
|---|---|---|---|
| New infections | 33.6k | −6.4k / −4.4k | −4.6k (−5.6, −3.6) |
| Sexual transmissions | 28.7k | −4.3k / −2.6k | −4.0k |
| Shares, % | 42.1 / 39.3 / 12.7 / 5.8 | −13.7 / −20.9 / +20.7 / +13.9 (set 0); −14.1 / −16.3 / +17.9 / +12.5 (sets 25, 50) | −13.4 / −25.8 / +19.2 / +19.9 |

2024: ODE −3.8k (set 0), agents −2.8k (−3.9, −1.7). The agents count more as interrupted than the ODE, because their ever-treated tracer also counts people sent back to uptake or care after ART. Reference ODE runs in 2020: without VMMC/PrEP flows 27.9k new infections (against 27.1k); copied-bug EMOD table 88.0k new infections and 65.7% interrupted.

## How much of the gap the cascade carries (2020)

Gap = Optima′ − EMOD′-fixed: undiagnosed +6.9, diagnosed never treated +34.4, on ART −33.2, interrupted −8.1 points; new infections −6.1k. Fraction = the swap's Δ over that gap (direction 2 with its sign reversed); 1 means the cascade alone closes it.

| Quantity | Dir. 1 main | Dir. 1s | Dir. 1sp | Dir. 2 ODE | Dir. 2 agents |
|---|---|---|---|---|---|
| Undiagnosed | 0.54 | 1.07 | 0.92 | 2.04 | 1.99 |
| Diagnosed never treated | 0.56 | 0.72 | 0.72 | 0.52 | 0.72 |
| On ART | 0.37 | 0.67 | 0.64 | 0.57 | 0.57 |
| Interrupted | 1.34 | 1.19 | 1.23 | 1.60 | 2.28 |
| New infections | −3.4 | −1.9 | −1.4 | −0.84 | −0.92 |

2024 looks the same.

**Symmetry:**
- **Diagnosed never treated and on ART: roughly symmetric.** With arm switching, the cascade carries about 0.6–0.7 of the gap whichever model it sits in.
- **Interrupted:** the swap closes the whole gap or overshoots it both ways; EMOD's cascade overshoots more inside Optima′ (1.6–2.3 against 1.2–1.3).
- **Undiagnosed: not symmetric** (about 0.5–1.1 into EMOD′, about 2.0 into Optima′).
- **New infections: same sign, different size, and opposite to the gap.** EMOD's cascade gives fewer infections in either host (sexual transmissions +9% for EMOD′ with Optima's cascade (1sp), −9 to −15% for Optima′ with EMOD's), yet EMOD′-fixed has more infections than Optima′.

## Analysis 23's ~1% gap: not real

At 64 seeds each, paired by seed, the 2000–2020 mean of table/own − 1: PLHIV −0.40% (−1.08, +0.27), new infections +0.01% (−1.08, +1.11), diagnosed, on ART and suppressed within 0.1%; |z| > 2 in 0–4 of 21 years per output. 2020 shares differ by at most 1 point, every CI covering 0. Analysis 23 compared against analysis 20's separately run seeds; the gap was noise.

## Caveats

1. **Optima's cascade on EMOD′ is not quite Optima's.** Without arm switching (analysis 23 provisional decision 5), only 69% of people on ART are on the suppressive arm against 93% in Optima′: the biggest artefact (arm switching cuts the rise in infections from +20.9k to +11.8k). PMTCT is skipped on other hosts (decision 7), adding +5.4k mother-to-child infections (not in the Figure 3A shares, which count sexual transmission). EMOD's VMMC and PrEP belong to EMOD's cascade and leave with it; putting them back takes off a further 3.2k in 2020 and 7.2k in 2024. The number on ART is Optima's Malawi target in people, applied to EMOD′'s population.
2. **PLHIV rises by +263k with Optima's cascade on EMOD′** (+235k more infections, −56k fewer HIV deaths over 1990–2020); infections are already higher by 2004. Not yet investigated.
3. The new-infection gap is small (6.1k), so fractions of it are unstable; use the Δs.
4. Our Optima end is the 2022 project. Against the published Optima, the gap in diagnosed never treated would be 24.5 points instead of 34.4.
5. EMOD′-fixed is uncalibrated; only two distinct EMOD cascade tables were tested.
6. EMOD′ counts are agents × 100, not rescaled to the census.

## Code (no default behaviour changes; tests 24 passed, `logs/tests.log`)

- `models/optima/v4/components.py` (new): Optima's cascade as a runner slot, on the host of the model supplying V, scaled by the D model's `SCALE`.
- `common/cascade/pathways.py`: `TableCascade(scale=None)` (EMOD′ runs with `pop_scale` 1, so Optima's ART target in people would otherwise be read as agents).
- `models/emod/v4/cascade_host.py`: `EmodHost` keeps `people.arm` itself, which the table cascade needs.
- Name clash: `models/emod/v4` and `models/optima/v4` both have an `extract_cascade.py`, and the runner puts both folders on `sys.path`; a hybrid with Optima's V and EMOD's cascade would import the wrong one (documented; `run.py` loads EMOD's by path).

## Next

- Robyn's decisions: arm switching on `EmodHost` as the default for table cascades (analysis 23 decision 5); Optima's PMTCT drawn per person on agent hosts (decision 7).
- Diagnose why infections are already higher in the early 2000s with Optima's cascade on EMOD′.
- Run the V, RN and D steps to find where the total-infection gap comes from.
