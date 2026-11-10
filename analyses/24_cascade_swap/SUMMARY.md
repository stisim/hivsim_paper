# 24: The cascade step (C) between EMOD′-fixed and Optima′, both ways — summary

**Rerun, 10 Nov (Robyn's decision 7: Optima's PMTCT on whenever Optima's cascade runs on agents):** direction 1 below is with Optima's PMTCT drawn per person, now the default of Optima's C slot. The numbers with PMTCT skipped are kept only as labelled superseded lines. Direction 2 and the open item don't run Optima's cascade on another model's agents and are unchanged.

**Answer:** in both directions, swapping the cascade moves the Figure 3A shares a long way towards the other model.
- It moves diagnosed never treated and on ART by about 0.5–0.75 of the gap between our two ends, and interrupted by the whole gap or more.
- On those shares the effect is roughly symmetric. Undiagnosed is not (0.7–1.1 of the gap into EMOD′, about 2.0 into Optima′).
- New infections don't follow the gap. EMOD′-fixed has more infections than Optima′ (39.7k against 33.6k in 2020), yet EMOD's cascade lowers infections in either host: −13 to −19% in Optima′, and Optima's cascade raises EMOD′'s by +43% as it stands (main arm), +11% (+6% sexual) once the other host artefacts are removed (1sp). So the cascade explains the direction of the share split but not the total: that comes from the other components.

**Ends** are ours, measured the same way at both: EMOD′-fixed with its own cascade (48 cells), and Optima′ ODE (the 2022 project). Published 2020: EMOD 38.5 / 3.7 / 40.6 / 17.1, Optima 45.2 / 28.2 / 17.6 / 9.0. Our Optima (analysis 19) is 42.3 / 40.7 / 11.4 / 5.7; Optima′ gives 42.1 / 39.3 / 12.7 / 5.8. Shares are undiagnosed / diagnosed never treated / on ART / interrupted, % of sexual transmissions.

**Runs** (Starsim /tmp/ss-nov09 f1273734, STIsim /tmp/sti-613 ebbc8e7, env hp_rc3; `run.py`, `tables.py`; all numbers in `outputs/report.txt`, `deltas.csv`, `gap.csv`, `open_item.csv`):
- **Direction 1:** EMOD′-fixed own cascade against Optima's table on EMOD′'s people (`models/optima/v4/components.py`, with Optima's PMTCT per person), EMOD sets 0/25/50 × 16 seeds = 48 cells, paired through the runner. Sensitivity arms (`variants/`): 1s, the ART arm switches on treatment failure and regain; 1sp, as 1s plus EMOD′'s own VMMC and PrEP. No errors. Rerun 10 Nov with PMTCT on (`logs/emod_pmtct_default.log`, 192 runs): `own` is identical to the stored own, and main and 1sp are identical to the addendum's `optima_pm` and `optima_swp_pm` (96 of 96 cell-years each).
- **Direction 2:** Optima′ own against EMOD′-fixed's condensed table (`p_restart` = 1, checked), ODE for sets 0/25/50 (25 and 50 have identical tables, so two distinct EMOD tables); agents at 100k with `Outputs(OptimaHost)`, 16 seeds per set against 64 of the own cascade (unpaired).
- **Analysis 23's open item:** Optima′ agents, Optima's table against own cascade, 64 seeds each, same code and environment.

## Direction 1 (Optima's cascade into EMOD′-fixed), 2020, paired Δ (95% CI), Optima's PMTCT on

| 2020 | EMOD′-fixed | + Optima's cascade | Δ main | Δ 1s | Δ 1sp |
|---|---|---|---|---|---|
| New infections | 39.7k | 56.6k | +16.9k (14.7, 19.1) | +7.3k (4.8, 9.9) | +4.3k (1.8, 6.8) |
| Sexual transmissions | 37.8k | 52.5k | +14.8k (12.7, 16.8) | +5.5k (3.2, 7.9) | +2.4k (0.1, 4.7) |
| Mother-to-child | 1.9k | 4.1k | +2.1k (1.9, 2.4) | +1.8k (1.5, 2.1) | +1.9k (1.6, 2.2) |
| Shares, % (undx / dx naive / on ART / interrupted) | 35.3 / 4.9 / 45.9 / 13.9 | 39.9 / 22.6 / 34.3 / 3.2 | +4.7 (3.3, 6.1) / +17.7 (16.3, 19.1) / −11.7 (−13.5, −9.9) / −10.7 (−11.3, −10.0) | +7.8 / +22.9 / −21.1 / −9.6 | +6.3 (4.8, 7.8) / +23.1 (21.2, 25.0) / −19.7 (−21.4, −18.1) / −9.7 (−10.3, −9.0) |
| PLHIV | 937.7k | 1151.8k | +214.1k (195.0, 233.2) | +209.9k | +194.2k (169.5, 219.0) |
| *PMTCT skipped (superseded)* | | *60.6k new infections* | *+20.9k (17.9, 23.9); sexual +15.5k; MTCT +5.4k; shares +3.7 / +19.3 / −12.2 / −10.8; PLHIV +263.4k* | *+11.8k; sexual +6.5k; MTCT +5.3k; +7.4 / +24.6 / −22.4 / −9.6* | *+8.6k (5.8, 11.4); sexual +3.5k (1.2, 5.7); MTCT +5.1k; +6.3 / +24.8 / −21.2 / −9.9; PLHIV +242.9k* |

2024, main arm: +20.8k (18.3, 23.3) new infections; shares +7.4 / +18.3 / −13.0 / −12.6; 1sp +3.1k (0.1, 6.2) (superseded, PMTCT skipped: main +27.2k, 1sp +8.6k). PMTCT barely moves the shares (Figure 3A counts sexual transmission); in 2020 it takes the extra mother-to-child infections from +5.4k to +2.1k, and sexual transmissions from +15.5k to +14.8k (main arm).

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
| Undiagnosed | 0.68 | 1.13 | 0.92 | 2.04 | 1.99 |
| Diagnosed never treated | 0.51 | 0.67 | 0.67 | 0.52 | 0.72 |
| On ART | 0.35 | 0.64 | 0.59 | 0.57 | 0.57 |
| Interrupted | 1.32 | 1.19 | 1.20 | 1.60 | 2.28 |
| New infections | −2.8 | −1.2 | −0.7 | −0.84 | −0.92 |

Direction 1 with Optima's PMTCT on. Superseded (PMTCT skipped), main / 1s / 1sp: undiagnosed 0.54 / 1.07 / 0.92, diagnosed never treated 0.56 / 0.72 / 0.72, on ART 0.37 / 0.67 / 0.64, interrupted 1.34 / 1.19 / 1.23, new infections −3.4 / −1.9 / −1.4. 2024 looks the same.

**Symmetry:**
- **Diagnosed never treated and on ART: roughly symmetric.** With arm switching, the cascade carries about 0.6–0.7 of the gap whichever model it sits in.
- **Interrupted:** the swap closes the whole gap or overshoots it both ways; EMOD's cascade overshoots more inside Optima′ (1.6–2.3 against 1.2–1.3).
- **Undiagnosed: not symmetric** (about 0.7–1.1 into EMOD′, about 2.0 into Optima′).
- **New infections: same sign, different size, and opposite to the gap.** EMOD's cascade gives fewer infections in either host (sexual transmissions +6% for EMOD′ with Optima's cascade (1sp), −9 to −15% for Optima′ with EMOD's), yet EMOD′-fixed has more infections than Optima′.

## Analysis 23's ~1% gap: not real

At 64 seeds each, paired by seed, the 2000–2020 mean of table/own − 1: PLHIV −0.40% (−1.08, +0.27), new infections +0.01% (−1.08, +1.11), diagnosed, on ART and suppressed within 0.1%; |z| > 2 in 0–4 of 21 years per output. 2020 shares differ by at most 1 point, every CI covering 0. Analysis 23 compared against analysis 20's separately run seeds; the gap was noise.

## Caveats

1. **Optima's cascade on EMOD′ is not quite Optima's.** Without arm switching (analysis 23 provisional decision 5), only 69% of people on ART are on the suppressive arm against 93% in Optima′: the biggest artefact (arm switching cuts the rise in infections from +16.9k to +7.3k; +20.9k to +11.8k with PMTCT skipped). Optima's PMTCT is now drawn per person (decision 7, decided); skipping it had added +5.4k mother-to-child infections in 2020 (+2.1k remain with it, real cascade content: Optima's PMTCT numbers start small, and EMOD′'s early ANC uptake is high, see the addendum). EMOD's VMMC and PrEP belong to EMOD's cascade and leave with it; putting them back takes off a further 3.0k in 2020 and 6.1k in 2024. The number on ART is Optima's Malawi target in people, applied to EMOD′'s population.
2. **PLHIV rises by +214k with Optima's cascade on EMOD′** (+108k more infections, −128k fewer HIV deaths over 1990–2020; +263k with PMTCT skipped, superseded). Diagnosed in the addendum below: the part PMTCT doesn't remove is Optima's larger ART numbers in 2005–2010 (deaths averted), a real cascade difference.
3. The new-infection gap is small (6.1k), so fractions of it are unstable; use the Δs.
4. Our Optima end is the 2022 project. Against the published Optima, the gap in diagnosed never treated would be 24.5 points instead of 34.4.
5. EMOD′-fixed is uncalibrated; only two distinct EMOD cascade tables were tested.
6. EMOD′ counts are agents × 100, not rescaled to the census.

## Code (as first run: no default behaviour changes; tests 24 passed, `logs/tests.log`. Since 10 Nov, Optima's C slot draws Optima's PMTCT by default: see Next)

- `models/optima/v4/components.py` (new): Optima's cascade as a runner slot, on the host of the model supplying V, scaled by the D model's `SCALE`.
- `common/cascade/pathways.py`: `TableCascade(scale=None)` (EMOD′ runs with `pop_scale` 1, so Optima's ART target in people would otherwise be read as agents).
- `models/emod/v4/cascade_host.py`: `EmodHost` keeps `people.arm` itself, which the table cascade needs.
- Name clash: `models/emod/v4` and `models/optima/v4` both have an `extract_cascade.py`, and the runner puts both folders on `sys.path`; a hybrid with Optima's V and EMOD's cascade would import the wrong one (documented; `run.py` loads EMOD's by path).

## Next

- Robyn's decision: arm switching on `EmodHost` as the default for table cascades (analysis 23 decision 5).
- Decided (Robyn, 10 Nov): Optima's PMTCT on whenever Optima's cascade runs on agents (analysis 23 decision 7). Optima's C slot passes `TableCascade(pmtct=True)` (`models/optima/v4/components.py`, `PMTCT = True`; False reproduces the superseded runs); `TableCascade`'s own default stays False, since analysis 23's HIVsim table runs and the Optima′ agent runs use it directly (on Optima′'s own agents it adds nothing: no one is pregnant, and Optima′ applies its PMTCT inside births). Test: `common/swaps/tests/test_runner.py::test_optima_c_on_emod_draws_pmtct` (fails with `PMTCT = False`).
- Run the V, RN and D steps to find where the total-infection gap comes from.

## Addendum: why PLHIV and infections are higher early (caveat 2)

**Answer:** the early excess is mother-to-child, from 2001, because EMOD's PMTCT leaves with its cascade and Optima's PMTCT is skipped on agent hosts (an artefact of decision 7). Drawing Optima's PMTCT number per person removes two thirds of the extra mother-to-child infections. The rest of the PLHIV rise is a real cascade effect: Optima's ART numbers are Malawi's programme data, larger than EMOD′-fixed's mechanistic ART in 2005–2012, so fewer people die.

**Runs** (`early.py`, `run.py emod_pm`; `outputs/early_report.txt`, `early_by_year.csv`, `early_cumulative.csv`, `early_cells.csv`, `early_track.csv`):
- Direction 1 again, sets 0/25/50 × 16 seeds, paired: own; `optima` as before; `optima_pm`, the same with Optima's PMTCT per person; `optima_swp_pm`, as 1sp (arm switching, EMOD's VMMC and PrEP) with Optima's PMTCT.
- New code: `TableCascade(pmtct=True)`, default False. At 12 weeks of pregnancy, a diagnosed woman not on ART starts PMTCT (efficacy 0.9, to the end of pregnancy) with Optima's share. That share is the table's number per year over the diagnosed women giving birth per year, counted on ART or not as Optima counts them, at most 1.
- Rerunning `own` and `optima` in the new code gives tables byte-identical to the stored ones (96 of 96 cells), so the defaults are unchanged.
- Tracked runs (set 0, seeds 1–8, to 2012, `early.py track`): VMMC, STI flag, PrEP, pregnant women with HIV (diagnosed, on ART, on PMTCT), children with HIV and on ART.

**First divergence: 2001, mother-to-child.**
- 1990–2000: Optima's table diagnoses far more people (294k against 11k in 2000) and has a few thousand on ART (EMOD′ has none before 2004). This barely changes infections (−3.7k, CI covers 0) and saves 11k deaths.
- 2001: EMOD's ANC testing and sdNVP start. Under EMOD's cascade about 28k women at any time are pregnant, have HIV and are on PMTCT by 2001, and about 35k in 2003. Under Optima's table, none are.
- 2001–2004: +41.0k (39.2, 42.8) mother-to-child infections; sexual transmissions −0.4k (−3.2, +2.3). Children with HIV: 95k against 78k in 2004.
- Sexual transmissions are lower in 2005–2010 (−27.7k), because more people are on ART. They are only higher from 2011 (+34k in 2011–15, +69k in 2016–20), from the larger PLHIV and from the arm not switching (caveat 1).

**What else EMOD's campaign carries, and when:**
- Testing (1990), symptomatic presentation (1990), ANC testing and PMTCT (2001), infant testing (2004), VMMC (2011), PrEP (2019).
- Traditional circumcision, the STI flag, condoms and seeding are in V and RN and stay. Circumcision among men 15–49 and the STI flag are the same in both arms, and PrEP is 0 before 2012.
- So before 2011, only testing, ART and PMTCT leave with the cascade.

**Cumulative 1990–2020, paired Δ (95% CI), people:**

| | New infections | Sexual | Mother-to-child | HIV deaths | PLHIV 2020 |
|---|---|---|---|---|---|
| `optima` (as before) | +234.6k (201.0, 268.3) | +71.9k | +162.7k (151.1, 174.3) | −55.8k | +263.4k (241.8, 284.9) |
| `optima_pm` | +107.9k (75.1, 140.6) | +56.0k | +51.9k (46.1, 57.6) | −127.6k | +214.1k (195.0, 233.2) |
| `optima_swp_pm` | +73.5k (32.3, 114.6) | +20.9k | +52.6k (46.0, 59.2) | −143.1k | +194.2k (169.5, 219.0) |

With PMTCT, HIV deaths fall further, because the extra infected children of `optima` died. The remaining PLHIV rise comes from deaths averted in 2005–2010 (−138k, against −110k without PMTCT): Optima's ART target is 29k in 2005 and 151k in 2008, against 15.5k and 92k in EMOD′-fixed's own cascade (set 0). That is a real cascade difference, not a scaling one: EMOD′'s population (13.1M in 2008, agents × 100) matches Malawi's census.

The extra mother-to-child infections that remain (+37k in 2001–2004) are also real cascade content, but they flatter EMOD′:
- Optima's PMTCT numbers start small: 1.4k in 2002, 7.3k in 2004.
- EMOD′'s ANC uptake in 2001–2003 is high (about 73% in 2001) because EMOD's sigmoid class ignores the campaign's Time_Value_Map, which gives 31% in 2001 (`inputs.FIXES['anc']`, off in EMOD′-fixed).

**Corrected direction-1 headline, 2020, paired Δ (95% CI):**

| 2020 | main (`optima`) | `optima_pm` | 1sp (stored) | `optima_swp_pm` |
|---|---|---|---|---|
| New infections | +20.9k (17.9, 23.9) | +16.9k (14.7, 19.1) | +8.6k (5.8, 11.4) | +4.3k (1.8, 6.8) |
| Sexual | +15.5k | +14.8k (12.7, 16.8) | +3.5k (1.2, 5.7) | +2.4k (0.1, 4.7) |
| Mother-to-child | +5.4k | +2.1k | +5.1k | +1.9k |
| Shares (undx / dx naive / on ART / interrupted) | +3.7 / +19.3 / −12.2 / −10.8 | +4.7 / +17.7 / −11.7 / −10.7 | +6.3 / +24.8 / −21.2 / −9.9 | +6.3 / +23.1 / −19.7 / −9.7 |
| PLHIV | +263.4k | +214.1k | | +194.2k |

- 2024: `optima_pm` +20.8k new infections (18.3, 23.3); `optima_swp_pm` +3.1k (0.1, 6.2).
- With every host artefact removed, Optima's cascade on EMOD′ raises 2020 new infections by about 11% (+4.3k), and sexual transmissions by about 6%.
- The shares barely move: PMTCT changes mother-to-child infections, which the Figure 3A shares don't count. The gap fractions in the table above stand.
- Direction 1's conclusion is unchanged: the cascade carries the share split, not the total.

The table above is now the headline (direction 1 at the top, rerun 10 Nov with PMTCT on by default; identical numbers).

**Slot-boundary questions for Robyn:**
1. Is PMTCT in C? Both models treat it as cascade, so here it goes with C. Decided 10 Nov: yes, on by default for Optima's C slot.
2. Decided 10 Nov: PMTCT on for Optima's cascade on agent hosts, with Optima's denominator (diagnosed women giving birth, on ART or not).
3. EMOD's VMMC and PrEP: C (they leave, as in the main arm) or prevention outside C (kept, as in 1sp)? They don't matter before 2011.
4. EMOD′'s high early ANC uptake: should the `anc` fix go into EMOD′-fixed? That is an EMOD′ question, not a swap one.

Tests: `tests/` and `common/swaps/tests` (`logs/tests_addendum.log`).
