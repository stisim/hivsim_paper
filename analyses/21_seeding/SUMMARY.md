# 21: Seeding options on the path HIVsim → EMOD′-fixed — summary

Measurements for decision D11 (Robyn's; provisional: "condition on take-off"). Nothing is decided here.

**Answer:** the Figure 3A share Δs don't depend on the seeding rule. Under every option, each step's share Δ stays within 0.03 of analysis 16's addendum, and no sign changes. The count Δs do depend on it, and conditioning on take-off doesn't fix that:
- On EMOD′'s network, more seeds give a larger epidemic by 2020, not only a higher chance of one. HEEE's runs that take off average 27.8k new infections in 2020 under EMOD's rule, against 43–47k when take-off is near certain.
- So the addendum's count Δs at RN and D, conditional or not, are partly seeding. With near-certain take-off: D's Δ new infections goes from +17.5k (+12.3k conditional) to between −6k and 0; RN's from −27k (−40k conditional) to between −20k and −23k; V's from +20k (+40k conditional) to +37k.
- Option (a), EMOD's rule scaled per capita, changes almost nothing (2.0 expected seeds against 1.7).

Tier 1, HIVsim provisional, one EMOD′ parameter set (0). Take-off means prevalence ≥ 0.5% at the end of 2020. 2400 runs (12 hybrids × 200 cells, commit 098eb29), 0 errors. The default is bit-identical: option 0's 1000 stored tables match the addendum's (`outputs/bitcheck.csv`). Tests pass (19).

**Code (all off by default):** `EmodHIV(seeding=None)` takes a rule `seeding(hiv) → uids` (None = EMOD's, same code path); `adapters.seed_per_capita(rate)` (option a; reuses EMOD's per-HIGH draws, so its seeds nest with option 0's in each cell; option c likewise); model files `components_fixed_scaled.py` (a) and `components_fixed_seed5.py` (c, `seed_coverage` × 5); option (b) is an explicit S slot: with `S=hivsim`, EMOD's `make_H` uses `hivsim.components.seeding_for()` (HIVsim's `init_prev_data` rule; on EMOD′'s network a `RISK_PREV` table of HIVsim's seeding probability by EMOD's Risk class and sex, measured on HIVsim's 1985 population). Combinations not built raise. See `common/swaps/README.md`, "Seeding".

**Seeds per run** (`outputs/rates.json`): EMOD′-fixed set 0 seeds 1.09% of Risk HIGH (not 2%) in 1980.92: 12.3 seeds (6–20) at 61,000 agents, 1.97e-4 per person. HIGH is 1.81% of everyone alive in EMOD′ and 1.48% in HIVsim's 1985 population (153 of 10,300).

| Rule | Expected seeds at 10,000 agents |
|---|---|
| EMOD's (option 0) | 1.7 |
| (a) scaled per capita | 2.0 |
| (c) 5× | 8.3 |
| (b) HIVsim's `init_prev_data` | 30.8 |

**Take-off** (95% Wilson intervals):

| Hybrid | 0 (EMOD rule) | (a) scaled | (b) S = HIVsim | (c) 5× |
|---|---|---|---|---|
| HHEE | 144/200, 0.72 (0.65–0.78) | 153/200, 0.77 (0.70–0.82) | 200/200, 1.00 (0.98–1) | 200/200, 1.00 (0.98–1) |
| HEEE | 164/200, 0.82 (0.76–0.87) | 155/200, 0.78 (0.71–0.83) | 200/200, 1.00 (0.98–1) | 200/200, 1.00 (0.98–1) |

HHHH, HHHE, EEEE and EEEE-c take off 200/200. Paired Δp per step: (a) V −0.24, RN +0.01 (ns), D +0.23; (b) and (c) 0 at every step on paths that don't pass through an option-0 hybrid; (b) with S after V, the S step is −0.28 (HHEE-S → HHEE); with S after RN, −0.18 (HEEE-S → HEEE). Taken-off HEEE runs reach 3.3% prevalence on average (as low as 0.5%), against 5.4–6.1% with near-certain take-off: HEEE's epidemic starts late from few seeds.

## Headline 2020, paired Δ per step, all cells where defined

Shares: undiagnosed / diagnosed never treated / on ART / interrupted (`outputs/headline_2020.csv`, `headline_2020_wide.csv`).

| Option, path | Output | C | V | S | RN | D |
|---|---|---|---|---|---|---|
| 0 (addendum) | New infections | +4.5k | +20.2k (14.9, 25.5) | | −27.2k (−31.9, −22.4) | +17.5k (15.1, 19.8) |
| | Shares | −.286 / +.047 / +.097 / +.142 | +.220 / −.001 / −.110 / −.109 | | −.281 / −.005 / +.211 / +.074 | −.023 / −.001 / +.013 / +.012 |
| (a) scaled | New infections | same | +22.1k (17.1, 27.1) | | −30.9k (−35.6, −26.1) | +19.3k (16.9, 21.7) |
| | Shares | same | +.210 / +.001 / −.102 / −.109 | | −.272 / −.011 / +.199 / +.084 | −.021 / +.002 / +.013 / +.005 |
| (b) C, V, S, RN, D | New infections | same | +37.2k (34.3, 40.2) | −17.0k (−22.8, −11.2) | as 0 | as 0 |
| | Shares | same | +.219 / −.014 / −.105 / −.100 | −.002 / +.010 / −.001 / −.007 | as 0 | as 0 |
| (b) C, V, RN, S, D | New infections | same | +37.2k | −23.7k (−26.3, −21.1) | −20.5k (−23.6, −17.3) | +17.5k |
| | Shares | same | as above | +.027 / .000 / −.014 / −.014 | −.301 / +.007 / +.212 / +.082 | as 0 |
| (b) S taken with D (C, V, RN, D+S) | New infections | same | +37.2k | | −20.5k | D+S −6.2k (−7.8, −4.7) |
| | Shares | same | as above | | as above | +.004 / −.001 / .000 / −.003 |
| (c) 5× | New infections | same | +36.7k (33.9, 39.6) | | −23.4k (−26.3, −20.4) | −2.9k (−4.6, −1.2) |
| | Shares | same | +.224 / −.014 / −.107 / −.103 | | −.297 / +.007 / +.210 / +.081 | −.005 / −.001 / +.005 / +.002 |
| (c) 5×, EEEE 5× too | New infections | same | as (c) | | as (c) | +0.1k (−1.4, +1.7) |

Under (b) and (c), every cell takes off at every step, so "all cells", "both arms took off" and "every hybrid of the path took off" are the same 200 cells.

## Alternative summaries of the existing results (option 0)

- **Take-off as its own output:** p = 1, 1, 0.72, 0.82, 1 along the path; Δp C 0, V −0.28, RN +0.10, D +0.18.
- **Both arms took off:** analysis 16's addendum (b): new infections V +39.6k, RN −40.5k, D +12.3k.
- **Every hybrid of the path took off** (133 cells, one common set, so the steps add up to the endpoint difference on those cells):

| Output | C | V | RN | D |
|---|---|---|---|---|
| New infections | +4.4k | +39.7k | −40.5k | +11.1k |
| Undiagnosed | −0.286 | +0.222 | −0.288 | −0.015 |
| On ART | +0.096 | −0.111 | +0.218 | +0.005 |
| Interrupted | +0.144 | −0.109 | +0.074 | +0.013 |

The shares here are within 0.003 of "both arms" and "all cells".

## Robust to the seeding choice (every option, every summary)

- **C:** −0.29 undiagnosed, +0.05 diagnosed never treated, +0.10 on ART, +0.14 interrupted, +4.5k new infections.
- **V:** undiagnosed +0.21 to +0.22, on ART −0.10 to −0.11, interrupted −0.10 to −0.11.
- **RN:** undiagnosed −0.27 to −0.30, on ART +0.20 to +0.22, interrupted +0.07 to +0.08. V and RN still cancel in undiagnosed but not in on ART.
- **D:** no share effect (|Δ| ≤ 0.023).
- **The S step itself:** share Δs ≤ 0.03; it carries the take-off loss (−0.18 to −0.28) and −17k to −24k new infections.

## Not robust

- **New infections at V, RN and D** (above).
- **D's infection excess:** most of option 0's +17.5k at D comes from HEEE's failed and late epidemics, not from EMOD′'s demography. With 5× seeding at both ends, D is +0.1k.
- **Conditioning removes the failed runs but not the late starts.** Taken-off HEEE has 27.8k new infections against HEEE-c's 43.2k and HEEE-S's 46.5k. On HIVsim's network the late-start effect is absent: taken-off HHEE has 69k against 67k.
- **EMOD′'s own counts:** EEEE-c has 43.3k new infections against EEEE's 40.3k, and prevalence 6.0% against 5.4%; shares the same within 0.001.

## Caveats

- On EMOD′'s network, (b) uses the `RISK_PREV` mapping of HIVsim's rule, not the rule itself.
- (b)'s seeds get EmodHIV's time since infection (fast-forwarded as EMOD's), since the natural history belongs to V.
- (b) seeds at the 1985 start, as option 0 does in the hybrids; EEEE seeds in 1980.9.
