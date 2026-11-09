# 17: The fixed EMOD′ variant (D10) — summary

**Answer:** with Daniel's two fixes (hivsim_paper PR #4: background mortality ×1000, and re-enrolment after ART dropout) and no recalibration, EMOD′ fits the common targets far better (score 7.5 against 108) and roughly reproduces the published EMOD's Figure 3A split (Bansi-Matharu 2025; `docs/core_question.md`). The dropout trap explains the whole excess of transmissions from people who interrupted ART; the mortality bug explains the population (2× WPP) and most of the excess of new infections. The other four copied bugs are within seed noise. This fits Daniel's account that his calibration used the correct mortality factor and that the trap came in later, so `emod_fixed` is probably closer to the model he calibrated than the copied-bug EMOD′.

**Switches** (`models/emod/v4/inputs.py`, `FIXES`, all off by default; `components.pars(pid, fixes)` or `harness.make_sim(..., fixes=)`):

| Switch | Bug | How |
|---|---|---|
| `mortality` | background mortality 1000× too low | as PR #4: mortality `ResultScaleFactor` 2.7397e-06 → 2.7397e-03 |
| `reenrol` | dropout trap | as PR #4: from 1990, on `OnART3`, `PropertyValueChanger` sets `InterventionStatus` to None; the condensed cascade reads it as `p_restart` = 1 |
| `art_male` | male non-suppressive ART mortality row (182–365 days, CD4 ≥ 25) ~7× too high | set to the female values |
| `risk_nodes` | Risk overlay's LOW-risk fractions on the wrong nodes | input-file edit |
| `anc` | ANC uptake table ignored | condensed cascade only |
| `pmtct` | Option B blocked for women sent to staging | condensed cascade only |

- **Default unchanged:** `build(all EMOD′)`, set 0, seed 1, to 2025, equals the parent commit (da52324) in all 162 result arrays and every analyzer table. `tests/` and `common/swaps/tests` pass (19/19).
- **Translation of `reenrol` checked:** the literal campaign interpreter with Daniel's event agrees with the condensed cascade (5 seeds each; 2020 OffART share 0.069 in both, ART coverage 0.795 vs 0.802; `outputs/literal_vs_condensed.csv`).
- **Runner:** `models/emod/v4/components_fixed.py` is EMOD′ with Daniel's two fixes (`models: {'emod_fixed': 'models/emod/v4/components_fixed.py'}`). It re-exports EMOD′'s factories, and the fixes travel with each slot's parameter set (`mortality` with D, `reenrol` with C), so a path `emod → emod_fixed` over D and C gives "bug fixed" as two steps (not yet run). Single-fix variants are in `variants/`.

**Runs:** 9 variants × Daniel's sets 0, 25, 50 × 10 seeds = 270 runs through the runner, each resumed from its own 1985 burn-in; no errors; every run took off. Outputs: `outputs/effects.csv`, `fig3a.csv`, `yearly.csv`, `scores.csv`, `literal_vs_condensed.csv`, `variants_set{000,025,050}.png`; `outputs/fig3a_published.csv` is a copy of `common/fig3a_malawi.csv`.

## Effects, set 0, 2020 (mean of 10 seeds; counts are agents × 100, not rescaled to the census)

| Variant | Pop / WPP | PLHIV | ART cov | OffART / PLHIV | Tx from interrupted | HIV deaths | Deaths on ART | Score (all / primary) |
|---|---|---|---|---|---|---|---|---|
| EMOD′ (bugs copied) | 1.99 | 2.26M | 0.41 | 0.42 | 0.51 | 109k | 0.22 | 108 / 102 |
| fix mortality | 0.91 | 0.97M | 0.41 | 0.42 | 0.52 | 47k | 0.21 | 55 / 99 |
| fix reenrol | 2.00 | 2.16M | 0.80 | 0.07 | 0.13 | 48k | 0.77 | 29 / 12 |
| fix art_male / risk_nodes / anc / pmtct (each) | ~1.99 | 2.2–2.3M | 0.41 | 0.42 | 0.50–0.51 | 106–112k | 0.20–0.23 | 106–114 / ~103 |
| **Daniel's fixes** (`emod_fixed`) | 0.92 | 0.94M | 0.80 | 0.07 | 0.13 (10–90%: 0.11–0.15) | 20k | 0.75 | 7.5 / 10.9 |
| all six | 0.91 | 1.02M | 0.81 | 0.07 | 0.13 | 21k | 0.74 | 7.5 / 11.0 |
| Data | 1 (WPP 19.5M) | 1.0M | 0.84 | | 0.08–0.20 (MIHPSA, across models) | 16k | 0.35–0.64 (across models) | |

Sets 25 and 50 give the same picture: `emod_fixed` scores 7.0–7.9 against 103–107 with the bugs copied.

## Figure 3A (undiagnosed / diagnosed never treated / on ART / interrupted; new infections, national)

| Year | Published EMOD (digitised) | EMOD′ (bugs copied) | Daniel's fixes, set 0 |
|---|---|---|---|
| 2020 | 0.385 / 0.037 / 0.406 / 0.171 | 0.319 / 0.029 / 0.144 / 0.508; 162k | 0.356 / 0.034 / 0.481 / 0.129; 38k (35–43k) |
| 2024 | 0.326 / 0.021 / 0.455 / 0.198; 30,582 | 0.297 / 0.020 / 0.128 / 0.554; 191k | 0.323 / 0.025 / 0.504 / 0.148; 39k |

- `reenrol` alone gives 0.359 / 0.039 / 0.473 / 0.129 in 2020 (94k); `mortality` alone leaves the shares where EMOD′ has them.
- Sets 25 and 50 with Daniel's fixes: 0.36 / 0.05–0.06 / 0.44–0.45 / 0.14 in 2020.
- The fixed variant has 2–4 points too few infections from undiagnosed and interrupted sources and 5–8 points too many from people on ART.
- The published total is rescaled to the census; rescaled the same way (×1/0.92), ours is about 42k in 2024 against the published 30.6k.

## Open

1. **`emod_fixed` is uncalibrated.** New infections after 2005 are too high (38k vs UNAIDS 21k in 2020) and fall too slowly; HIV deaths 20k vs 16k; deaths on ART 75%, above the 35–64% across MIHPSA models. A refit would rerun Daniel's calibration (`optim_script.py` and his form) with both fixes in the inputs; the census rescale would then be ~1, so population can carry real weight; it would free the transmission, network and linkage/dropout parameters, and draw a fresh posterior rather than reweighting the old resampled sets (129 distinct vectors among 250).
2. Daniel's PR #4 is not merged; `reference/source` is untouched. Which files produced the published results is asked of Daniel (9 Nov).
3. `anc` and `pmtct` work in the condensed cascade only; the compiled cascade hasn't been tried with the new event.
4. The `emod → emod_fixed` path over D and C hasn't been run.
