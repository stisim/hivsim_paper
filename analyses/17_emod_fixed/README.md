# 17: EMOD′ with its copied bugs fixed

`docs/swap_design.md` D10 (Robyn: a fixed variant of EMOD′ is needed, as extra steps on EMOD paths). EMOD′ copies the published EMOD-HIV Malawi model's bugs (decided 3 Nov), so it reproduces Daniel's published model. This analysis gives each bug that matters for the outputs an off-by-default fix, exposes "EMOD′ fixed" to the runner as a model, and measures what each fix does.

**Question:** what does each fix, and the fixed variant, do to the outputs and to the fit to the common targets, with the calibrated parameters left as they are?

**The fixes** (`models/emod/v4/inputs.py`, `FIXES`; `components.pars(pid, fixes)`, `harness.make_sim(..., fixes=)`). The first two are Daniel's own fixes (hivsim_paper PR #4, branch `mwi_demo`, not merged), reproduced as edits of the input files as he made them:

| Fix | Bug (`docs/emod_inventory.md`, Problems found) | What the fix does |
|---|---|---|
| `mortality` | 2: background mortality 1000 times too low | `ResultScaleFactor` of `MortalityDistributionMale/Female` 2.7397e-06 → 2.7397e-03 (1/365: deaths per person per year), as in the PR |
| `reenrol` | 1: the dropout trap (from 2009 dropouts stay OffART) | Daniel's new campaign event: from 1990, on `OnART3` (ART dropout) a `PropertyValueChanger` sets `InterventionStatus` to `None`, so the ART-start choices at `OnART1` are no longer disqualified. In the condensed cascade this is `dropout.p_restart` = 1 |
| `art_male` | 3: male non-suppressive ART mortality, 182–365-day row about 7× | that row (CD4 ≥ 25) as the female table has it |
| `risk_nodes` | 5: the Risk overlay's LOW fractions on the wrong nodes | each region's fraction on its own node |
| `anc` | 6: the ANC testing uptake table is ignored | `anc.p_test` from the table (condensed cascade only) |
| `pmtct` | Daniel note 3: combination PMTCT blocked by ART staging | `pmtct.blocked_if_staged` = 0 (condensed cascade only) |

Not fixed: the calibrated parameters with no effect (4: which value was meant is unknown), the ART-mortality draw rule (7: EMOD behaviour), the 1% truncation of NChooser targets (8: a scale choice), test-and-treat from 2006/7 (9: a modelling choice), and PrEP's exclusion (about 10 agents).

**The variant:** `models/emod/v4/components_fixed.py` is EMOD′ with Daniel's two fixes (`mortality`, `reenrol`). It re-exports EMOD′'s factories, so the runner and `own()` treat it as EMOD′, and each slot's fixes come with its parameter set: `mortality` is in D's, `reenrol` in C's. So a path `{'from': 'emod', 'to': 'emod_fixed', 'slots': ['D', 'C']}` with `models: {'emod_fixed': 'models/emod/v4/components_fixed.py'}` swaps "bug fixed" in as steps (D: mortality; C: re-enrolment). For the single fixes this analysis has one small factory file each (`variants/`).

**Checks first:**
- Default unchanged: `build(all EMOD′)` from the runner, set 0 seed 1 to 2025, every result array and analyzer table identical to the parent commit (`tests/` and `common/swaps/tests` pass).
- `reenrol` in the condensed cascade against Daniel's event run by the literal interpreter (`emod_cascade.py`): OffART share and ART coverage, 5 seeds each, set 0.

**Runs** (`run.py`, env `hp_rc3`, `PYTHONPATH=/tmp/sti-rs-v4-art`): every variant (EMOD′, each fix, Daniel's two = `emod_fixed`, all six) on its own (all slots), condensed cascade, parameter sets 0, 25 and 50, 10 seeds, stop 2025, through the runner (`swaps.runner.start`: resumed from each hybrid's 1985 burn-in, stored under `raw_results/swaps/`), with OffART read from `EmodReport` as well.

**Outputs** (by year, mean and 10–90% over seeds): population against WPP, PLHIV, ART coverage, OffART share of PLHIV, the share of transmissions from people who interrupted ART (MIHPSA Malawi across models: 8.4–20.1%, analysis 09) and the full split by cascade stage (`common/outputs.py`), HIV deaths, share of HIV deaths on ART, and the score against `common/targets.csv` (`common/targetfit.py`, primary and secondary). Take-off (D11) is reported, not conditioned on, as these are EMOD′'s own D.

**Caveat:** EMOD's parameters were calibrated with the bugs in, so fixing them un-calibrates it. The fixed variants are reported uncompensated; the SUMMARY says what a refit would need. Daniel believes (`emod_responses.md` in PR #4) his original calibration used the correct mortality factor and that the dropout trap came in later, with the move to `ARTMortalityTable`; if so, the fixed variant is closer to the model he calibrated than EMOD′ is.

**What would count as done:** the default is identical; `reenrol` agrees between the condensed and the literal cascade; the effect table per fix with uncertainty over seeds.
