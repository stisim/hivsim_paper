# 18: Does our EMOD reproduce the published EMOD in Figure 3A?

**Question:** before decomposing EMOD vs Optima HIV in Bansi-Matharu et al. 2025 Figure 3A (Malawi; `docs/core_question.md`), does our rerun of Daniel's EMOD reproduce the published EMOD panel: the share of new infections by the source partner's treatment status, and the total (30,582 new infections in 2024, caption)? Which EMOD input files do: the PR #1 files as merged (analysis 09: interrupted 50% in 2020, 54% in 2024, against about 17–20% in the figure), or with Daniel's fixes from PR #4 (`origin/mwi_demo`, not merged)?

**Figure 3A digitised:** `digitise_fig3a.py` writes `common/fig3a_malawi.csv` (EMOD and Optima HIV panels, 2000–2040, four shares a year). Each bar is sampled down its centre and each pixel assigned to the nearest legend colour; shares are good to about ±1–2 points (1 point ≈ 1.9 px). Against the values quoted in the text it reads 0.6–0.8 points off (EMOD 2024 undiagnosed 32.6 vs 33.3, on ART 45.5 vs 46.1; Optima 2024 on ART 35.6 vs 34.8, interrupted 8.5 vs 8.4).

**Categories:** the paper's four onto analysis 09's (`common/transmission.py`):

| Paper | Analysis 09 |
|---|---|
| Undiagnosed | `undiagnosed` (never diagnosed, never on ART) |
| Diagnosed but treatment naive | `diagnosed_naive` (diagnosed, never on ART; includes pre-ART care and those lost before ART) |
| Receiving treatment | `art_first_lt6m`, `art_reinit_lt6m`, `art_supp`, `art_nonsupp` |
| Not receiving treatment (after having started) | `off_first`, `off_later_lt6m`, `off_later_ge6m` |

The paper: "diagnosed and ART-naive (ie, never started treatment)", "having ART interrupted at the time of transmission", so the mapping is one to one.

**Variants (EMOD binary, emod-hiv 2.37.0; input files in a scratch copy of `models/emod/reference/source/InputFiles`, read with the new `run_emod.py --inputs`; the repo's source files are not changed):**

| Variant | Inputs | Runs |
|---|---|---|
| `bug` | PR #1 files as in the repo | analysis 09's `raw_results/emod/tr10x5` (reused) |
| `mort` | PR #4's `Demographics_Malawi.json`: mortality `ResultScaleFactor` 2.74e-6 → 2.74e-3 | `fig3a_mort` |
| `drop` | PR #4's campaign: at `OnART3` (dropout) `InterventionStatus` → `None`, so dropouts can re-enrol | `fig3a_drop` |
| `both` | both files | `fig3a_both` |

Each: `run_emod.py --transmission --seeds 1 2 3 4 5 --sets 0 25 … 225 --inputs <dir> --out fig3a_<v>`, to 2025, as analysis 09 (10 sets × 5 seeds). The calibration (the 250 resampled sets) is not redone, so the fixed variants are the same parameter sets with corrected inputs.

**Totals (national):** EMOD runs at 1% (`x_Base_Population` 0.01, as `optim_script.py`). Daniel's post-processing (`source/Analysis/mwi_transmission.Rmd`) scales by 17,563,749 (census 2018) over the sim's whole population at 2018.5, and counts `Newly Infected` from `ReportHIVByAgeAndGender` at ages 15–64. Both that and the plain ×100 are reported, for 2020 and 2024, median over runs, with sexual transmissions (all ages) for comparison.

**Compared:** for each variant, 2020 and 2024, the four shares (median and 5–95% over 50 runs) against the figure, and the totals against 30,582 (2024).

**Expected:** the dropout fix lets people restart ART, so interrupted should fall and on ART rise towards the figure (17% / 41% in 2020); the mortality fix shrinks the population (about half by 2008) and so the scaled total less than the raw one.

`run.py` classifies the runs into `outputs/`; `run.py --plot` draws `figures/`.
