# 18: Does our EMOD reproduce the published EMOD in Figure 3A? — summary

**Answer:** the shares match with Daniel's dropout fix (PR #4's campaign: `InterventionStatus` → `None` at `OnART3`), and they don't match without it. With the fix, every 2020 and 2024 share is within 5 points of the digitised figure, and each figure value sits inside our 5–95% range over runs. On ART is 3–5 points high and interrupted 3–5 points low. With the PR #1 files as merged, interrupted is 50% in 2020 (figure 17%) and on ART 14% (figure 41%). The mortality fix (`ResultScaleFactor` ×1000) leaves the shares as they are, but it changes the population: at 1% the sim's 2018 population is then the census (scaling factor 105 instead of 48). The total comes out too high in every variant. With both fixes, it is 41,300 new infections at ages 15–64 in 2024 (median; 5–95% over runs 26,500–59,400), against the caption's 30,582. That is 35% high, though inside the run-to-run range, and one of the 10 sets (175) gives 30,600. With the PR #1 files it is 79,100, 2.6 times the caption. **For Figure 3A, use the both-fixes files** (`raw_results/emod/fig3a_both`): the shares come from the dropout fix, and the mortality fix makes the population and the scaling right. The total is still about a third high. It isn't explained here; it could be the parameter sets (10 of the 250) or the counting, see Caveats.

Runs: emod-hiv 2.37.0, the 10 sets of analyses 04–09 × seeds 1–5, to 2025, `run_emod.py --transmission --inputs <copy of InputFiles with PR #4's file(s)>` (new `--inputs DIR` option; the repo's source files are unchanged). 150 new runs, none failed.

## Figure 3A digitised

`common/fig3a_malawi.csv` (model, year, category, share; EMOD and Optima HIV, 2000–2040; caption totals in `#` header lines: read with `pd.read_csv(path, comment='#')`). Each bar is read down its centre by pixel colour (`digitise_fig3a.py`); precision ±1–2 points (shares under ~1.5 are edge blur). Against the values quoted in the paper's text, the readings are 0.6–0.8 points off.

Category mapping onto analysis 09's groups: undiagnosed = `undiagnosed`; diagnosed never treated = `diagnosed_naive` (incl. pre-ART care and people lost before starting ART); on ART = `art_lt6m` + `art_supp` + `art_nonsupp`; interrupted = the three `off_*` groups.

## Shares (%) of new infections by the source's status, mean over 50 runs (5–95%), against Figure 3A (`outputs/shares_vs_figure.csv`)

| Year | Category | Figure 3A EMOD | PR #1 files | Mortality fix | Dropout fix | Both fixes |
|---|---|---|---|---|---|---|
| 2020 | Undiagnosed | 38.5 | 32.0 (26.7–36.7) | 31.8 | 37.0 (31.4–41.9) | 36.7 (29.8–42.8) |
| | Diagnosed, never treated | 3.7 | 3.2 (2.2–4.3) | 3.2 | 4.2 (2.3–6.0) | 4.4 (2.5–7.0) |
| | On ART | 40.6 | 14.3 (13.1–15.6) | 14.2 | 45.4 (41.4–49.5) | 45.0 (40.0–50.7) |
| | Interrupted | 17.1 | 50.5 (46.8–56.1) | 50.7 | 13.3 (10.8–15.8) | 13.9 (11.1–17.7) |
| 2024 | Undiagnosed | 32.6 (text 33.3) | 30.7 (25.2–34.7) | 29.3 | 34.3 (27.8–41.6) | 33.7 (25.7–43.4) |
| | Diagnosed, never treated | 2.1 | 2.2 (1.2–3.3) | 2.1 | 2.8 (1.5–4.4) | 2.9 (1.3–5.1) |
| | On ART | 45.5 (text 46.1) | 12.6 (11.0–14.6) | 12.5 | 48.1 (43.2–52.5) | 48.3 (42.3–54.3) |
| | Interrupted | 19.8 (text ≤ 20.1) | 54.5 (49.3–60.8) | 56.1 | 14.8 (11.9–18.9) | 15.1 (11.0–20.1) |

Over time (`figures/shares_vs_figure.png`), the dropout-fix variants follow the figure's shape from 2000 to 2024. In 2010–2015 they are a few points short on undiagnosed (2015: 45 vs 49) and on interrupted (10 vs 13), and over on on-ART (34 vs 28). The PR #1 and mortality-fix variants diverge from 2012, when interrupted starts to take over. Diagnosed, never treated agrees in every variant (it comes before ART, so the dropout trap doesn't touch it).

## New infections, national (median over 50 runs; `outputs/totals_by_variant.csv`)

Daniel's scaling (`mwi_transmission.Rmd`): × 17,563,749 / sim population (all ages) at 2018.5 (our report step 2018.42); `Newly Infected` at ages 15–64. In 2024 the runs stop at 2025.0 and miss the last step(s), so 2024 is annualised.

| Variant | Scale factor | 2020 | 2024 | 2024 ×100 (1% sim, unscaled) |
|---|---|---|---|---|
| PR #1 files | 48.4 | 68,000 | 79,100 (59,400–113,000) | 163,300 |
| Mortality fix | 105.3 | 65,800 | 70,200 | 66,100 |
| Dropout fix | 48.1 | 43,400 | 46,500 (31,800–65,900) | 95,000 |
| Both fixes | 104.9 | 39,700 (29,800–55,500) | 41,300 (26,500–59,400) | 39,600 |
| Figure 3 caption | | | 30,582 | |

With both fixes, the scaled and unscaled totals agree because the population is now right. Sexual transmissions (`TransmissionReport`, all ages) give the same numbers within 2%. By parameter set (both fixes, 2024, mean of 5 seeds), the totals run from 30,600 (set 175) to 60,500 (set 225).

## Caveats

- Figure 3A was digitised at ±1–2 points.
- The published EMOD presumably used more parameter sets (or another selection) and possibly other report windows; Daniel's Rmd reads an experiment called `...recalib-transmission___2023_11_04`. 10 sets × 5 seeds is a sample, and the totals spread widely across sets. The 35% excess is about 1 SD of the run-to-run spread.
- The caption doesn't say which ages its total counts. Ours is 15–64 recipients, as in Daniel's Rmd; all ages would be within 2%.
- The fixed variants use the same calibrated parameter sets: nothing was recalibrated with the fixes in. That the shares match anyway supports Daniel's view that the dropout trap came in after the calibration used for MIHPSA.
- The EMOD binary and EMOD′ (`emod_fixed`, analysis 17) agree: 2020 shares 36.7 / 4.4 / 45.0 / 13.9 here vs 35.6 / 3.4 / 48.1 / 12.9 there (set 0 only).
