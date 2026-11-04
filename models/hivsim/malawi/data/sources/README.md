# Sources for the calibration targets

Read by `../../make_targets.py`, which writes `results/targets.csv`.

| File | What it is |
|---|---|
| `emod_MalawiNationalPrevalence.csv` | EMOD's Malawi prevalence targets, copied from `hivsim_paper's models/emod/reference/calib_targets_data/` (PR #1). By sex and 5-year age group: DHS 2004 (`Year` 2004), DHS 2010 (2010), MPHIA 2015-16 (2015, ages to 64) and DHS 2015-16 (2015.5), with beta-distribution effective counts and 95% intervals |
| `emod_Malawi_Incidence_new.csv` | EMOD's incidence targets: MPHIA 2015-16 by sex and age, and MPHIA 2020-21 for 15-49 by sex |
| `emod_Malawi_OnART.csv` | EMOD's numbers on ART by sex, ages 15-64, 2004-2015 (programme data via Spectrum) |
| `mphia_2020_21.csv` | Typed from the MPHIA 2020-21 report (`resources/241122_Mphia_Foreword.pdf`): Table 6.3 (prevalence by sex and age, with the number tested), Table 5.1 (incidence, 15-49), Table 9.1.B (conditional 95-95-95, 15+) |
| `optima_fsw_prevalence.csv` | Prevalence among female sex workers from Optima's data (`hivprev`, FSW, in `hivsim_paper's models/optima/reference/Malawi_20220912.prj`); bounds only for 2019 |
