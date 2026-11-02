# Data for HIVsim Malawi

`../malawi.py` reads these files. It stops with a list of any required file that is missing.

## Required

| File | Format | Source |
|---|---|---|
| `malawi_age_1985.csv` | `Time,Sex,AgeStart,Value`: population by single year of age in 1985, `Sex` = `Both sexes` | UN WPP |
| `malawi_asfr.csv` | `Time,AgeStart,Value`: births per 1,000 women by single year of age, per year | UN WPP |
| `malawi_deaths.csv` | `Time,Sex,AgeStart,Value`: annual death probability by sex (`Male`/`Female`) and age, per year, excluding AIDS deaths, which HIVsim models itself | `make_data.py`: UN WPP all-cause rates (`malawi_deaths_all_cause.csv`) with AIDS deaths removed by `sti.data.data_utils.dedup_deaths`, anchored at 1985 and 2025 |
| `init_prev_hiv.csv` | `risk_group,sex,sw,init_prev`: HIV prevalence at the start, by risk group (0-2), sex (`female`/`male`) and sex work (0/1) | Set to give prevalence near 1% in 1985; adjusted in calibration |
| `n_art.csv` | `year,n_art`: number of people on ART, per year | `make_data.py`: UNAIDS from 2010, Optima's data (`files/Malawi_20220912.prj`) before |

The STIsim repo has examples of the first three in `tests/test_data/zimbabwe_*.csv`.

## Optional

| File | Format | Used for |
|---|---|---|
| `condom_use.csv` | `partnership,<year>,<year>,...`: condom use by partnership type `"(i,j)"` (risk groups of the two partners), per year | the sexual network; default is STIsim's |
| `n_vmmc.csv` | `year,n_vmmc`: number of men circumcised (VMMC) per year | VMMC; none if absent |
| `vls.csv` | `year,value`: share of people on ART who are virally suppressed | ART suppression target; none if absent |

## Calibration targets (P3)

Not read by `malawi.py`. The target set is to be agreed first. It will probably include adult (15-49) prevalence, PLHIV, new infections, HIV deaths and people on ART, from the data that Optima and EMOD were fitted to.
