# Exp 02 — Do condom use and circumcision bend the epidemic down?

**Question.** In `../01_coverage_check/SUMMARY.md`, every target was inside the range of runs on its own, but no run had both Malawi's early epidemic and its decline: in the median run, new infections in 2020 were 1.45 times their 1990-95 peak, against 0.16 in the data. Exp 01 had no condom use (STIsim's default is none) and no circumcision. Do they let runs reach both?

**Plan.**
- The same nine parameters and ranges as exp 01, plus `condom_scale` (0.5-1.5, multiplying all condom use, capped at 1), since the level of condom use is uncertain and only its pattern is meant to carry over (Robyn).
- New inputs, made by `models/hivsim/make_data.py`:
  - `condom_use.csv`: condom use by partnership type and year, with the pattern of STIsim's Zimbabwe file from the sti_notification project (none before 1990, rising through the 1990s; lowest between two low-risk partners, highest in sex work), scaled so that partnerships not between two low-risk people are at 49.6% in 2020, MPHIA 2020-21's condom use at last sex with a nonmarital partner.
  - `vmmc.csv`: share of men circumcised by age, MPHIA 2020-21 Table 12.5: the nonmedical (traditional) share up to 2010, rising linearly to the total (medical and nonmedical) in 2020.5.
- 200 parameter sets (Latin hypercube) and 10 seeds of the central set, 10,000 agents, 1985-2022, as in exp 01.
- Plots and tables from `models/hivsim/targetfit.py`, now with a check of whole trajectories: the share of runs within a factor of 1.5 of each UNAIDS series in every year, and new infections in 2020 over their 1990-95 peak.

**Success criteria.** Some runs have new infections in 2020 at 0.11-0.24 times their early peak (within a factor of 1.5 of the data's 0.16), and some runs are within a factor of 1.5 of PLHIV and new infections in every year, while prevalence by age stays covered as in exp 01. If the decline still isn't reached, the missing piece is elsewhere (testing and ART too weak, the network, or the epidemic starting too late), and exps 03-04 from exp 01's SUMMARY come next.
