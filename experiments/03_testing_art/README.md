# Exp 03 — Can testing by sex and over time get ART coverage and viral suppression right?

**Question.** After [exp 02](../02_condoms_and_circumcision/SUMMARY.md) the later epidemic is close (best run within 4-16% of UNAIDS PLHIV and new infections from 2010), and Robyn wants ART coverage and viral suppression to match well. In exp 02 the ART input was not reached in many runs (median 82% of it in 2020), because too few people were diagnosed, and men on ART were below the data in nearly every run. Testing was the same for men and women and scaled up from 2004 in every run. If men's testing and the timing of the scale-up can vary, do runs reach the ART input and match ART coverage and suppression, by sex?

**Plan.**
- New targets, all primary (`results/targets.csv`, added at the end so earlier rows keep their numbers): ART coverage (share of PLHIV on ART) from UNAIDS, 2010-2021, all ages; and from MPHIA 2020-21 (Table 9.1.A, 15+), ART coverage and viral suppression among all PLHIV, by sex. Suppression among people on ART is still an input (`vls.csv`), so suppression among PLHIV tests diagnosis and ART together.
- Model change: general testing is now a function of sex and year (`make_testing` in `models/hivsim/malawi.py`): it rises linearly from `test_start` to 2020, and men's is `test_rel_m` times women's. Low-CD4 and ANC testing are as before.
- Ranges: exp 02's ten, plus `test_rel_m` 0.5-3 (log) and `test_start` 1995-2010. 200 parameter sets and 10 seeds of the central set, 10,000 agents, 1985-2022.
- Plots from `models/hivsim/targetfit.py`, with a new figure for ART coverage and the MPHIA cascade; the trajectory checks now start in 2000 and include ART coverage (within 10 percentage points of UNAIDS every year).

**Success criteria.** Some runs are within 10 points of UNAIDS ART coverage every year 2010-2021, within MPHIA's intervals for ART coverage and suppression among PLHIV for both sexes, and reach the ART input (over 95% of it) through 2020, while keeping exp 02's later epidemic. If the ART input is reached but men stay low, the split between the sexes needs more than testing (for example linkage or loss from care, `sti.HIVCare`).
