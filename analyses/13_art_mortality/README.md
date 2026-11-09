# 13: HIVsim's mortality on ART

**Question:** is HIVsim's mortality on ART too low, and why? In the Malawi calibration (`hivsim_malawi`, exp 03), runs with the right ART coverage have AIDS deaths collapsing after 2015 (best run: 0.10× UNAIDS in 2021). `docs/plan.md` lists this as a model gap to fix in STIsim before the swaps, since deaths on ART are a paper output. This analysis measures it; it changes no STIsim code.

**Code:** STIsim PR #612 (`rs-v4-art-split`, 6058d5c), where ART is a declared module (`sti.ART`: `art_death`, `dropout`, `post_art_death`; CD4 on treatment in `ART.update_cd4`), with Starsim from `/tmp/ss-nov09`, env `hp_rc3`:

    PYTHONPATH=/tmp/sti-g:/tmp/ss-nov09 conda run -n hp_rc3 python run.py [hivsim|malawi|optima|emod|plot]

**Method:**

1. **HIVsim's death rate on ART, by CD4 at ART start and time on ART** (`run.py hivsim`). A plain `hivsim.Sim` (default modules, 2000-2030, monthly, 30,000 agents, 6 seeds) with ART coverage of PLHIV rising linearly 0 (2004) → 30% (2010) → 60% (2015) → 80% (2020) → 85% (2030), and testing (general 20%/yr, CD4 < 200 80%/yr) so that the coverage target can be met. An analyzer follows everyone from their first ART start: person-time and deaths each step by state (on suppressive ART, on non-suppressive ART, off ART after starting it), CD4 at first start (<50, 50-200, 200-350, 350-500, 500+), and time since first start (0-6 m, 6-12 m, 1-2 y, 2-5 y, 5 y+); and the cause of each death (`art_death`, `post_art_death`, `hiv_death`, background). Also the mean CD4 by time on ART, and HIV deaths by year split into never treated, on ART and off ART after starting.
2. **Malawi** (`run.py malawi`). The exp 03 best run (`hivsim_malawi` `make_sim`, its parameters and seed 19, 50,000 agents, 6 seeds) with the same analyzer, against UNAIDS AIDS deaths 2010-2021; then variants that test the candidate causes (chosen after step 1, listed in SUMMARY), to give the expected effect of each on deaths 2015-2021.
3. **Comparators.** Optima′ (`run.py optima`): HIV deaths per person-year by cascade state and CD4 band, from Optima's calibrated Malawi run (`models/optima/reference/optima_malawi_ref.obj`, `raw.death` and `raw.people` by state). EMOD′ (`run.py emod`): its `ARTMortalityTable` rates (from the Malawi campaign file), by time on ART and CD4 at ART start. Literature: IeDEA / ART-CC cohort mortality on ART in southern Africa, cited in SUMMARY.

**Expected:** HIVsim's on-ART rate is the off-ART rate at the *current* CD4 count × 0.25 (suppressed), and CD4 recovers within about two years, so after the first year the rate should fall to about 0.1%/yr, well below cohort estimates (~1-2%/yr in later years). If so, the gap is structural (mortality tied to recovered CD4, capped by the off-ART rate), not a single parameter.
