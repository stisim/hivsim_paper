# 07: HIVsim with EMOD's condensed cascade

**Question:** can the condensed cascade (`common/cascade/`) drive HIVsim (STIsim's `sti.HIV`) through the same `CascadeHost` interface it drives EMOD′ with, and what changes when HIVsim's own care (coverage-forced testing and ART, VMMC to survey coverage) is swapped for EMOD's mechanistic cascade with Malawi's EMOD parameters (set 0)? This is the first component swap.

**What's run:**
- Starsim v4 from `~/starsim` (commit ebd52c1a, 4.0.0rc2) and STIsim branch `rs-cascade-host` (from `rs-art-split` c815b0a; adds `HIV.who_stage()`, `HIV.new_symptoms()`, `HIV.start_pmtct()`, `ART.start_art(effective=)`, and PrEP courses that HIV ends), on `PYTHONPATH`.
- HIVsim Malawi from hivsim_malawi commit 803b2f7 (`models/hivsim/malawi/`, copied unchanged), with the best run of its experiment 03 (calibration is still under way, so this is a provisional parameter set), 10,000 agents, 1985-2022.
- Three variants (`models/hivsim/swap.py`), 50 seeds each (the condensed cascade without PMTCT, 20):
  - `own`: HIVsim as calibrated: general, low-CD4 and ANC testing, ART to the programme's numbers on ART (`n_art.csv`) with suppression to MPHIA's share (`vls.csv`), VMMC to survey coverage by age.
  - `condensed`: HIVsim's testing, ART and VMMC removed; the condensed cascade with EMOD's Malawi table (`common/cascade/tables/malawi_emod_set000_*.csv`), its VMMC and PrEP flows summed over EMOD's nodes and scaled by HIVsim's `pop_scale`; `sti.ART` kept only to hold ART's states (it starts and stops no one); circumcision before the programmes held at HIVsim's pre-2010 survey level (as EMOD keeps traditional circumcision outside the cascade).
  - `condensed_nopmtct`: as `condensed`, with PMTCT efficacies 0, to check that PMTCT reaches HIVsim.

**Measures:**
1. Sanity check of the host (`condensed`): the cascade's events against HIVsim's (12-week pregnancies against conceptions, 6-week-olds against births, debuts, first symptoms), consistency of states (on ART if and only if in care on ART; everyone on ART or in staging diagnosed), ART starts and stops against the cascade's deliveries, the suppressive share against the table's arms, PMTCT (mothers covered, and infections in children with and without it).
2. The swap, `own` against `condensed`: the common targets (`common/targets.csv`, measured by hivsim_malawi's `targetfit`), and the cascade over time (PLHIV, on ART, diagnosed, ART coverage, off ART after starting it, suppressed on ART, new infections, HIV deaths, infections of children, circumcision), with EMOD′'s condensed runs from analysis 06 as a reference for the shares.

`run.py` runs and writes `outputs/`; `run.py --plot` draws `figures/`.

**What to expect:** the two HIVsim variants share HIVsim's biology and network, so differences come from care: HIVsim's ART follows the programme's numbers by construction, while the condensed cascade starts and stops ART mechanistically, with EMOD's dropout trap (from 2009 most who drop out never return), so ART coverage should be lower and HIV deaths higher with the condensed cascade. EMOD's parameters were fitted with EMOD's disease, not HIVsim's, so neither variant is expected to fit the targets as well as its own calibration.
