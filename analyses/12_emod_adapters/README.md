# 12: EMOD′ components on HIVsim's slots (adapters for steps 6 and 7)

**Question:** do EMOD′'s HIV module (EmodHIV) and its network and risk (EmodNetwork, EmodRisk) run on HIVsim's slots (STIsim's demography, and for EmodHIV also `sti.StructuredSexual`) through the adapters in `models/emod/v4/adapters.py`, and does each swapped-in component still behave as it does in EMOD′? These are the adapters of `docs/swap_design.md` §5, steps 6 and 7.

**Hybrids** (`adapters.make_sim`; parameter set 0, monthly from 1985 to 2025, 10,000 agents in 1985 as HIVsim's demography is run (Decided, 9 Nov)):
- **6 (V on HIVsim):** EmodHIV on `StructuredSexual` (with Malawi's condom data) and STIsim's demography (`ss.Pregnancy`, `ss.Deaths`, the 1985 age structure, from `models/hivsim/malawi/data`). No care cascade.
- **7 (VRN on HIVsim's D):** EmodHIV, EmodNetwork and EmodRisk on STIsim's demography; and the same without HIV, for the network checks.
- **Home:** EMOD′ (`harness.make_sim(pid=0, hiv=True)`, and without HIV), from 1960.5 at 1%.

**What's bridged, and how** (listed in full in SUMMARY):
- D → EMOD′'s people fields (`EmodPeople`, named `emoddemog`): `node` drawn from EMOD's node shares at the start, a newborn's from its mother; `emod_age` as `(age + dt)·365`; the step's births and mothers from the pregnancy's deliveries; women reaching 12 weeks of gestation.
- HIVsim's R → EMOD's (step 6): Risk HIGH while FSW or client (so EMOD's seeding of 2% of HIGH reads them), LOW for `risk_group` 0, else MEDIUM; the STI flag drawn at the network's debut with EMOD's probabilities by that Risk.
- HIVsim's N → EmodHIV's contacts (step 6): `acts` (whole acts per step, fixed per partnership) as they are; `condoms` as the chance of a condom per act (same meaning in both); types mapped for the infection log.
- Starting EMOD's R and N in 1985 (step 7): the starting population gets the sex work and STI flag it would have had from its debut; the network runs its 15-year burn-in as a warm-up at the start, with time frozen, so entry rates are set and frozen before 1985's first step.
- Numba's global random stream: EmodNetwork reseeds it from its own generator each step.

**Tests** (`models/emod/v4/test_adapters.py`, pytest): each hybrid builds and runs 1985→2025 at 2,000 agents; a run repeats exactly; the network's draws don't move when other code draws from Numba's global stream first (and do without `own_stream`, showing the risk); EmodHIV in EMOD′ is unchanged by its `contacts()` hook (same results as before it, checked once against the parent commit).

**Sanity checks** (`run.py`, 6 seeds of each; `run.py --plot`):
- EmodHIV's natural history: Kaplan–Meier survival from infection to HIV death (censored at other deaths and 2025), by age at infection (<15, 15–24, 25–34, 35+), home vs hybrids 6 and 7. Expected: the same curves within noise, since the prognosis depends only on age at infection; a difference would mean `emod_age` or deaths are read wrongly.
- EmodNetwork: in 1990, 2000, 2010 and 2020, current relationships of each type per person by sex and age group (15–24, 25–34, 35–49), new relationships per person-year by type, the man's age minus the woman's at formation by type, and the share at Risk HIGH; home vs hybrid 7 without HIV. Expected: similar but not equal. The demographies differ (EMOD's reads mortality 1000 times too low; HIVsim's has Malawi's), and EmodNetwork's entry rates are per person in each age bin, set from the population at its burn-in (1960.5–1975.5 at home, 1985 in the warm-up), so differences that follow from age structure are expected; differences inside an age group and sex that are large, or in the age gaps, would not be.
- The epidemic (not a check of fidelity): prevalence 15–49 and new infections per year in each, and how many runs die out after seeding (seeding 2% of HIGH gives few seeds at 10,000 agents).
