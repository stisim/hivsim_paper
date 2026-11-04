# 03: EMOD′ HIV transmission and natural history against EMOD

**Question:** with no testing or treatment, does EMOD′ (`models/emod/v4/`: demography, network, Risk and sex work, and now `hiv.py`) give the same epidemic as EMOD? That covers per-act transmission with EMOD's multipliers, the STI flag, traditional circumcision, transmission at birth, the seeding, and the untreated course of infection.

**What's run:**
- **EMOD, reduced to what EMOD′ has:** sets 0, 100 and 200, with the campaign cut down to sex work, STI co-infection, traditional circumcision and the seeding. Testing, care, ART, VMMC, PrEP and PMTCT are all dropped. The run reports transmissions (`TransmissionReport.csv`), deaths (`HIVMortality.csv`), and infections by HIV stage in `ReportHIVByAgeAndGender` from 1960:
  `models/emod/reference/run_emod.py --sets 0 100 200 --cores 3 --out nh3 --natural-history`
- **EMOD′:** `harness.make_sim(pid, seed, hiv=True)` for the same 3 sets, 5 seeds each.
- `run.py` reduces both to the same tables (`outputs/`), and `run.py --plot` draws `figures/`.

**Compared, over time:** prevalence by age, sex and node; new infections; people with HIV; HIV deaths; the HIV stage of people with HIV; years from infection to death from HIV; transmissions by relationship type and by the stage of the partner who transmitted.

**Caveat:** this EMOD is not the published model. Without ART, the epidemic is larger and runs on, so it is a check of the mechanisms, not of the fit. EMOD's random numbers can't be matched, so runs agree only in distribution. EMOD has one run per set, and EMOD′ five.
