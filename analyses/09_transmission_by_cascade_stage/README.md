# 09: Transmissions by the source's cascade stage, EMOD against EMOD′

**Question:** what share of new infections comes from sources who were undiagnosed, diagnosed but never on ART, on ART (by time since starting and by arm), or off ART after starting (interrupted), year by year, in EMOD and in EMOD′ with the condensed cascade; and do the two agree within run-to-run noise? This is the policy endpoint "transmissions by the source's cascade stage" (`docs/swap_design.md` 2.2), and the share from interrupted treatment is what the treatment-interruption scenario acts on.

**Categories:** Daniel's (`models/emod/reference/source/Analysis/mwi_transmission.Rmd`), made into one partition of the transmissions, so they add up to all of them; defined in `common/transmission.py` and below. For a transmission at time t, by the source's history up to t (events in the same time step count as before it, because in both models the cascade runs before exposure in a step):

| Category | Definition |
|---|---|
| `undiagnosed` | Never diagnosed and never on ART |
| `diagnosed_naive` | Diagnosed, never on ART |
| `art_first_lt6m` | On ART, in the first ART episode, for less than 6 months |
| `art_reinit_lt6m` | On ART, in a second or later episode (restarted after stopping), for less than 6 months |
| `art_supp` | On ART for 6 months or more in this episode, on the suppressive arm |
| `art_nonsupp` | On ART for 6 months or more in this episode, on the non-suppressive arm |
| `off_first` | Off ART after one episode (first interruption) |
| `off_later_lt6m` | Off ART after two or more episodes, stopped less than 6 months ago |
| `off_later_ge6m` | Off ART after two or more episodes, stopped 6 months ago or more |

Starting ART counts as being diagnosed. Coarse groups for the paper: undiagnosed; diagnosed, never on ART; on ART under 6 months; on ART suppressive arm; on ART non-suppressive arm; interrupted (the three `off_*`). Neither model has viral load: "suppressed" is the ART arm (EMOD's `OnEffectiveART`, transmission × 0.04 after a 183-day ramp; non-suppressive × 0.65). Daniel's viral-load columns were a fixed 92/8 split of on-ART counts, not model output, so they are not reproduced. His CD4-at-initiation split of the under-6-month categories is not reproduced either.

**What's run:**
- **EMOD:** `models/emod/reference/run_emod.py --transmission --seeds 1 2 3 4 5` for parameter sets 0, 25, …, 225 (the 10 sets of analyses 04–06), to 2025: 50 runs, emod-hiv 2.37.0. `--transmission` turns on Daniel's reports (`config_transmission.json`): `TransmissionReport` (each sexual transmission with source and recipient IDs and their individual properties, including `ARTStatus`), `ReportEventRecorder` (here `HIVNewlyDiagnosed`, `StartedART`, `StoppedART`, all ages), `ReportHIVART` and `ReportHIVMortality` (not `ReportHIVInfection`: a row per infected person per step, about 1.5 GB a run). Daniel's start year (2016) and minimum age (15) are changed to 1960 and 0 so each source's history is complete. Output only.
- **EMOD′:** `harness.make_sim(pid, seed, cascade=True, condensed=True)` with the recorder `common.transmission.TransmissionLog`, the same 10 sets × 5 seeds: 50 runs. Starsim v4 at ebd52c1a (`PYTHONPATH=/tmp/starsim-ebd52c1a`), as analyses 05–08.
- **Classification:** `common/transmission.py`. For EMOD (`from_emod()`), each transmission's source is joined to its diagnosis and ART start and stop events; on/off ART and the arm are read from the source's `ARTStatus` at the transmission and checked against the events. For EMOD′ (and any model later), the recorder writes one row per transmission with the source's state at that moment, from `CascadeHost.cascade_stage()` (diagnosed, on ART, suppressive arm) and an ART history it keeps by watching who is on ART each step. `classify()` is the same for both.
- Sexual transmissions only (EMOD's `TransmissionReport` has no mother-to-child transmissions). All source ages and both sexes; years are calendar years of the transmission.

**How runs compare:** in distribution, as in analyses 05–06: for each year and category, the share of that year's transmissions, median and 5–95% range over the 50 runs of each model, and the difference of means with its standard error. Years 2010, 2015, 2020 and 2024 (2025 is a partial year) are tabulated. Shares of the interrupted group are also set against MIHPSA Malawi's published range across models, 8.4–20.1% (`docs/prompts.md`).

**Expected:** EMOD′ matched EMOD within noise on every target in analyses 04–06, so the shares should agree too. EMOD's dropout trap (`docs/emod_inventory.md`, problem 1) means few people restart ART after 2009, so `art_reinit_lt6m` and `off_later_*` should be small, and `off_first` large and growing.

**Also checked:** that the categories partition the transmissions; Daniel's "diagnosed, not on ART" count, which went negative in his Rmd, recomputed his way on our runs to show why.

`run.py` runs EMOD′ and classifies both (`outputs/`); `run.py --plot` draws `figures/`.
