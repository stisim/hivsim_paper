# 12: EMOD′ components on HIVsim's slots — summary

**Answer:** both adapters work. Hybrid 6 (EmodHIV on `StructuredSexual` and STIsim's demography) and hybrid 7 (EmodHIV, EmodNetwork and EmodRisk on STIsim's demography) build and run 1985→2025. EmodHIV's natural history matches EMOD′'s on both demographies. EmodNetwork on STIsim's demography is within about 8% (median) of EMOD′'s network by type, sex and age, with the same age gaps. The exception is commercial relationships, which run 20–55% above EMOD′'s from 2000 on (open item below).

**EmodHIV (natural history; `survival_vs_home.csv`, `figures/survival.png`)**
- Kaplan–Meier survival from infection to HIV death, by age at infection (<15, 15–24, 25–34, 35–44, 45+) at 2, 5, 10 and 15 years. Pooled runs: EMOD′ (6 seeds, 314,000 infections), hybrid 6 (24 seeds, 30,000) and hybrid 7 (24 seeds, 41,000).
- 4 of the 40 comparisons are beyond 2 SE (max |z| 2.95, at 15 years). The points along each curve are correlated, so this is within noise. The curves overlap in every age group, children infected at birth included.
- The H side of the adapters is therefore right: `emod_age`, births with mothers (MTCT at birth), and HIV deaths requested through STIsim's people.

**EmodNetwork (`network_vs_home.csv`, `agegap_vs_home.csv`, `formation_vs_home.csv`, `warmup.csv`, `figures/network.png`)**
- Current relationships per person by type, sex and age group (1990, 2000, 2010, 2020): hybrid 7 / EMOD′ has a median |log ratio| of 0.05–0.09.
  - For scale, EMOD′ against itself with its mortality read as labelled (`fix_mortality`, which changes only the age structure) gives 0.04–0.07. Most of the residual is the demography, as expected.
  - Transitory, informal and marital relationships are within 0.78–1.27 in every cell.
- Man's age minus woman's at formation: transitory 3.27 vs 3.31 years, informal 4.56 vs 4.49, marital 3.58 vs 4.23, commercial 9.03 vs 9.42.
- Risk HIGH is 4.9–5.0% of 15–49 year olds in both (`epi_summary.csv`). So the starting population's sex-work history (`EmodRiskAt.start_history()`) puts the right share there from 1985.
- **The warm-up length matters.** With EMOD's own burn-in length (15 years), relationships were 17% from EMOD′'s (median), marital 1.45× in 1990, and the marital age gap 2.5 years instead of 4.2.
  - Cause: with ages frozen, long relationships pile up in the frozen age groups, and the entry rates set from that population stay wrong.
  - 5 years gives the best fit overall (2 years fits 1990 better but later years worse). It is the default (`warmup.csv`).

**The epidemic (not a fidelity check; `seeds.csv`, `figures/epidemic.png`)**
- Seeding 2% of Risk HIGH gives 0–5 seeds at 10,000 agents. The epidemic took off in 14 of 24 runs of hybrid 6 and 21 of 24 of hybrid 7. EMOD′ at 1% gets 7–23 seeds.
- In hybrid 6, HIVsim's FSW and clients are 2.0% of 15–49 year olds, against EMOD's 5% at HIGH. This is part of what the swap measures, not an error.
- 15–49 prevalence in 2020 is 4–5% in the hybrids against 11% in EMOD′. Neither hybrid has a cascade, and the spread across seeds is large.

**Mismatches bridged (and how)**

| # | Mismatch | Bridge (`models/emod/v4/adapters.py`) |
|---|---|---|
| 1 | EMOD′ modules read `sim.demographics.emoddemog` (births, mothers, 12 weeks) and `people.node`, `people.emod_age` | `EmodPeople`, an `ss.Demographics` named `emoddemog`, after `ss.Pregnancy` and `ss.Deaths` |
| 2 | Nodes: EMOD has 3 separate nodes; STIsim is national | `node` drawn from EMOD's InitialPopulation shares at the start; newborns take their mother's |
| 3 | `emod_age`: EMOD's age in days after the step's ageing (single precision); Starsim ages at the end of the step | `(age + dt)·365` in float64, 0 in the step of birth. A newborn is one step older than in EMOD from its second step on, which is negligible |
| 4 | Births with the mother: `ss.Pregnancy` (rc3) gives no per-step list in the step; callbacks run at the end of the step | Newborns are people aged 0 with a parent whose `ti_delivered` is this step (tested equal to the pregnancy's `birth` count). No prenatal network is given, so fetuses aren't agents and EmodHIV's transmission at birth is the only MTCT |
| 5 | 12 weeks pregnant (for the cascade later) | From `pregnancy.P` and `ti_P` (as `models/hivsim/cascade_host.py`) |
| 6 | Risk (step 6): EMOD's LOW/MEDIUM by node and HIGH while FSW/client; HIVsim's `risk_group` 0/1/2 and `fsw`/`client` | HIGH = current `fsw` or `client`; else LOW for `risk_group` 0, MEDIUM otherwise; set each step before the diseases. So EMOD's seeding (2% of HIGH, at the 1985 start since 1983 is before it) reads HIVsim's sex workers |
| 7 | STI co-infection flag (EMOD's R; HIVsim has none) | Drawn once at `StructuredSexual`'s debut (day 1 for those past it) with EMOD's probabilities by the Risk above |
| 8 | Acts: EMOD gives Poisson acts per step with a forced first act; STIsim gives whole acts per step, drawn once per partnership | STIsim's `acts` taken as they are (they are N's output) |
| 9 | Condoms | Same meaning (the chance an act has a condom); read as is |
| 10 | Relationship types | stable → MARITAL, casual → INFORMAL, one-time → TRANSITORY, sex work → COMMERCIAL; only the infection log reads them |
| 11 | EmodHIV read its network's edges directly | `EmodHIV.contacts()` hook (default: its network's edges; tested identical to before on EMOD′); `EmodHIVOnSti` overrides it |
| 12 | EmodNetwork's burn-in (1960.5–1975.5) sets entry rates from the population then; a 1985 start has no relationships | Warm-up of 5 years at the first step, with time and ages frozen, then rates frozen (see above) |
| 13 | EmodRisk gives sex work only to debuts from 1961, as they happen; a 1985 start has no sex-work history (and would treat everyone past debut as debuting in 1985) | `EmodRiskAt.start_history()`: for those past debut, the STI flag, and sex work drawn from their debut year (if from 1961), with uptake and dropout times, so they are HIGH, MEDIUM after it, or still to start |
| 14 | Seeding reads Risk HIGH at the first step's `step_state`, before networks and interventions | The warm-up and start history run from `EmodPeople` in the first step's demography, before the diseases |
| 15 | Numba's global random stream: EmodNetwork seeded it once at init, so any other Numba code drawing from it (or a second sim initialised in the process) moves its draws | `EmodNetworkAt` reseeds it from its own generator each step (`own_stream`). Tested: other Numba draws don't change a run with it, and do without it. No Starsim or STIsim code draws from Numba's global stream in these hybrids, so the risk is latent |

**Tests** (`models/emod/v4/test_adapters.py`, 6 tests, 14 s): both hybrids run 1985→2025 at 2,000 agents; births and nodes; repeatability; the Numba stream.

**Upstream friction (proposals, not worked around silently):**
- `stisim/sim.py:166-170, 278-349` (sti.Sim, rs-v4-art). A location's demographics (`demographics='malawi'`) can't be given together with another demographic module: either the location string or modules, not both. `adapters.sti_demography()` repeats its 10 lines (Pregnancy, Deaths and People from the data folder). Proposed: accept `demographics=['malawi', module, ...]` (the location's modules first), or a public `sti.location_demographics(location, datafolder, start, n_agents)` that `process_demographics` calls.
- `starsim/demographics.py` `Pregnancy` (rc3, ~line 1170, where `ti_delivered` is set). There is no in-step accessor for "this step's newborns and their mothers": callbacks (`diseases.py:846`) fire at the end of the step, after the diseases. Proposed: `preg.newborns` and `preg.delivered` (UIDs of the step), set at delivery.
- `stisim/networks/base.py:73`: `edge_type` is a float array (`ss_float`). Proposed: an int dtype, as `ss.SexualNet` types are categories.

**Open:**
- Commercial relationships in hybrid 7 drift above EMOD′'s from 2000 on (men 15–34 and women 15–34: 1.2–1.55×; commercial formations per 15–49 year old 1.20–1.23× from 2000 vs 0.97–1.11× in 1990), while the HIGH share is the same.
  - The frozen per-person entry rates (set from the 1985 warm-up) may not hold for commercial as the population's age structure shifts, or the starting FSWs' history may mix ages differently.
  - Check next: commercial queue lengths and FSW age at uptake against EMOD′'s.
- Seeding gives few seeds at 10,000 agents (1–3 out of 24 runs fail even in hybrid 7). For the swap path this is the decided scale, but Δs will need many seeds or conditioning on take-off.
- Expect to rebase on `swap-factories` (HIVsim Malawi on rc3): the hybrids here use a test setup (STIsim's Malawi demography and `StructuredSexual` with Malawi's condom data), not HIVsim's `make_sim`.
