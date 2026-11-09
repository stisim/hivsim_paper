# 12: EMOD′ components on HIVsim's slots — summary

**Answer:** both adapters work. Hybrid 6 (EmodHIV on `StructuredSexual` and STIsim's demography) and hybrid 7 (EmodHIV, EmodNetwork and EmodRisk on STIsim's demography) build and run 1985→2025. EmodHIV's natural history matches EMOD′'s on both demographies. EmodNetwork on STIsim's demography is within about 8% (median) of EMOD′'s network by type, sex and age, with the same age gaps. The exception was commercial relationships, 20–55% above EMOD′'s from 2000 on: an adapter error in the commercial entry rates, now fixed (addendum). The tables and figures below are from the rerun with the fix, except where marked as before it.

**EmodHIV (natural history; `survival_vs_home.csv`, `figures/survival.png`)**
- Kaplan–Meier survival from infection to HIV death, by age at infection (<15, 15–24, 25–34, 35–44, 45+) at 2, 5, 10 and 15 years. Pooled runs: EMOD′ (6 seeds, 314,000 infections), hybrid 6 (24 seeds, 30,000) and hybrid 7 (24 seeds, 35,000; 41,000 before the fix).
- 2 of the 40 comparisons are beyond 2 SE, both in hybrid 6 (max |z| 2.95, at 15 years; 4 before the fix, 2 of them in hybrid 7). The points along each curve are correlated, so this is within noise. The curves overlap in every age group, children infected at birth included.
- The H side of the adapters is therefore right: `emod_age`, births with mothers (MTCT at birth), and HIV deaths requested through STIsim's people.

**EmodNetwork (`network_vs_home.csv`, `agegap_vs_home.csv`, `formation_vs_home.csv`, `warmup.csv`, `figures/network.png`)**
- Current relationships per person by type, sex and age group (1990, 2000, 2010, 2020): hybrid 7 / EMOD′ has a median |log ratio| of 0.06–0.10 (0.05–0.09 before the fix).
  - For scale, EMOD′ against itself with its mortality read as labelled (`fix_mortality`, which changes only the age structure) gives 0.04–0.07. Most of the residual is the demography, as expected.
  - Transitory, informal and marital relationships are within 0.74–1.32 in every cell (0.81–1.27 before the fix; the other types differ only by the fix's new random draws).
- Man's age minus woman's at formation: transitory 3.49 vs 3.31 years, informal 4.54 vs 4.49, marital 3.49 vs 4.23, commercial 9.67 vs 9.42.
- Risk HIGH is 4.9–5.3% of 15–49 year olds in hybrid 7 and 4.9–5.2% in EMOD′ (`epi_summary.csv`). So the starting population's sex-work history (`EmodRiskAt.start_history()`) puts the right share there from 1985.
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
- (Resolved: see the addendum.) Commercial relationships in hybrid 7 drift above EMOD′'s from 2000 on (men 15–34 and women 15–34: 1.2–1.55×; commercial formations per 15–49 year old 1.20–1.23× from 2000 vs 0.97–1.11× in 1990), while the HIGH share is the same.
  - The frozen per-person entry rates (set from the 1985 warm-up) may not hold for commercial as the population's age structure shifts, or the starting FSWs' history may mix ages differently.
  - Check next: commercial queue lengths and FSW age at uptake against EMOD′'s.
- Seeding gives few seeds at 10,000 agents (1–3 out of 24 runs fail even in hybrid 7). For the swap path this is the decided scale, but Δs will need many seeds or conditioning on take-off.
- Expect to rebase on `swap-factories` (HIVsim Malawi on rc3): the hybrids here use a test setup (STIsim's Malawi demography and `StructuredSexual` with Malawi's condom data), not HIVsim's `make_sim`.

## Addendum (9 Nov): the commercial drift was the entry rates, set from 1985's clients per FSW

**Answer:** an adapter error, now fixed. EMOD′ freezes each type's queue-entry rates per person at the end of its burn-in (1975.5). Commercial has one age bin and all HIGH hold its flag, so its rates depend only on the available clients per FSW in each node: an FSW's rate is 0.075·(1 + clients/FSW) per day. In 1975.5 sex work had run only since 1961's debuts, and clients (who start up to 10 years after debut and stay 2–30 years) were still accumulating: 2.8 HIGH men per HIGH woman, against 3.4 in 1985 and 4.1–4.4 from 1995. The warm-up set the rates from 1985's population instead, so FSW got 1.25× EMOD′'s rate. Clients outnumber FSW about 4 to 1 and their queue entries pile up unmatched (EMOD′ in 2020: 0.4 million men's commercial queue entries, from 2% of its agents, and no women's), so the FSW side sets the number of commercial relationships, and the higher FSW rate went straight into them. The HIGH share was right all along; the rate per FSW was not. (`commercial.py`, `outputs/commercial*.csv`, `figures/commercial.png`; EMOD′ and hybrid 7 without HIV, 6 seeds.)

**Evidence (before the fix, `outputs/commercial_before.csv`, `commercial_rates_before.csv`):**

| | EMOD′ | 7 (before) | 7 (after) |
|---|---|---|---|
| Rates frozen at | 1975.5 | 1985 (warm-up) | 1985, commercial from the Risk at 1975.5 |
| Available clients per FSW at the freeze (rate_f / rate_m), nodes 1–3 | 1.86, 2.26, 3.41 | 3.98, 3.17, 4.71 | 1.51, 3.04, 3.35 |
| Commercial formed per FSW per year, 1985–2024 | 102.6–103.9 | 122–132 | 99–108 |
| Commercial formed per 15–49 year old, 2000 / 2010 / 2020 (ratio to EMOD′) | 1.06 / 1.00 / 0.99 | 1.14 / 1.23 / 1.24 | 0.96 / 0.95 / 0.98 |

- Commercial formations per FSW are flat in both (the rate is frozen and FSW are the limiting side), so the 1.25× is constant from 1985. It looked like a drift from 2000 because in 1990 the hybrid has 24% fewer FSW per 15–49 year old than EMOD′ (0.81% vs 1.06%, below), which hid it then (formations 0.97× EMOD′'s in 1990).
- The mean age of those at HIGH is the same in both (women 21.1–21.6 years, men 27–29), so the starting FSWs' ages (`start_history()`) were not the cause.

**Fix (`models/emod/v4/adapters.py`, `network.py`):** `EmodRiskAt.start_history(then=1975.5)` also returns the Risk each person's drawn history gives them at the end of EMOD's burn-in. After the warm-up, `EmodNetworkAt` (`burnin_risk=True`, the default) sets that Risk, redraws flags and maxima for it, sets the commercial rates alone from those available (`EmodNetwork.set_rates(eligible(), types=[COMMERCIAL])`), and puts the start's Risk back (flags and maxima redrawn again). `set_rates()` and `eligible()` are the burn-in block of `EmodNetwork.form()` moved into methods, with nothing else changed.
- EMOD′ is bit-identical: its rows of `net.csv`, `form.csv`, `epi.csv`, `gap.csv` and `surv.csv` (all three EMOD′ models, 6 seeds each) are the same as from the parent commit. Hybrid 6 is unchanged too. Tests pass (6).
- **Tried first and dropped:** running the whole warm-up with the Risk at 1975.5 and setting all types' rates from it. Commercial came right, but young marital and informal relationships rose (marital, women 15–24, 1.58× in 1990; men 25–34 1.20–1.27× in 2000–2020; non-commercial cells 0.70–1.58× EMOD′'s, from 0.81–1.27). The other types' rates depend on the age structure as well as the Risk, and 1975's sex work on 1985's ages fits worse than 1985's on both. Only commercial's rates depend on nothing but the sex work's history.

**Rerun of the network comparison with the fix (`network_vs_home.csv`, `formation_vs_home.csv`, `agegap_vs_home.csv`, `figures/network.png`, `figures/commercial.png`):**
- Commercial relationships per person (men 15–24, 25–34, 35–49; women 15–24, 25–34) in 2000, 2010 and 2020: 0.85–1.09× EMOD′'s (|z| ≤ 1.1) except 1.28× for men 15–24 and 1.36× for women 25–34 in 2020 (|z| 1.8–1.9); before the fix, 1.03–1.56×, |z| up to 6.7.
- Commercial formations per 15–49 year old: 0.96, 0.95, 0.98× EMOD′'s in 2000, 2010, 2020 (before: 1.14, 1.23, 1.24×).
- Median |log ratio| over all cells: 0.10, 0.07, 0.07, 0.06 in 1990, 2000, 2010, 2020 (before: 0.09, 0.05, 0.08, 0.08). Without commercial, 0.09, 0.08, 0.07, 0.06 (before: 0.08, 0.05, 0.07, 0.06): the other types move only by the new random draws (median |z| 1.1–1.8 both times). Non-commercial cells are 0.74–1.32× EMOD′'s (before 0.81–1.27).
- Age gaps: commercial 9.67 years vs EMOD′'s 9.42 (before 9.03); transitory 3.49 vs 3.31, informal 4.54 vs 4.49, marital 3.49 vs 4.23 (as before).
- Risk HIGH 4.9–5.3% of 15–49 year olds (EMOD′ 4.9–5.2%). `warmup.csv` (the choice of warm-up length) is from before the fix and wasn't rerun; the fix only touches commercial's rates, set after the warm-up. Survival in hybrid 7: no comparison beyond 2 SE (max |z| 1.85).

**Still open: fewer FSW from 1986 to 1997.** In hybrid 7, FSW (Risk HIGH women) are 0.98% of 15–49 year olds in 1985 (EMOD′ 1.08%). The share then falls to 0.86% in 1992 and is back to EMOD′'s by 1998–2000 (30 further seeds: 0.98, 0.95, 0.88, 0.86, 0.91, 0.98, 1.00% in 1985, 1986, 1990, 1992, 1995, 1998, 2000; SE 0.02; EMOD′ 1.08, 1.10, 1.07, 1.05, 1.02, 1.02, 1.01%). It isn't specific to commercial relationships and was there before the fix, but with the fix it shows: commercial formations per 15–49 year old are 0.76× EMOD′'s in 1990 and 0.88× in 1995.
- What's known: new FSW per year (8/yr at 10,000 agents, 1985–1990) are what women's debuts give. The starting FSWs' remaining time as FSW (mean 3.0 years) is as expected for EMOD's Weibull. The starting women still to take up sex work are 0.42 per FSW, against EMOD′'s 0.58 in 1985. The difference may be only demographic (STIsim's population grows 2.7% a year from 1985, EMOD′'s 4.6%; 47 women's debuts a year per 1000 women 15–49, against EMOD′'s 54), but the cause isn't pinned down.
- Check next, if it matters for the swap: the starting population's debut years against those that 1985 onwards implies, and the stock that `start_history()` gives against a run of EmodRisk alone on STIsim's demography from 1961.

