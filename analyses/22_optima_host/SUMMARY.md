# 22: Optima's population label, and `CascadeHost` on Optima′'s agents — summary

**Answer:** both pieces work. The **population label** (`models/optima/v4/poplabel.py`) gives each EMOD′ or HIVsim agent an Optima population from age, sex and current sex work, updated every step. On HIVsim's people it reproduces Optima's age and sex structure to about a point. It can't reproduce Optima's key populations: MSM are absent from both models, HIVsim has half Optima's FSW, and EMOD′ has 1.3 times its clients. The ways people enter and leave them also differ (list below). **`OptimaHost`** (`models/optima/v4/host.py`) implements `CascadeHost` on Optima′'s declared states as agents, including `cascade_stage()` with Figure 3A's categories (D8). With it, `common/outputs.py`'s `Outputs` runs on Optima′ as agents, and adding it doesn't change the run: all 16 seeds are identical, year by year, to analysis 20's. The source shares agree with the Optima′ ODE's attribution (analysis 19's method): mean |agents − ODE| is 0.35 points over 2000–2020, with 6% of state-years at |z| > 2. In 2020 the D8 split is 43.3 / 39.4 / 12.2 / 5.1 (undiagnosed / diagnosed never treated / on ART / interrupted), against 42.1 / 39.3 / 12.7 / 5.8 for the ODE and 42.3 / 40.7 / 11.4 / 5.7 for Optima's tracer. The logged sexual infections times the scale equal the disease's own results exactly.

## 1. The population label

**Definitions** (`poplabel.definitions()`, checked against the project's `popkeys`): FSW (women), Clients and MSM (men), all 15–49; and F/M by 5-year band, 0–4 to 85+. `label(age, female, fsw, client, msm)` is the definition read the other way: a woman aged 15–49 who sells sex now is FSW, a man aged 15–49 who buys it is a client, MSM likewise, and everyone else goes into their sex and band. `stisim_label(sim)` uses STIsim's `StructuredSexual.fsw` and `.client` (currently selling or buying). `emod_label(sim)` uses `emod_age` and Risk HIGH (a current FSW or client, by sex). `PopLabel` is an analyzer that keeps `optima_pop` up to date each step and records the sizes at each whole year.

**Sizes, 1990–2020** (`outputs/sizes.txt`, `key_pops.csv`, `age_structure.txt`). Optima′'s own sizes agree with Optima's: the ODE follows its forced `popsize` exactly, and the agents (analysis 20, 48 seeds at 100k) are within 0.3% for every population. EMOD′: parameter set 0 at 1% from 1960.5; HIVsim: provisional reference set at 100,000 agents from 1985.

| | Optima | EMOD′ | HIVsim |
|---|---|---|---|
| FSW, % of women 15–49 | 2.11 (constant: forced) | 1.9–2.2 | 0.9–1.1 |
| Clients, % of men 15–49 | 6.36 | 7.7–8.3 | 6.1–6.7 |
| MSM, % of men 15–49 | 0.99 | 0 | 0 |
| Total population, 2020 / Optima | 1 | 1.96 (1.41 in 1990) | 1.01 |
| Women 15–49 (general), % of women, 2020 | 47.8 | 44.7 | 47.5 |
| Women 75+, % of women, 2020 | 1.5 | 5.3 | 1.2 |
| Largest band difference, % of sex, 2020 | | 3.7 (M85+) | 1.1 |

**Mismatches: where the label can't match Optima**
1. **No MSM** in either model: Optima's MSM population (1% of men 15–49, 45k in 2020) stays empty.
2. **FSW size:** HIVsim has about half Optima's FSW (1.0% against 2.1%). EMOD′ matches to within 0.2 points, by coincidence: neither model is fitted to Optima's sizes.
3. **Client size:** EMOD′ has 1.3 times Optima's share (~8% against 6.4%); HIVsim matches.
4. **Key populations have no ages in Optima.** Nobody ages into or out of FSW, Clients or MSM; leaving FSW (`risktransit`: 32 years to each of the 7 bands, ~4.6 years in sex work) puts the person in any of the seven 15–49 bands alike, Clients likewise (47 years each, ~6.7 years). The label leaves sex work in the person's own band, and everyone leaves a key population at 50.
5. **Entry, exit and duration are each model's own.** Optima: balanced flows, fixed (forced) shares. EMOD′: lifetime fate at debut by node, entry after a uniform delay (≤10 years women, ≤5 men), exit after a Weibull time for women (~5.4 years) or uniform 2–30 years for men (~16). HIVsim: lifetime share (5.5% of women, 23% of men), starting ages normal(20, 3) and normal(25, 5), durations lognormal(5, 3) and lognormal(10, 5) years. Different turnover means different HIV histories in the key populations even at the same size.
6. **Risk within the general population is lost.** Optima has none, so EMOD′'s LOW/MEDIUM and HIVsim's risk groups 0/1/2 have nowhere to go. In the other direction, `OptimaHost.risk()` gives HIGH to FSW and Clients and MEDIUM to everyone else, with no LOW; an adapter would need a parameter for EMOD's LOW fraction or HIVsim's shares (swap_design §1.4).
7. **EMOD′'s demography itself:** 1.4× Optima's size in 1990 and 2× by 2020, 75+ 5.3% of each sex against 1.5%: EMOD′'s copied mortality bug (fixed in `emod_fixed`, analysis 17), not the label. HIVsim's age structure is within about a point of Optima's in every band.
8. **EMOD′'s nodes** (3, no mixing) have no Optima counterpart; the label ignores them.
9. **When the label updates:** `PopLabel` is an analyzer (end of step). An Optima′ component running on these people in a hybrid should call `stisim_label`/`emod_label` itself, or the label should move into a module that runs earlier.

## 2. `OptimaHost` and `Outputs` on Optima′ as agents

**Mapping** (`host.py` docstring): sex from the population; age a fixed value within the band (key populations 15–49); CD4 a fixed value within its band, WHO stage from CD4, first symptoms on first reaching CD4 < 200. `test` moves undx → dx; `start_art` moves any untreated state with HIV to svl or usvl, by arm; `stop_art` moves to lost; `circumcise` moves susreg → progcirc (all via `ss.move(src, dst, uids=)` per CD4 band). Debuts are ageing from 10–14 into 15–19. Pregnancy, the 12-week and 6-week events, PMTCT and PrEP can't be per agent (Optima's births have no mothers; its PrEP is a coverage): the host gives no events for these and raises if asked to act. `OptimaTracker` (after Optima's rules) keeps an ever-treated flag, first symptoms and the step's events (from `Disease.add_callback` on infection, birth and HIV-death transitions).

**Sources:** Optima's force of infection has no pairs. For each infection the host draws the source cell (partner population × infected state) from Optima's own parts of the force of infection on the recipient's population and circumcision status in that step (`Transmission.parts`, kept by a one-line addition to `optima_prime.py`), then an agent in that cell, with its own random numbers. `Outputs` takes transmissions, new infections and HIV deaths from a host that provides them, and `TransmissionLog` reads the source's sex from the host (7 lines changed).

**Check** (`run.py outputs`, 16 seeds at 100,000 agents; `run.py compare`, `outputs/compare.txt`), sexual transmissions summed over each year, 2000–2020:
- **Doesn't change the run:** in all 16 seeds, sexual infections are identical every year to analysis 20's runs with the same seed; analysis 20's seed 1 rerun with this code is identical in every output.
- **Totals:** logged sources × scale equal the disease's results exactly. Mean +0.13% against the ODE; largest gap 2020 (+8.4%; these 16 seeds sit higher than analysis 20's 48, at +5.8%): analysis 20's known 2020 gap from Optima's target-driven ART at 100k agents.
- **By source state:** against the ODE, mean |difference| 0.35 points, largest 3.0; 6% of state-years at |z| > 2, mostly 2015–2020 (2020: dx 35.8% vs 38.8%, z −4.5; care 2.9% vs 0.5%, z +2.2): analysis 20's random walk in the care stock under the ART target. The sources follow the states.
- **Figure 3A (D8):** the host's ever-treated flag and the transmission log's own ART history give the same shares (within 0.03 points). 2020: agents 43.3 / 39.4 / 12.2 / 5.1; ODE 42.1 / 39.3 / 12.7 / 5.8; Optima's tracer 42.3 / 40.7 / 11.4 / 5.7.
- **One systematic difference:** "interrupted" in 2000 is 1.0 ± 0.1% as agents against 0.2% in Optima's tracer. Likely (not checked) the ART target's churn as agents: random overshoots move people from ART back to care, and once flagged they stay ever treated. From 2010 on within noise (1.4 vs 1.1, 3.1 vs 3.0, 5.1 vs 5.7).

**Unchanged:** Optima′'s ODE acceptance identical to analysis 20's (`logs/acceptance_after.log`). Existing tests pass (25); new `tests/test_optima_host.py` (3 tests). ~190 s per seed at 100k; `Outputs` adds ~10%.

## 3. What's left: Optima's cascade as a table

- **Rate pathways** in `common/cascade`: a per-step probability from a rate by group (sex, age, CD4, Optima population) and year, for people in a given place in care: Optima's `hivtest`, `aidstest`, `linktocare`/`aidslinktocare`, `leavecare`/`aidsleavecare`, `returntocare` (and HIVsim's two testing rates). On agents it needs the population label as a group key.
- **Stock-target pathways:** number on ART by year (`numtx`), filled by a priority rule (lowest CD4 first before 2018) and removed when over; and shares held from a year (`fixpropdx`, `fixpropcare`, `fixproptx`, `fixpropsupp`). Also HIVsim's `sti.ART(coverage=, vls_coverage=)`. `agents.CascadeTargets` has the arithmetic; port it to act through `CascadeHost`.
- **Rows for Optima's table**, extracted from the project (as `extract_cascade.py` does for EMOD), including treatment failure (`treatfail`) and viral-load regain (`regainvs`, `numvlmon`).
- **Interface gaps:** `CascadeHost` has no "into care" (pre-ART) or "lost from care" move; Optima's PMTCT and PrEP are coverages; no mothers for 12-week or infant events. Each needs a decision (per-person status drawn from coverage, or the pathway as a rate).
- **For `rates.py`:** both new pathway types also need a compartmental form, so the Optima′ ODE can run the same table.

## Starsim friction (rc3), with proposed upstream fixes

1. **`people.female` (and `people.age`) disagree with a declared model's own groups.** Optima′'s `pop` encodes sex, but `people.female` is random and agrees for only 50% of agents; anything reading it is silently wrong. Proposed: let `ss.Disease.stratify()` declare a group's sex (as `groups.set_sex`, `internals/dimensions.py:1231`, does for a `sex` dimension). Workaround: the host takes sex from `pop`.
2. **No per-agent "which state" for a declared disease** (the host loops over 36–38 `BoolState`s). Proposed: `disease.state_of(uids)`.
3. **`ss.move(uids=)` moves from one source state** (`transitions.py:973`; 6 calls per move across CD4 bands). Proposed: `ss.move({src: dst, ...}, uids=)`.
4. **Callbacks work well:** `Disease.add_callback` on declared transitions gave exact per-agent events with no change to the run; worth documenting for declared models as agents.

**Files:** `README.md`, `run.py` (`label`, `sizes`, `outputs`, `ode`, `compare`), `age_structure.py`; `outputs/` (sizes, key_pops, age_structure, ode_by_state, ode_by_source_steps, per-seed outputs, compare.txt); `logs/`. New: `models/optima/v4/poplabel.py`, `host.py`, `tests/test_optima_host.py`. Changed: `models/optima/v4/agents.py` (`make_sim(interventions=)`), `optima_prime.py` (`Transmission.parts`, output only), `common/outputs.py`, `common/transmission.py`.
