# Component swaps: design (P5)

*Draft, 4 Nov 2026; status of the interfaces and paths updated 10 Nov 2026.*

* This is the design for the paper's core analysis: swapping components between HIVsim, EMOD′ and Optima′ one at a time, and attributing the gap in policy outputs to them. It also feeds P4 (the component map and interfaces, sprint 4), whose boundaries need approval by Robyn, Daniel and Cliff.*

Throughout, **Decided** cites the doc and date where it was decided. Everything else is **proposed** here and needs a decision. The open decisions are listed at the end.

## Summary

- Each model is cut into the same **slots**: demography, risk and sex work, sexual contact, transmission, natural history (with the effects of ART), care cascade (with VMMC and PrEP), seeding, and measurement. Optima has one more, the solver (ODE or agents).
- Some slots are **fused** in code today: EMOD′ and HIVsim each have transmission, natural history and seeding in one disease module. In Optima, demography and risk are one thing (its 39 populations), and so are contact and transmission (its mean-field force of infection). The first swaps use these coarse slots; finer ones need the refactor in sprint 6.
- **The cascade is the one slot with a shared interface already** (`CascadeHost`, `common/cascade/`). This design makes the condensed cascade the common currency for the cascade: each model's cascade becomes a parameter table, and two new pathway types (rates, and stock targets) let HIVsim's and Optima's cascades fit it.
- **First path: HIVsim → EMOD′** (proposed in place of HIVsim ↔ Optima′ in `plan.md`, sprint 7), with four coarse slots in the order cascade → HIV disease → network and risk → demography. Beyond what analysis 07 is building, it needs two adapters, the outputs analyzer and the runner.
- **Attribution:** with four coarse slots, all 2⁴ = 16 hybrids cost about 16 core-hours, so exact Shapley values (the average over all 24 orders) are affordable. Sequential paths are what the figures show; Shapley is the summary number.
- **Calibration:** the headline is compensated (tier 2): after each swap, refit two knobs, transmission scale and ART uptake, to the primary targets. Uncalibrated (tier 1) results are reported alongside. **Endpoint identity** is checked bit for bit; it holds by construction if each model's own `make_sim` is built from the same slot factories as the runner.

## 1. Components

### 1.1 Slots

Each slot has three layers, swapped separately: mechanism, inputs and parameters. By default all input swaps are one step (**Decided**, `prompts.md`, 27 Oct).

| Slot | What it does | Reads | Writes |
|---|---|---|---|
| **D** Demography | Initial population and start year, births and pregnancy, background deaths, ageing; Optima's population forcing; EMOD's nodes | Its inputs | People: alive, age, sex, node; events: pregnancy at 12 weeks, birth (with the mother), 6 weeks old |
| **R** Risk and sex work | Risk group, who is a sex worker or client and for how long; EMOD's STI flag; Optima's risk transitions between populations | People | `risk`, `fsw`, `client`, `sti`, and the Optima population label (1.3) |
| **N** Sexual contact | Partnerships, their types, coital acts and condom use | People, R | Edges: `p1` (man), `p2` (woman), type, acts this step, share of acts with a condom |
| **T** Transmission | Per-act probability by direction and age; cofactors (condoms, circumcision, STI, PrEP); how acts combine; mother-to-child | N's edges, H's infectiousness, susceptibility from C and R | New infections, each with its source and route |
| **H** Natural history and ART biology | Stages, CD4, infectiousness by stage, HIV deaths; what ART does (mortality on ART, CD4 recovery, infectiousness when suppressed or not) | Infections from T; ART starts and stops from C | CD4, stage, WHO stage, symptom onset, relative infectiousness, HIV deaths; the HIV side of `CascadeHost` |
| **S** Seeding | How the epidemic starts: EMOD infects 2% of HIGH-risk people in 1983; HIVsim sets prevalence by risk group in 1985; Optima sets each population's state in 1990 | People, R, H | Initial infections |
| **C** Care cascade and prevention programmes | Testing, staging, linkage, pre-ART and care, ART start and arm, failure and regain, dropout and return; PMTCT; VMMC and PrEP numbers | Through `CascadeHost` only | Through `CascadeHost` only |
| **M** Measurement | The policy outputs (2.2), and how transmissions are attributed to cascade stages | T's infection log; H; C's place in care | Outputs |
| **X** Solver (Optima only) | Optima′ as an ODE or as agents | | |

**Decided:** the ODE-to-ABM switch is its own component, swapped first or last on every Optima path, and all hybrids run as agents (`prompts.md`, 27 Oct).

### 1.2 What each model has in each slot

| Slot | EMOD′ (`models/emod/v4/`) | HIVsim (`models/hivsim/`, STIsim; `sti.ART` its own module) | Optima′ (`models/optima/v4/optima_prime.py`) |
|---|---|---|---|
| D | `EmodDemography`: from 1960.5, 3 nodes, 1% scale, mortality read 1000× too low (copied) | STIsim demographics (`malawi_*` files), from 1985, 10,000 agents | Births, background deaths, ageing between 5-year bands, risk transfers; `ForcePopsize` (writes the engine's compartments); from 1990 |
| R | `EmodRisk`: Risk LOW/MEDIUM by node, HIGH while FSW or client; STI flag | Inside `StructuredSexual`: `risk_group` 0/1/2 by sex; `fsw`, `client` (lifetime fate, durations) | The 39 populations themselves, with `risktransit` |
| N | `EmodNetwork`: pair formation with queues, 4 types, frozen rates after 1975.5 | `StructuredSexual` (stable, casual, one-time, sex work) and `MaternalNet` | None: acts per year between pairs of populations (regular, casual, commercial) |
| T | Inside `EmodHIV`: Poisson acts, forced first act, multipliers on probabilities | Inside `sti.HIV`: `beta_m2f`, multipliers on `rel_trans`/`rel_sus` | `Transmission`: force of infection per population from compartment counts, with `force` and `inhomo` |
| H | Inside `EmodHIV`: Weibull prognosis by age, CD4 deterministic from it, acute/latent/AIDS; ART mortality tables | `sti.HIV` (acute, latent, falling CD4, deaths by CD4) and `sti.ART`'s biology | 6 CD4 bands, death by band, CD4 moves on ART |
| S | `EmodHIV` (1983, HIGH risk, fast-forwarded) | `sti.HIV` `init_prev_data` (1985) | Initial shares by population and state (1990) |
| C | Condensed cascade with EMOD's table (`common/cascade/`, analysis 06; the literal and compiled campaigns are its references) | `sti.HIVTest` ×2, `sti.ANCTest`, `sti.ART(coverage=n_art, vls_coverage=vls)`, `sti.VMMC` | Cascade transitions and `CascadeTargets` (number on ART, holds from 2020; writes the engine's compartments) |
| M | Infections logged with source and stage (`EmodHIV.infections`) | `common/outputs.py` on `StisimHost` (tracked source) | ODE: compartment flows; agents: `Outputs(OptimaHost)`, source drawn by each stage's share of the force (analysis 22) |

### 1.3 Interfaces

| Interface | Between | Exists? | Notes |
|---|---|---|---|
| `CascadeHost` | C ↔ H, D, R | **Yes** for all three: EMOD′ (`cascade_host.py`), STIsim (`StisimHost`, analysis 07), Optima′ as agents (`OptimaHost`, analysis 22) | Per agent. For Optima′ agents it moves people between declared states with `ss.move()` (Starsim #1492) |
| Cascade table | C's parameters | **Yes** for EMOD (`common/cascade/params.py`, analysis 06) and Optima (`tables/malawi_optima_params.csv`, analysis 23, with the rate and target/hold pathway types). HIVsim's: **to write** | Needed for EMOD′ → HIVsim |
| Cascade as compartment rates | C → Optima′ ODE | **Yes** (`common/cascade/rates.py`, analysis 08) | Used for the ODE side of Optima paths (2.3) |
| Contact edges | N → T | De facto: both agent networks are `ss.SexualNet`s with `p1`, `p2`, `acts`, `condoms` | Field meanings differ (acts per step or per year; condom probability or share). **To write:** one spec and an adapter on each side |
| Infectiousness and susceptibility | H → T, C/R → T | Starsim's `rel_trans` and `rel_sus` exist; EmodHIV computes its own (`infectiousness()`) | **To write:** EmodHIV to expose relative infectiousness per agent in the same form |
| Life events | D → H, C, N | Partly inside `CascadeHost` (debut, 12 weeks pregnant, 6 weeks old) | **To write** as a D interface: pregnancy, birth with mother, debut age, `age_years` |
| Population label | R → T, C (Optima) | **Yes** (`models/optima/v4/poplabel.py`, analysis 22; no MSM in EMOD′ or HIVsim) | Every R maps its people to Optima's 39 populations (sex, 5-year band, FSW, client); Optima′'s R provides `risk`, `fsw`, `client` back |
| Seeding rule | S | **To write** | Each S seeds at its own year, or at the sim start if later (1.4) |
| Outputs | M | **Yes** (`common/outputs.py`, analyses 15–16; `cascade_stage()` on every host) | Needs T to log each infection's source and H/C to give each person's cascade stage (`CascadeHost.cascade_stage(uids)`, a new method) |

Starsim doesn't enforce interfaces, so these are our own layer, enforced by our own tests (**Decided**, `prompts.md`, Starsim v4 notes).

**Two new pathway types for the condensed cascade (proposed).** Without them only EMOD's cascade fits the table.
- **Rate pathways:** a per-step probability from a rate by group (sex, age, CD4, Optima population) and year, for people in a given place in care. This covers HIVsim's general and low-CD4 testing, and Optima's `hivtest`, `aidstest`, `linktocare`, `leavecare` and `returntocare`.
- **Stock targets:** a number (or share) on ART by year, filled each step from those eligible by a priority rule, and removed when over; and shares held from a year. This covers HIVsim's `sti.ART(coverage=..., vls_coverage=...)`, and Optima's `numtx` (lowest CD4 first before 2018), `fixproptx` and `fixpropsupp`.

With these, all three cascades run on one mechanism, and swapping the cascade's mechanism (event-driven versus target-driven) becomes a parameter-layer step that can be measured on its own.

### 1.4 Mismatches

| Mismatch | EMOD′ | HIVsim | Optima′ | Proposed bridge |
|---|---|---|---|---|
| **Population groups** | Risk LOW/MEDIUM by node; HIGH = current FSW or client; 3 nodes with no mixing | `risk_group` 0/1/2 by sex; `fsw`/`client`; national | 39 populations: FSW, Clients, MSM (15–49) and sex × 5-year band; no general-population risk | Population label (1.3). Agents from EMOD or HIVsim get Optima's label from sex, age and sex work, with no MSM. Optima′ agents get `risk` from EMOD's LOW fraction, or HIVsim's shares, as a parameter of the adapter. On Optima paths D and R move as one step |
| **Relationship types** | Transitory, informal, marital, commercial | Stable, casual, one-time, sex work | Regular, casual, commercial acts per population pair | No mapping is needed for the agent networks: each T reads acts and condoms per edge, and type only for its own per-type parameters (an adapter maps types) |
| **CD4 against viral load** | CD4 deterministic from the prognosis; stage by time to death; suppressive or non-suppressive arm fixed at ART start | CD4 continuous; acute/latent/falling; effective or non-suppressive ART, share set by `vls_coverage` | 6 CD4 bands; suppressed or unsuppressed on ART; 63% suppressed at start (the step artefact, kept) | None of the three tracks viral load. The interface is CD4 (continuous; Optima′ gives a value inside the band, fixed per person) and suppressed or not on ART. EMOD's 250 threshold (2008–2011) falls inside Optima's 200–350 band, so it is decided within the band by that value |
| **WHO stage and symptoms** | WHO stage; `NewlySymptomatic` event | Neither (07 derives them) | AIDS = CD4 < 200; `aidstest` is a rate | Derived from CD4 where missing. This is part of H, so it moves with H |
| **Cascade drive** | Flows: probabilities and delays; no target | Stock: ART corrected to a number each step (re-engagement only through that) | Stock: `numtx`, then shares held from 2020 | Rate and stock-target pathways (1.3). This is expected to be a large part of the gap in transmissions from people who interrupted treatment |
| **Interrupted treatment** | `ARTStatus` OffART (the dropout trap, copied) | `art_discontinued` | `lost` mixes people lost before and after ART | `cascade_stage()` has an "interrupted" stage; Optima′ needs `lost` split into two states (same rates) to report it. It's a change to measurement, not to dynamics |
| **Infectiousness after diagnosis** | None | None | `effdx` lowers it for everyone diagnosed | Part of Optima's T, read from `cascade_stage()` |
| **Transmission arithmetic** | Multiplies probabilities per act | Multiplies hazards | Whole and fractional acts per step of 0.2 years | Kept: each T keeps its own arithmetic |
| **Prevention inputs** | VMMC and PrEP as numbers (NChooser, truncated at 1% scale) | VMMC coverage by age | `propcirc`, `prep` as coverage in the force of infection | In agents, coverage becomes per-person status drawn by population each year (a flow row in C's table) |
| **Start year and initialization** | 1960.5, network burn-in to 1975.5, seeding 1983 | 1985, prevalence set at the start | 1990, state set by population | **Decided** (Robyn, 9 Nov): all hybrids start in 1985, because 1960.5 is too slow. EMOD′ runs a burn-in from 1960.5 to 1985 once per parameter set and seed, saves the sim, and every hybrid with EMOD's slots continues from it. The 1985 state is part of D (and S), so it swaps with them. A resumed EMOD′ must equal a continuous one exactly (a test). Slots swapped in at 1985 that carry state (network, natural history) start from the snapshot's people; the network may need a short warm-up with transmission off. Optima's native 1990 state: to settle with the Optima paths |
| **Scale** | 1% (36,600 agents in 1960); NChooser truncation depends on it | 10,000 agents in 1985 | Counts | **Decided** (Robyn, 9 Nov): agent numbers vary and move with D: hybrids with HIVsim's D run at 10,000 agents, with EMOD's at 1%. Checked once: HIVsim at 10,000 agents against HIVsim at 1%, in distribution. EMOD's absolute numbers (VMMC, PrEP) are scaled to agents with random rounding when its C runs on HIVsim's D |
| **Time step** | Monthly (30.42 days) | Monthly | ODE at 0.02 years, from rates per 0.2-year step | Agent hybrids are monthly. Optima′ as agents uses the same rates; the sprint 5 dt study checks it |
| **Compartments against agents** | Agents | Agents | Compartments; `ForcePopsize`, `CascadeTargets`, `Births` and `Transmission` read or write `engine.y` | Agent versions of these four: **done** (`models/optima/v4/agents.py`, analysis 20). They are also what any Optima slot needs to run in a hybrid |

## 2. Swap paths

### 2.1 Steps, orders and what they visit

A path from A to B in a given order swaps one slot (or one layer of a slot) per step. Each step is a hybrid: B's slots so far, and A's slots for the rest.

Running the same order backwards, from B to A, visits **the same hybrids**. So it adds nothing. Order effects come from running **the same order from the other end** (B to A, swapping c₁ first), or from Shapley values (2.5). `prompts.md` (open swap design) proposes each single swap both ways, one full path, and "the same path in reverse". Proposed reading: the reverse is that same order, started from B.

**Order, proposed: cascade → HIV disease → network and risk → demography.** Seeding moves with the HIV disease while the two are one module.
- The cascade first, because its interface already exists in both directions (EmodHost, and STIsim's host from 07), and its effect from the first position is the headline result: how much of the gap is the cascade.
- Demography last, because it carries the start year, the scale and EMOD's mortality bug (twice the census by 2008). Swapped early, it would colour every later step. Swapped last, its effect on shares (not counts) should be small, and that is testable.
- Each step needs one adapter: B's new slot running on A's remaining slots (2.3).
- The decided starting point in `prompts.md` mentions an "outside-in" order. It is read here as programme → biology → behaviour → population. If Robyn meant the reverse, the lattice in 2.5 covers both.

### 2.2 Outputs decomposed

From `vision.md`, with years proposed here:

| Output | Definition | Year |
|---|---|---|
| Transmissions by the source's cascade stage | Share of new infections whose source was undiagnosed, diagnosed and never on ART, on ART and suppressed, on ART and not suppressed, or interrupted (off ART after starting) | 2024 (as MIHPSA), and 2020 |
| Deaths on ART | Share of HIV deaths among people on ART | 2024, and 2020 |
| HIV incidence | Per 1000 person-years, 15–49, by sex; and new infections, all ages | 2024, and 2020 |
| Intervention ranking | Better testing against better retention, under the treatment-interruption scenario from 2026 (**Decided**, `prompts.md`, 27 Oct) | 2026–2035; January (`plan.md`) |

Along with these, every hybrid also records the common targets (`common/targetfit.py`), to show what the data pin down and what they don't.

Optima′ is checked only for 2000–2020 (**Decided**, 30 Oct), and its projections are ill-conditioned (`optima_inventory.md`). Endpoint identity doesn't depend on this, but Optima's 2024 values are as much its post-2020 rules (shares held from 2020) as its fit. HIVsim's ART input ends in 2021, so it needs a projection rule to 2025 (proposed: hold coverage at 2021).

**Measurement** (**Decided**, `prompts.md`, 27 Oct): each endpoint keeps its model's own attribution method, and the switch of method is one explicit step. Proposed:
- On every agent hybrid, record **both** cheap methods: tracked source, and cascade-weighted (each stage's share of infectiousness-weighted exposure). This is option (c) of the recorded options, at no extra run cost.
- Where T is Optima's mean-field force of infection, there is no individual source. Each infection draws its source's stage in proportion to the stages' contributions to its force, which equals the cascade-weighted share in expectation. Tracked source is then defined on every hybrid.
- On Optima paths the measurement step sits next to the ODE ↔ agents step. On HIVsim ↔ EMOD′ both models use tracked source, so there is no measurement step.

### 2.3 The three paths

Coarse slots, as they can be swapped with the code as it is: **D**, **RN** (risk, sex work and network), **V** (the HIV disease: transmission, natural history, ART biology, mother-to-child and seeding), **C** (cascade and prevention programmes). They are split further (V into T, H and S) after the sprint 6 refactor.

**Status, 10 Nov:** HIVsim → EMOD′ is run (all five hybrids, analysis 16). EMOD′ ↔ Optima′ is now the **core path** (`core_question.md`), and its C step is done both ways (analysis 24). Its V, RN and D steps wait on a design decision, because Optima's transmission has no network. D6 puts X in the middle of Optima paths. The "Swappable when" column below is the 4 Nov plan.

| Path | Steps (proposed order) | Adapters needed for this order | Swappable when |
|---|---|---|---|
| **HIVsim → EMOD′** | C, V, RN, D | EMOD's condensed cascade on `sti.HIV` (**07**); EmodHIV on `StructuredSexual` and STIsim demography; EmodNetwork and EmodRisk on STIsim demography | C once 07 lands; the rest after two adapters (sprint 6–7). The **first milestone** (section 5) |
| **EMOD′ → HIVsim** (same order, from the other end) | C, V, RN, D | HIVsim's cascade as a condensed table (needs rate and stock-target pathways) on EmodHIV; `sti.HIV` on `EmodNetwork` and EMOD's demography; `StructuredSexual` on EMOD's demography | After the two pathway types and three more adapters: sprint 8 to January |
| **EMOD′ ↔ Optima′** (the published gap) | X (Optima′ ODE → agents) and measurement, then C, H, NT, DR | Optima′ as agents (the four engine-reading pieces); population label; `CascadeHost` on Optima′'s declared states; Optima's cascade as a table (both pathway types); EmodHIV giving infectiousness to Optima's force of infection | After Optima′ as agents (sprint 5) and the population label: January, or sprint 8 for C alone |
| **HIVsim ↔ Optima′** | As EMOD′ ↔ Optima′ | The same Optima pieces, with STIsim's side in place of EMOD's | Becomes cheap once either Optima path exists |

**Optima, ODE side (proposed, needs a decision).** Analysis 08 maps the condensed cascade into Optima′'s compartments. That allows a short **compartmental segment** next to the Optima endpoint, with the swaps that make sense in an ODE (EMOD's or HIVsim's cascade table and natural history parameters, as rates). Those runs are deterministic and take seconds, with no noise to pair. The order on Optima paths would then be: Optima′ ODE → [compartmental swaps] → X (to agents) → agent swaps → the other model. This puts X in the middle, not first or last, so it changes the 27 Oct decision. The alternative keeps X first and uses 08 only as a check of the cascade swap's ODE form.

### 2.4 How hybrids are named

A hybrid is a choice of (model, layer) for each slot, plus a pair of parameter sets and a seed. Example: `C=emod, V=hivsim, RN=hivsim, D=hivsim`. With layers split: `C.mech=emod, C.inputs=hivsim, C.params=hivsim, ...`.

### 2.5 Attribution

Let *v*(*S*) be an output for the hybrid with slots *S* from B and the rest from A. Each (parameter set pair, seed) is one **cell**; every hybrid in a cell uses the same seed, so differences are paired.

| Method | What it gives | Hybrids needed (n slots) | Cost at n = 4 / 8 |
|---|---|---|---|
| **Sequential** (one order) | Δₖ = *v*(first k) − *v*(first k−1); sums to *v*(B) − *v*(A) exactly, by endpoint identity | n + 1 | 5 / 9 |
| Both ends, same order, plus single swaps each way (`prompts.md` starting point) | Each slot's effect when swapped first and when swapped last; the difference measures its interactions | about 4n | 16 / 32 |
| **Shapley** (all orders averaged) | φᵢ = Σ over S ⊆ N∖{i} of \|S\|!(n−\|S\|−1)!/n! · [*v*(S ∪ {i}) − *v*(S)]; sums to *v*(B) − *v*(A); independent of order | 2ⁿ | 16 / 256 |
| Sampled orders (permutation sampling) | Unbiased estimate of Shapley | about m(n − 1) for m orders | for n ≥ 10 |

Proposed:
- **Shapley on the coarse slots.** With four slots, all 16 hybrids are the full lattice. With measured run times (section 4) that is about 16 core-hours for 10 pairs × 3 seeds, so about 10 minutes on the 160-core machine. Shapley is the summary number, and the figures show the sequential paths (waterfalls) as well. The real cost of the lattice is not compute: every slot of one model has to run with every combination of the other's, which needs the adapters of both directions in 2.3.
- **Layers and finer slots** (8–12 steps) are swapped along one or two sequential paths. Alternatively, they use nested (Owen) values within a coarse slot whose Shapley value is already known, which is cheap.
- **Uncertainty:** each cell's φ (or Δ), then the mean over cells and a bootstrap interval over cells. Report the spread over parameter-set pairs (parameter uncertainty) separately from the spread over seeds.
- **Pairing parameter sets:** HIVsim has ≥10 accepted runs (the stopping rule, `targets.md`), EMOD has 10 sets (those of analysis 04), and Optima has one. Pair them 1:1 at random, fixed once. A random pairing is an unbiased estimate of the average over all pairs. All Optima pairs use its single set.

## 3. Calibration and endpoint identity

**Tier 1, uncalibrated** (superseded 9 Nov: tier 2 throughout, below). A slot's parameters move with it, as the parameter layer. Tier 2 (compensated: refit one to three shared scalars, such as the transmission scale or ART uptake, to the primary targets after each swap, and report how far they move) is for January (`plan.md`). Full recalibration is out of scope.

**Decided** (Robyn, 9 Nov): the compensated results are the headline, and the uncalibrated ones and the difference (what calibration absorbs) are reported alongside. The knobs are a transmission scale and ART uptake. Proposed: ART uptake multiplies the rates into ART (testing and linkage) in flow cascades; a stock-target cascade (HIVsim's, Optima's) meets the number on ART by construction, so there it scales testing. Both endpoints are refitted the same way, so the refitted path runs between refitted endpoints; the difference between each published model and its refit is reported as its own component (calibration choices). The targets are the P3 stopping-rule targets (Robyn, 9 Nov). Every hybrid is compensated throughout, including for the 15 December interim (Robyn, 9 Nov: "calibrated pieces throughout"), so the refit is part of the runner from the first path.

**Endpoint identity** (**Decided**, `prompts.md`, 27 Oct): a fully swapped hybrid equals the rebuilt target model exactly. Check it by comparing results, every result array equal, with the same seed and module order (**Decided**, `prompts.md`, Starsim v4 notes; `sim.hash()` is gone). This is a unit test of the runner, separate from fidelity (each rebuilt model against its original, P1/P2). Proposed tests:

1. **Both endpoints**, for every path: `build(all A)` equals A's own `make_sim`, and `build(all B)` equals B's, for 3 cells. Short horizon (to 1995) on every change, and the full horizon once per milestone.
2. **Adapters are inert:** an adapter is inserted only when a consumer and its provider come from different models. Test 1 catches any that leak.
3. **Order of modules:** the runner places modules in one canonical slot order (D, R, N, V/T/H, S, C, M), and each model's `make_sim` is refactored to the same order and the same factories. Then identity holds by construction, and the test guards against regressions. The refactor itself is checked once for fidelity, in distribution, against the old `make_sim`.
4. **CRN check** (planned for sprint 3, not done): the variance of Δ₁ (the first step) with paired seeds against independent seeds, over 30 cells.

**Endpoints used in the swaps** (proposed): EMOD′ with the condensed cascade (within noise of the literal campaign, analyses 04–06); HIVsim native (from 1985, 10,000 agents); Optima′ as an ODE (checked against Optima, P2), with X to agents checked in sprint 5.

**Where common random numbers hold:**

| Model or module | Random numbers | CRN across hybrids |
|---|---|---|
| STIsim modules (`sti.HIV`, interventions) | Starsim distributions, per agent (per-agent CRN streams: Starsim PR #1538, `prompts.md`) | Yes, while agents and distribution names are the same |
| `StructuredSexual` | Starsim | Not when people are born or die (`plan.md`, risks) |
| Condensed cascade deliveries | Starsim (`coverage=`, `delay=`) | Yes; but `Holder.flow()` (VMMC, PrEP) uses its own NumPy stream, so not for those |
| EMOD′ (`EmodDemography`, `EmodRisk`, `EmodHIV`, `EmodCascade`) | One NumPy stream per module (`rand_seed + 1…7`) | Reproducible and separate by module, but sequential, so draws shift to other people as soon as the population differs. Weak CRN |
| `EmodNetwork` | Numba's global `np.random` (`nb_seed`) | Weak. Anything else using Numba's global stream would interleave with it (to check) |
| Optima′ ODE | Deterministic | Not needed |
| Optima′ as agents | Declared engine; draws are keyed by transitions | Changes whenever the set of transitions changes (`prompts.md`) |
| EMOD itself | AES counter per node | None: compared only in distribution |

Consequences:
- Endpoint identity needs only reproducibility, which every module has.
- CRN only lowers the variance of the Δs. If the CRN check shows little gain on EMOD′ steps, the fallback is more seeds (cheap, section 4). Moving EMOD′'s module streams to Starsim distributions is optional and would change EMOD′'s draws, so it would need fidelity rerun.

## 4. Runner and bookkeeping

Proposed layout: `common/swaps/` for the runner, `models/<model>/components.py` for each model's slot factories, and one analysis folder per path.

**Configuration** (a dict or YAML file per path):

```yaml
path: hivsim_to_emod
from: hivsim
to: emod
slots: [C, V, RN, D]            # order of the sequential path
lattice: true                   # also run every subset (Shapley)
pairs: {hivsim: [...10 accepted run ids], emod: [0, 25, ..., 225]}  # paired 1:1, fixed
seeds: 3                        # per pair; the cell seed is used by every hybrid in the cell
scale: 0.01
start: 1960.5
stop: 2025
scenario: null                  # 'interruption' etc. from January
```

**Build:** `build(choice, cell)` calls each chosen slot's factory, `models/<model>/components.py: make_<slot>(pars, context)`, in the canonical order. It inserts an adapter only between slots from different models, and sets `rand_seed` from the cell. Each model's own `make_sim` calls the same factories.

**Run:** a process pool, one Numba thread per process, as in analyses 05–06. Each hybrid is identified by a hash of (choice, pars versions, Starsim and STIsim pins, repo commit). A hybrid already run, for example an endpoint shared by two paths, is not rerun.

**Store:**
- Per run: the M outputs and the common-target measures, as one row per run and year, in `raw_results/swaps/<hash>/<cell>.parquet` (gitignored).
- A manifest (`raw_results/swaps/manifest.csv`): hash, choice, cell, seed, commit, pins, run time.
- Each path's analysis copies what its figures need into `analyses/NN_<path>/outputs/` (by the repo's rule).

**Compare:** per output, Δ per step and φ per slot with intervals over cells; waterfall plots per output; first-against-last position effects; the target scorecard of every hybrid.

**Run times and budget:**

| Model | Measured | Source |
|---|---|---|
| EMOD′, condensed cascade, 1% scale, to 2025 | 74 s per run (30 at once) | analysis 06 |
| EMOD′, compiled / literal campaign | 96 s / 60 s | analysis 05 |
| HIVsim, 10,000 agents, 1985–2021 | 10 s alone, 20–30 s with 150 at once | hivsim_malawi experiments 01–02 |
| EMOD′, 1960.5–1985 burn-in | 10 s (13 MB saved; 0.4 s to load; then 60 s to 2025) | analysis 11 |
| Optima′ ODE, dt 0.02 | Seconds (plus 30–75 s of R0 at `verbose>=1`) | `optima_prime_plan.md` |
| Optima′ as agents | See `docs/readouts/2026-11-10_sprint3/model_stats.md` | analysis 20 |

Budget per path, assuming about 2 min per hybrid run:

| Design | Hybrids | Runs (10 pairs × 3 seeds) | Core-hours | Wall time, 160 cores |
|---|---|---|---|---|
| Sequential, 4 coarse slots | 5 | 150 | 5 | a few minutes |
| Lattice, 4 coarse slots (exact Shapley) | 16 | 480 | 16 | about 10 min |
| Lattice, 8 slots | 256 | 7,680 | 256 | about 2 h |
| Lattice, 12 slots | 4,096 | 122,880 | 4,096 | too much: sample orders |
| Any of the above with 2 scenarios (January) | | × 3 | × 3 | |

Compute doesn't constrain this design; adapters do.

## 5. First milestone

**Plan** (`plan.md`): first path with endpoint identity by sprint 7 (1–7 December), with the swap runner and endpoint test passing on a trivial swap in sprint 6.

**Proposed: the path HIVsim → EMOD′, with the coarse slots in the order C, V, RN, D, uncalibrated.** This replaces HIVsim ↔ Optima′ as the first path.
- Since the plan was made, EMOD′ has gone from the critical path to done: it reproduces EMOD in analyses 02–06, with the condensed cascade and `CascadeHost`.
- Both EMOD′ and HIVsim are coded agent models.
- Any Optima path first needs Optima′ as agents, agent versions of four engine-reading pieces, the population label, `CascadeHost` on declared states, and Optima's cascade as a table. None of that exists yet.
- HIVsim ↔ Optima′ moves to January (or sprint 8 for its cascade step alone), together with EMOD′ ↔ Optima′, which shares all the Optima pieces.

**What has to be true:**
- HIVsim has reference parameter sets in `models/hivsim/reference/` (provisional ones are fine: endpoint identity doesn't depend on fit).
- One Starsim and STIsim pin runs both models. EMOD′'s condensed cascade needs Starsim after 4.0.0rc1 (ebd52c1a, with #1538), but `envs/starsim.yml` pins 4.0.0rc1.
- Both models are built from slot factories, and `build(all HIVsim)` and `build(all EMOD′)` equal their own `make_sim` exactly.
- The five hybrids of the path run, and the outputs analyzer gives every output in 2.2 except the ranking.

**Work, in order:**

| # | Work | Who (proposed) | Depends on |
|---|---|---|---|
| 0 | Decisions below | Robyn | Done, 9 Nov |
| 1 | Re-pin Starsim (post-#1538) and STIsim for both models; rerun the checks (`acceptance.py`, analysis 06 at a short horizon, an HIVsim run) | Cliff, Robyn | |
| 2 | HIVsim's reference sets copied in, with their hivsim_malawi commit | Robyn | P3 stopping rule (9 Nov) |
| 3 | Slot factories: `models/emod/v4/components.py` (from `harness.make_sim`, already modular) and `models/hivsim/components.py` (from `malawi.make_sim`); each `make_sim` calls them; identity tests | agent | 1 |
| 4 | EMOD′ burn-in to 1985, saved and resumed: resumed equals continuous exactly; run time measured. HIVsim at 10,000 agents against 1%, once | agent | 1, 3 |
| 5 | Analysis 07 lands: EMOD's condensed cascade on `sti.HIV`, with the `node` and `risk` the flow rows use (node = all) | other agent | 1 |
| 6 | Adapter: EmodHIV on `StructuredSexual` and STIsim demography (acts per step, condom share, type map; `emod_age`; pregnancy and birth events; seeding of HIGH risk read as `fsw` or `client`) | agent | 3 |
| 7 | Adapter: EmodNetwork and EmodRisk on STIsim demography (`emod_age`, nodes assigned by EMOD's node shares, debut) | agent | 3 |
| 8 | Outputs analyzer `common/outputs.py`, and `cascade_stage()` on both hosts; check STIsim's infection log gives sources | agent | 5 |
| 9 | Runner (`common/swaps/`): configuration, cells, pool, store, manifest; endpoint tests for both ends | agent, with Cliff | 3 |
| 10 | CRN check on step 1 (C) | agent | 8, 9 |
| 11 | Run the path, analysis `NN_hivsim_to_emod` (09 if 07 and 08 take those numbers): README before, SUMMARY after | agent, Robyn | all |

Then, for sprint 8: the other direction (EMOD′ → HIVsim, which needs the rate and stock-target pathways for HIVsim's cascade, and three more adapters), the 16-hybrid lattice and Shapley values, and the first Optima pieces.

**Risks:**

| Risk | Fallback |
|---|---|
| HIVsim's calibration isn't finished (deaths on ART too low, `targets.md` and exp 03) | Run on provisional sets; identity holds; rerun when final, since everything is scripted |
| Analysis 07 slips, or STIsim lacks pieces EMOD's table needs (WHO stage, symptoms, PMTCT per mother) | Derive WHO stage and symptoms from CD4; start the path at V with C held as HIVsim's (one step fewer) |
| Resuming a saved Starsim v4 sim doesn't reproduce a continuous run | Fix upstream in Starsim; meanwhile run EMOD′ from 1960.5 |
| EmodNetwork depends on EMOD-only people fields (single-precision ages, node, Risk-driven flags) more than seen | Swap R and N together with D (three coarse slots) |
| HIVsim at 10,000 agents differs from HIVsim at 1% beyond noise | Make the agent number its own step at the end of HIVsim paths |
| Weak CRN for EMOD′ modules makes Δs noisy | More seeds (cheap); optionally move EMOD′'s streams to Starsim distributions |
| Numba's global random stream is shared between EmodNetwork and other Numba code | Check in step 3; give EmodNetwork its own Numba generator state |
| Starsim re-pin changes HIVsim or Optima′ results | Rerun the fidelity checks (step 1); pin deliberately (**Decided**, `prompts.md`) |

## Decisions

| # | Decision | For | Proposed | Decided (Robyn, 9 Nov) |
|---|---|---|---|---|
| D1 | One time frame for all hybrids (from 1960.5), with HIVsim's and Optima's initial states applied as events in 1985 and 1990, so their swap endpoints are re-expressions checked in distribution | Robyn | Yes | Start all hybrids in 1985, with an EMOD′ burn-in from 1960.5 saved and resumed (1.4) |
| D2 | HIVsim's cascade expressed as a condensed-cascade table (with new rate and stock-target pathways), so all cascades share one mechanism; this also serves Optima's | Robyn, Cliff | Yes | Yes |
| D3 | First path HIVsim → EMOD′ instead of HIVsim ↔ Optima′ | Robyn | Yes | Yes |
| D4 | Order C, V, RN, D, and what "outside-in" meant | Robyn | As here; the lattice makes it a presentation choice | Yes |
| D5 | Shapley values on the coarse slots as the summary attribution, sequential paths for the figures | Robyn | Yes | Yes |
| D6 | Optima: a compartmental segment (with 08) next to the ODE endpoint, which moves X from first or last into the middle | Robyn | Discuss | Yes |
| D7 | Slot boundaries and interfaces (1.1, 1.3), as input to P4 | Robyn, Daniel, Cliff | For the sprint 4 checkpoint | Yes, for the sprint 4 checkpoint |
| D8 | Optima's `lost` split into lost before and after ART, for reporting only | Robyn | Yes | Yes |
| D9 | Scale 1% for all agent hybrids | Robyn | Yes | No: agent numbers vary with D; one check of HIVsim at 10,000 agents against 1% (1.4) |
| D10 | Whether EMOD's copied bugs (the dropout trap, mortality scale) also get a "fixed" variant as extra steps on EMOD paths. They are parameters or switches already, so the cost is small | Robyn, Daniel | After the first path | Yes: a fixed variant is needed. Built: `models/emod/v4/components_fixed.py` (Daniel's two fixes, PR #4); switches in `inputs.FIXES`; analysis 17 |
| D11 | Runs where the epidemic doesn't take off (EMOD's seeding on HIVsim's demography takes off in 2–3 of 5 runs at 10,000 agents; analyses 12, 15) | Robyn | Condition on take-off, make seeding its own step, or seed more | Condition on take-off (decided 9 Nov), **provisional, to revisit**: keeping only runs that take off can skew Δs. Take-off differs between hybrids (e.g. 5/5 vs 2/5), so the kept seeds differ between arms, which breaks the pairing, and whether the epidemic takes off is itself part of the swap's effect. Options when revisited: report take-off probability as its own output next to conditional Δs; condition on a cell's runs taking off in both arms; seed enough that take-off is near certain; or make seeding its own step (S). **Evidence (analysis 21):** share Δs are robust to the seeding rule (within 0.03, no sign changes); count Δs at V, RN and D are not, and conditioning doesn't fix it (it drops failed runs but keeps late starts: taken-off HEEE 27.8k new infections vs 43–47k with near-certain take-off). HIVsim's seeding as an S step, or 5× EMOD's, gives take-off 1 in every hybrid |
