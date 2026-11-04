# Prompts

Two parts:

1. **Plan:** the prompts we expect to run, in order, each with its inputs, outputs, acceptance check and human checkpoint. Drafts, refined before running.
2. **Record:** prompts as actually run, copied exactly as typed, each with a short summary of the response and the commit. Follows `starsim/v4_design/admin/prompts.md`.

## Part 1: Plan (draft)

**Decided (27 Oct):** two separate checks. *Fidelity:* each rebuilt model matches its original, within tolerance. *Endpoint identity:* a fully swapped hybrid equals the rebuilt target model exactly.

**Decided (27 Oct): the decision scenario is ART treatment interruption.**
- It is defined at interface level from 2026, so every hybrid implements it through its own cascade component.
- It makes interruption and re-engagement central components in P4.
- HIVsim currently has no explicit re-engagement rate. Re-engagement only happens through coverage correction, so this needs building.
- MIHPSA Malawi gives a published comparison point: 8.4–20.1% of transmission from people who interrupted treatment.

**Decided (27 Oct): attribution is computed at the path endpoints only.**
- Each endpoint gets the attribution method its own model used.
- Measurement is reported as one explicit step, so the share of the gap that is purely definitional is a result in its own right.

**Options recorded, if the endpoints show a large gap from definitions alone:**
- (a) compute all three methods (tracked source, cascade-weighted, counterfactual zeroing) on every hybrid;
- (b) compute all three on a few key hybrids only;
- (c) compute only tracked source and cascade-weighted on every hybrid. Both are cheap; counterfactual zeroing needs one extra run per cascade stage.

**Open: swap design.** The proposed starting point for each pair:
1. Each single B component swapped into A.
2. Each single A component swapped into B.
3. One full path from A to B, in outside-in order.
4. The same path in reverse.

Revisit once underway.

**Decided (27 Oct): each functional block has three layers:**
- **mechanism** (the code);
- **inputs** (historical program and demographic data);
- **parameters** (calibrated values).

Each layer is swapped separately. By default, all input swaps form one step, which is split out only where it matters.

Optima's costing module was not used in MIHPSA and is out of scope, so Optima is treated as an epidemiological model only.

**Decided (27 Oct): rebuild agents (P1, P2) run on their own in a test-and-fix loop,** with the fidelity checks as pass/fail. Hard rules:
1. **No changing calibrated parameter values to reach tolerance.** Only structure and logic may be fixed.
2. **No changing the shared interfaces.** A model that doesn't fit an interface stops the agent and escalates to Robyn, Daniel and Cliff.
3. **Log every judgement call** wherever the original code or config is ambiguous.
4. **At most ~5 fix attempts per component,** then escalate.

The model owner signs off at the end: Daniel for EMOD, Robyn for Optima.

**Decided (27 Oct):** all hybrids run as agent-based models. Rebuilt Optima is verified as an ODE against the original Optima, and as an ABM against its own ODE form. The ODE-to-ABM switch is treated as a component of its own, swapped first or last on every Optima path.

**Decided (27 Oct):** the component map comes before any re-expression, so that every model is built against the same component interfaces.

| # | Prompt (short) | Inputs | Output | Acceptance | Human checkpoint |
|---|---|---|---|---|---|
| P0 | **Humans, not agents.** Gather the originals and export their reference outputs | Daniel: EMOD-HIV Malawi (config, campaign, demographics, binary version, calibrated parameter sets). Robyn: Optima HIV Malawi project | `raw_results/originals/` (gitignored); reference outputs summarised into `results/` | both originals rerun and match their MIHPSA outputs | Done by Daniel and Robyn |
| P4 | List the atomic components of each model and propose shared interfaces | EMOD source + Malawi config, Optima source + Malawi project, HIVsim/stisim, Starsim v4 | component map + interface spec | each model's behaviour accounted for by exactly one component | **Approve boundaries:** Robyn, Daniel (EMOD), Cliff (v4) |
| P1 | Re-express EMOD-HIV (Malawi) in Starsim v4, against the interfaces | P4 spec, EMOD source + Malawi config, Starsim v4, stisim | `models/emod/` | for each output and year 2000–2040, the rebuilt model's median lies inside the original's 90% interval across seeds | TBD |
| P2 | Re-express Optima HIV (Malawi) in Starsim v4, against the interfaces | P4 spec, Optima source + Malawi project, Starsim v4, stisim | `models/optima/` | within 2% of the original Optima, for each output and year 2000–2020 (the years with treatment data; decided 30 Oct) | Robyn signs off (sprint 4). First pass done 30 Oct: `models/optima/acceptance.py` |
| P3 | Build a common target set from the data EMOD and Optima actually used, then calibrate HIVsim Malawi to it (calib plugin); score the EMOD and Optima originals on the same set | P0 calibration inputs, Bansi-Matharu 2025 + supplement, calib plugin | `common/targets.csv`, priors, fitted HIVsim, a fit scorecard for all three models | HIVsim fits the common set; scorecard complete | **Approve target set** (Robyn) before calibrating. Approved 2 Nov; calibration in `experiments/` |
| P5 | Refactor HIVsim to the interfaces; design and run the sequential swaps | P1–P4 outputs | swap plan, path results | endpoint identity is exact (unit test, with common random numbers) | TBD |

### P1 (draft)

```text
Implement EMOD in the Starsim architecture.
```

- Decided (Robyn, 3 Nov): EMOD′ copies the bugs in Daniel's model (background mortality scale, ART dropout trap, the male non-suppressive mortality-table row, the ignored ANC table, the no-op calibrated parameters), so it reproduces the published model. Fixes can be switches, off by default.

### P2 (draft)

```text
Implement Optima HIV in the Starsim architecture.
```

### P3 (draft)

```text
Extract all the calibration targets and adjustable parameters from the paper. Use the calib plugin and find a fit.
```

### P4 (draft)

```text
List out the atomic components of each model and design the sequential swapping.
```

## Starsim v4 notes (reviewed 28 Oct, `origin/v4` = 1ce69a13)

- Upstream from the EMOD′ work (4 Nov): PR #1538 (cycles through delayed deliveries; per-agent CRN streams), merged. Issue #1542 lists the remaining Delivery gaps the campaign compiler worked around (ordering between deliveries, leaving a wait when ineligible, public forward references, capacity by step, and smaller items).

**Update (3 Nov).**
- The pin is now Starsim 4.0.0rc1 (PyPI). It needs Python 3.12 or later (on 3.11 its Rust kernels silently fall back to Numba) and works with sciris 3.7.0.
- The declared/coded spike is done (`spikes/seam/`, 29 Oct; removed 3 Nov, last in commit a751284): a declared cascade runs on HIVsim's coded HIV, via `ss.move()` along a transition (Starsim #1492, merged).
- HIVsim's cascade is coded, not declared: STIsim #607 splits ART from HIV and adds optional Optima-style care, so swapping the cascade is a change of settings.

- **Base.**
  - Build on `origin/v4`, pinned to a commit, and rerun the endpoint-identity test on every re-pin.
  - Since 27 Oct `v4` has Rust kernels (same numbers as Numba, bit for bit; no effect on agent runs) and a reorg: `ss.interop` is now `ss.io`, with `ss.import_from()`/`ss.export_to()`.
  - It needs `sciris>=3.5` (with 3.4 or older, `ss.Sim` crashes reading its own results), although `pyproject.toml` still says `>=3.2.8`.
  - Cliff's STIsim port is stisim `origin/v4` (HIVsim is the `hivsim` package in the same repo). It is written as code (`ss.Infection`), not declared: see Performance.
- **Which models need the declared form.**
  - Only Optima′: it is the one verified as an ODE (against the original) and as agents (against its own ODE), and the declared form gives both from one model.
  - HIVsim and EMOD′ can be written as code. EMOD is individual-based anyway.
  - What every path needs is declared and coded components working together behind the P4 interfaces in one agent sim (e.g. Optima′'s declared cascade on HIVsim's coded natural history). This is untested, so spike it before sprint 3.
  - Optima is legacy Optima HIV ([optimamodel/optima](https://github.com/optimamodel/optima)), not Atomica, so the Atomica importer doesn't apply: Optima′ is a hand rebuild.
- **Interfaces.** Starsim doesn't enforce module interfaces (no reads/writes declarations). The P4 interface contract is our own layer, enforced by our own tests. State it in agent terms (e.g. CD4 category and ART status per person, `start_art(uids)`, a death hazard), so declared and coded components can each implement it.
- **ODE ↔ ABM.**
  - Declared models without individual networks run under every solver.
  - `SexualNet` now has a compartmental approximation (partner-change rates between groups), marked as approximate.
  - Clocks, `ss.attr`, queues, `delay=`, and `.count`/`.since` of events are agents only, as is any coded module that keeps a value per person (`define_states()`).
  - Optima′ should use mixing pools or `ss.contacts`, with CD4 × cascade as declared states.
  - Agents take whole time steps, so the ABM-vs-ODE check needs a dt study, not a tight tolerance.
- **CRN.** Probably not a major issue: HIV outputs should vary little between runs at a reasonable population size. Still worth the cheap check in sprint 3 (paired vs independent variance on one real swap). `SexualNet` isn't CRN-safe when people are born or die, and in the declared engine draw keys follow the transitions and rounds of moves, so engine changes alter results between pins.
- **Endpoint identity.** `sim.hash()`, `ss.diff()`, `to_spec()` and `to_code()` were removed on 27 Oct. Check identity by comparing results: every result array equal, same seed and module order.
- **HIV pieces.**
  - Directional `ss.peract` is exact.
  - Modifiers multiply the hazard, not the probability (about a 1% difference from EMOD's convention for acute-phase multipliers).
  - `SexualNet` has one partnership type per network; STIsim's `StructuredSexual` has several.
  - The stock ART target (`target=`) and `ss.Resource` were removed on 20 Oct, and there is no `ss.program()`. A target is `deliver=` or `capacity=` as a function of the sim.
  - A coded module's parameters can't be targeted by `ss.modify`/`ss.change`, only its `rel_sus` and `rel_trans`.
  - Starsim 3's `ss.Calibration` is gone, replaced by `sim.fit()`. Check that the calib plugin (P3) works with it.
  - There is no EMOD importer; it's listed in `v4_design/todo.md`.
- **Performance** (28 Oct, STIsim gonorrhea, 5,000 agents, 20 years monthly).
  - Declared: 0.80 s; coded: 0.58 s; the network alone is 0.35 s. So the declared disease costs about twice the coded one, which is why STIsim's diseases are written as code.
  - The overhead is per-batch Python in the agent engine, not the random draws (the same in number and cost), and is not fixed by the Rust kernels. Optima′ as agents would use mixing pools, a different (faster) transmission path, so measure it on an Optima-like model before planning engine work for the paper.
- **Ask Cliff:**
  - Is the EMOD importer planned?
  - Can transition and distribution naming be kept stable, for CRN?
  - Would he consider optional interface declarations?
  - Could we have a frozen API window for the paper?
  - *Answered:* stock ART targets stay in STIsim (removed from Starsim, 20 Oct); `SexualNet` has a compartmental approximation.
- **Malawi files:** ask Rowan for the Optima Malawi project.

## Part 2: Record

*None yet.*
