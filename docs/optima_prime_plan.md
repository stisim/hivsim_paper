# Optima′: build plan

Optima′ re-expresses the Optima HIV Malawi project (`models/optima/original/Malawi_20220912.prj`) as a declared Starsim v4 model that runs compartmentally. It keeps Optima's handling of time (decided 30 Oct; see `optima_inventory.md`, "Time-step dependence"): Optima's step of 0.2 years is read as rate × 0.2, so each rate in Optima′ is Optima's per-step probability at dt = 0.2 divided by 0.2, with Optima's product structure.

## Acceptance

- Within 2% of the calibrated Optima run (dt = 0.2) for each output and year, **2000–2020** (the years with treatment data). Outputs: PLHIV, new infections, HIV deaths, diagnosed, on ART, suppressed, population size, overall and by population.
- Run as an ODE at a small step (0.02 y), with Optima's dt = 0.2 input series (`models/optima/original/optima_malawi_ref.obj`, interpolated between 0.2-year points).
- The step-independent harness (`models/optima/original/model_dt.py`, `FIX={'linear','inputs'}`) is the intermediate reference. At dt = 0.01 it is within 1.4% of the calibrated run for 2000–2019 and 2.5% at 2020. Optima′ at the same step should match it closely, component by component, because both compute the force of infection once per step from the start-of-step state.

## Structure

One declared disease, `optima`, stratified by `pop` (39 labelled groups, shares from Optima's 1990 population sizes):

- 38 states, Optima's own: `susreg`, `progcirc`, and 6 cascade states (undiagnosed, diagnosed, in care, on ART unsuppressed, on ART suppressed, lost) × 6 CD4 stages (acute, >500, >350, >200, >50, <50). 1,482 compartments.
- Transitions generated in a loop from Optima's per-step probabilities (`model.py` lines 79–355 and 776–829), including direct transitions where Optima allows two moves in one step (e.g. undiagnosed >500 → diagnosed >350).
- Rates that vary by group and year: `ss.by(pop={g: ss.timeseries(...)})`.

## Components, in build order

Each is checked against the harness before the next is added, with the components not yet built held as in Optima (e.g. treatment off, by zeroing `numtx` in both).

| # | Component | v4 mechanism | Notes |
|---|---|---|---|
| 1 | Demography: background deaths, births into F0-4/M0-4, ageing between bands, risk transfers FSW/Clients ↔ general | `ss.death`, weighted-count inflows, `ss.transfer` (time-varying) | Births per woman by mother group can't be declared (`engine.py:415`); use weighted counts. Risk transfers are balanced by group size (a function rate). The ageing band-skip of Optima's old code needs a direct transfer |
| 2 | Population forcing (`forcepopsize`, and reconciling groups with no inflows) | Rule each step | No public API for rescaling compartments: needs `ss.scale` upstream, or internals for now |
| 3 | Natural history: CD4 progression, HIV deaths, initial distribution | Transitions, `init=` | Progression reduced by the chance of dying in the step (product structure); every stage lasts at least one step at dt = 0.2 (acute 0.35 y) |
| 4 | Transmission: Optima's act-based force of infection with `force`, `inhomo`, condoms, circumcision, STIs, PrEP/PEP | Function rate `ss.by(pop=fn)` reading compartment counts | Held once per step, as in Optima. Reads per-group counts through engine internals (`engine.count`); a public per-group reader would help |
| 5 | Cascade: testing, linkage, loss to follow-up, return to care, AIDS testing, treatment failure | Transitions | `r*dt` rules are rates as given |
| 6 | ART: number on ART from `numtx` (lowest CD4 first before 2018), 63% suppressed at ART start, cascade proportions frozen from 2020 | Rule each step via `ss.move` / `ss.Delivery` | `ss.move` records no flow; new-on-ART flows need recording |
| 7 | Mother-to-child transmission and PMTCT | Inflows into infected states; rule for diagnosing mothers | |

## Starsim changes needed (upstream)

1. **Fixed** (Starsim #1498, merged). **Bug:** a rate given as a function by group (`ss.by(pop=fn)`) used as an inflow is multiplied by each group's share; a static `ss.by` is not (`engine.py:577`, also the check at `engine.py:418`). Reproduced: [100, 0, 50] becomes [20, 0, 25].
2. **Rescaling compartments:** a public way to scale the people in some states and groups by a factor, recording the change (for population forcing).
3. **Per-group counts in functions:** a public reader of compartment counts by group (now `engine.count(..., total=False)`, internal).
4. **Slow initialization:** the R0 shown at `verbose>=1` takes 30–75 s for a model of this size (`Compartmental.next_generation()`); make it cheap or skippable.
5. Nice to have: `ss.by` from a (groups × years) array.

## Open

- The 2020 step in the harness (new infections −1.0% to −2.0%, coinciding with the cascade proportions being frozen) is unexplained.
- Treatment beyond 2020 in the published MIHPSA runs: waiting on Rowan.
