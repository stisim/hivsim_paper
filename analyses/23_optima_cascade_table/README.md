# 23: Optima's cascade as a table

**Question:** can Optima's own cascade be written as a `common/cascade` table, so that the cascade slot (C) can be swapped between EMOD′ and Optima′ through the shared `CascadeHost` interface? This is the last piece analysis 22 (§3) lists before the EMOD′ ↔ Optima′ path (`docs/core_question.md`: Bansi-Matharu 2025 Figure 3A, Malawi 2020), and decision D2 of `docs/swap_design.md` (a cascade as a table, with rate and stock-target pathways).

**What's added** (`common/cascade`, minimal and general):
- **Rate pathways:** a rate per person-year, by group (sex, ages, CD4 band, Optima population) and year, for people in a given place in care (undx, dx, care, lost, svl, usvl), moving them to another place; or a number a year spread over the people in some places (Optima's viral-load monitoring).
- **Stock-target pathways:** the number on ART by year, filled from those in care by a priority rule (lowest CD4 first before a year) and removed back to care when over; and shares held from a year (diagnosed, in care, on ART, suppressed), as Optima's `CascadeTargets`.
- Both in agent form (`common/cascade/pathways.py`, through `CascadeHost`, with new host moves `link_to_care`, `lose_from_care`, `switch_arm` and a `care_place` reader, defaults no-op so EMOD′ and HIVsim are unchanged) and compartmental form (`rates.py`, for the Optima′ ODE).
- `models/optima/v4/extract_cascade.py`: Optima's Malawi rows from the project's inputs into `common/cascade/tables/malawi_optima_params.csv`.

**Checks** (`run.py`, env `hp_rc3`, dt = 0.02, 1990–2021):
1. **ODE:** Optima′ running the table against Optima′ with its own cascade: PLHIV, diagnosed, on ART, suppressed, new infections, HIV deaths 2000–2020 and the Figure 3A source shares. Expected: exact to rounding (the table holds Optima's series at its time points), at least within 1%.
2. **Agents:** Optima′ as agents (100,000 agents, 8+ seeds) with the table through `OptimaHost` against analysis 20/22's runs with its own cascade (same seeds): within Monte Carlo noise (|z| < ~2 for most years).
3. **Others unchanged:** EMOD′ and HIVsim identity tests (`tests/`, `common/swaps/tests`), and Optima′'s ODE acceptance with its own cascade.
4. If time allows: a smoke run of the table on HIVsim's or EMOD′-fixed's people (runs; not interpreted).
