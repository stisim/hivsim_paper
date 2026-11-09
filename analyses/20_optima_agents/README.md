# 20: Optima′ as agents

**Question:** does Optima′ (`models/optima/v4/`, verified as an ODE against Optima within ~2%, 2000–2020) run with Starsim's agent solver, and do the agents agree with its ODE? This is the X slot of the swap design (`docs/swap_design.md` §1.1: the ODE-to-ABM switch is a component of its own, and every hybrid on an Optima path runs as agents), and the first of the pieces §5 lists for the EMOD′ ↔ Optima′ path.

**What needed changing.** Six pieces of Optima′ read or write the ODE engine's compartments (`engine.y`) directly:
- Readers (rates worked out once a step from the counts at its start): `Transmission` (force of infection), `Births` (births and mother-to-child transmission), `Regain` (regaining suppression, from the number on ART), `VMMC` (table cascade only). They now call `optima_prime.state_counts(sim, names)` (states × populations), which as an ODE reads the compartments as before and as agents counts the agents by population (`sim.groups.count()`, public), once a step (`agents.counts()`).
- Writers (rules applied after each step): `ForcePopsize` (Optima's population forcing) and `CascadeTargets` (number on ART from `numtx`, lowest CD4 first before 2018, proportions held from their fix years). Agent versions in `models/optima/v4/agents.py`:
  - `CascadeTargets`: the same arithmetic, the people moved with `ss.move(src, dst, p=share)` (public; each agent moves with the probability, so the expected number is Optima's and the split among populations is in proportion, as Optima's is).
  - `ForcePopsize`: as agents people can only enter or leave through declared transitions (`ss.move()` refuses `ss.outside` for people), so the forcing is declared flows: `ss.immigration` into and `ss.emigration` out of each state not on ART, by population, sized from the counts at the start of the step so that each population ends the step at its target (as the ODE's rescaling after each step makes it).
  Both use only public mechanisms, so the agent version also runs as an ODE (`agents.make_sim(solver='ode')`): that run separates what the agent versions of the pieces change (deterministic, against the original ODE) from what running as agents changes (stochastic, against it).

**What's run** (`run.py`, env `hp_rc3`, Starsim 4.0.0rc3 from PyPI; Optima's own cascade, all components, dt = 0.02 as the ODE):
1. `run.py ode`: the original ODE (`optima_prime.make_sim`) and the agent formulation as an ODE (`agents.make_sim(solver='ode')`), with run times. Expected: within a fraction of a percent of each other.
2. `run.py abm`: 20+ seeds as agents (scale set by run time; target ~1% of Malawi, ~100k agents), each saved to `outputs/abm_<seed>.obj`, with run time.
3. `run.py compare`: for 2000–2020, each whole year: PLHIV, new infections (sexual + mother-to-child, the step at each whole year and the year's mean), HIV deaths (as Optima reports them), diagnosed, on ART, population, overall and by population. ABM mean ± its standard error (over seeds) against the ODE; z = (ABM mean − ODE)/SE. Agreement: |z| < ~2 for most points, or a systematic gap explained (e.g. integer agents in small populations among the 39, stochastic rounding, the forcing's one-step lag).
4. Optima′'s ODE acceptance (`models/optima/v4/acceptance.py`) rerun after the change to `optima_prime.py`: must be unchanged.

Flows as agents are counted in whole agents per step; annual numbers are compared as the mean over the year's steps as well as the step at the whole year, because the single step is noisy.
