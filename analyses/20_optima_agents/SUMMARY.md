# 20: Optima′ as agents, against its own ODE — summary

**Answer:** yes. Optima′ runs as agents and agrees with its ODE: at 400,000 agents (16 seeds), every output's mean over 2000–2020 is within 0.35% of the ODE and no year is more than 1.6% off; at 100,000 agents (48 seeds) the means are within 0.65%, with systematic gaps in 2017–2020 (2020: new infections +4.6%, HIV deaths +6.5%) that come from Optima's own target-driven ART rules and shrink with more agents. Optima′'s ODE acceptance is unchanged (`logs/acceptance_{before,after}.log`). The ODE ↔ agents switch (X) is now `agents.make_sim()` / `agents.make_sim(solver='ode')`.

## What needed agent versions (`models/optima/v4/agents.py`, new)

Six pieces read or wrote `engine.y` directly:
- **Four readers** (`Transmission`, `Births`, `Regain`, `VMMC`) now call `optima_prime.state_counts(sim, names)`: as an ODE it reads the compartments exactly as before; as agents it counts by population with `sim.groups.count()`, once a step. The only other change to `optima_prime.py` is a hook, `make_disease(extra=)` (17+/9−).
- **`CascadeTargets`** (number on ART from `numtx`, lowest CD4 first before 2018, proportions held from 2020): the same arithmetic, moving people with `ss.move(src, dst, p=share)`. Bug found and fixed: splitting ART starts into suppressed and unsuppressed as two moves in a row capped the second share above 1 when a whole CD4 stage started at once (64.2% suppressed instead of 63.2%, unsuppressed on ART 4–11% too low in 2015–20); each person now gets their own chances, the second conditional on the first.
- **`ForcePopsize`**: as agents, people can only enter or leave through declared transitions, so the forcing is declared `ss.immigration`/`ss.emigration` flows for each state not on ART, by population, sized from the counts at the start of the step so each population ends the step at its target.
- Everything uses public mechanisms, so the agent version also runs as an ODE: within 0.02% of the original ODE on every output. So the differences below come from running as agents, not from the rewritten pieces.

## Agents vs ODE, 2000–2020, every whole year (dt = 0.02)

Mean % = average over years of (agent mean / ODE − 1); z = (agent mean − ODE) / SE across seeds (`outputs/tables_100000.txt`, `outputs/tables_400000.txt`).

| Output | 100k agents, 48 seeds: mean % (max \|%\|), mean z | 400k agents, 16 seeds: mean % (max \|%\|), mean z |
|---|---|---|
| Population | 0.00 (0.00), −0.16 | 0.00 (0.00), −0.03 |
| PLHIV | +0.15 (0.31), +0.48 | −0.05 (0.16), −0.26 |
| Diagnosed | +0.19 (0.70), +0.68 | −0.10 (0.42), −0.52 |
| On ART | −0.06 (0.83), −0.18 | −0.04 (0.28), −0.48 |
| Suppressed | −0.07 (1.40), +0.02 | −0.34 (1.20), −0.58 |
| New infections, summed over each year | +0.50 (4.59 in 2020), +0.50 | −0.11 (1.32), −0.29 |
| Sexual infections, summed over each year | +0.65 (5.75 in 2020), +0.60 | −0.16 (1.28), −0.33 |
| HIV deaths, summed over each year | +0.37 (6.45 in 2020), +0.19 | −0.23 (1.61), −0.26 |

- **By population:** population size, PLHIV, diagnosed and on ART for all 39 populations: none has |z| > 3 in more than a quarter of the years, at either scale.
- **Single-step flows** (Optima's convention: the flow in the one step at each whole year) are very noisy as agents (SE 4–6% for infections, 5–9% for deaths); compare yearly sums.

**Systematic gaps at 100k, and why.** On ART −0.8% in 2017 and −0.7% in 2020; in 2020 new infections +4.6% (z 3.8) and HIV deaths +6.5%. They come from Optima's rules: ART is set by a target and taken from the people in care, so the in-care stock absorbs every random fluctuation with nothing pulling it back, and drifts across seeds like a random walk (2015: 51k–120k at 100k agents, against 88k in the ODE). When care runs out (2016–17), seeds that are short put fewer people on ART, but seeds with a surplus can't exceed the target, so the agent mean falls below the ODE. In 2020 Optima freezes its cascade proportions, carrying that spread into infections and deaths. At 400k the 2015 care spread narrows to 74–96k, on ART is −0.22% in 2017 and −0.16% in 2020, and 2020 infections and deaths are both +0.25%.

## Run time

ODE 47–55 s (1990–2021). Agents, 100k (~1% of Malawi in 1990; 208k alive by 2021): 180 s alone (22 s init), 245 s with 48 at once. Agents, 400k: ~415 s each with 16 at once. Agents, 100k at a monthly step: 85 s (not yet checked against dt = 0.02: the sprint 5 time-step study). Most of the time is per step (1,482 compartments and their time series), not per agent.

## Starsim friction (rc3), with proposed upstream fixes

1. **No count by group of a declared state that works under every solver.** As agents, `engine.count(states, total=False)` ignores `total` for states of people (`solvers/agents.py:1185`), and `sim.groups.count()` only works as agents. Proposed: `disease.count(state, by='pop')` under every solver. Workaround (marked): `state_counts()` checks which engine runs.
2. **No way to add people from code.** `ss.move()` refuses `ss.outside` for people (`transitions.py:1071`), so a rescaling rule can't be an intervention as agents. Proposed: allow `ss.move(ss.outside, state, n=, who=)` for people, or an `ss.scale(states, factor, who=)`. Workaround: declared flows in `agents.ForcePopsize`.
3. **`ss.move(n=)` rounds to whole agents** (`solvers/agents.py:730`): fewer than half an agent's worth moves nobody. Used `p=` instead. Proposed: stochastic rounding.
4. **Splitting one state across two destinations needs a conditional second share** (easy to get wrong). Proposed: `ss.move(src, {dst1: p1, dst2: p2})` with one draw per agent.

## What's left for the path (EMOD′ ↔ Optima′)

- **The population label:** as agents, Optima′ carries the 39 populations as a group dimension (`pop`), moved by ageing and risk transfers, with no actual ages. Hybrids with EMOD′'s agents need each agent's Optima population worked out from age, sex and risk group, and kept up to date.
- **`CascadeHost` on declared states**, per agent via `ss.move(uids=)` (exists, agents only).
- **Optima's cascade as a table**, for both pathway types.
