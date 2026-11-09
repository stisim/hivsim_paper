# 22: Optima's population label, and `CascadeHost` on Optima′'s agents

**Question:** two of the three pieces analysis 20 left for hybrids on the EMOD′ ↔ Optima′ path (`docs/swap_design.md` §1.3, §1.4, §5):

1. **The population label.** As agents, Optima′ carries its 39 populations as a group (`pop`), moved by ageing and risk transfers, with no ages. When Optima′'s components sit on another model's people (EMOD′'s or HIVsim's demography and risk), each agent needs its Optima population from its age, sex and sex work, kept up to date. Can that label reproduce Optima's own population sizes, and where can't it?
2. **`CascadeHost` on Optima′'s declared states**, per agent, with `cascade_stage()`, so that `common/outputs.py` (`Outputs`) runs on Optima′ as agents and gives Figure 3A's categories (D8: lost split into before and after ART by whether the person has ever been on ART).

**Code:**
- `models/optima/v4/poplabel.py`: the 39 populations' definitions (Optima Malawi project, `models/optima/reference`), `label(age, female, fsw, client, msm)`, `stisim_label(sim)` and `emod_label(sim)`, and `PopLabel`, an analyzer that keeps each agent's label and records the derived sizes by population each whole year.
- `models/optima/v4/host.py`: `OptimaHost`, the `CascadeHost` on Optima′'s declared states (moves with `ss.move(uids=)`, agents only); ages and sex from each agent's population; `cascade_stage()`. Optima's transmission has no pairs (a force of infection by population), so the source of each new infection is drawn from the sources of the force of infection on the recipient's population (partner population × infected state, Optima's own split, as analysis 19's attribution), then an agent in that state and population.

**What's run** (`run.py`, env `hp_rc3`, Starsim 4.0.0rc3 from PyPI; HIVsim with `PYTHONPATH=/tmp/sti-rs-v4-art` as analyses 14–16):
1. `run.py label emod|hivsim`: EMOD′ (parameter set 0, 1% scale, from 1960.5) and HIVsim (its reference set, 10,000 agents, from 1985) to 2021 with `PopLabel`; derived sizes by population, 1990–2020, against Optima's (its forced `popsize`, which the ODE follows exactly, and Optima′ agents from analysis 20's 100k runs), scaled to people and as shares within sex. Expected: the general populations close where the models' demography agree; the key populations (FSW, clients, MSM) not, and the differences listed.
2. `run.py outputs SEED` (detached, a few seeds at 100,000 agents): Optima′ as agents (`agents.make_sim`) with `Outputs(OptimaHost)`; `run.py compare`: sexual transmissions by Fig 3A category in 2000–2020, against analysis 19's Optima′ ODE attribution (`analyses/19_fig3a_optima/outputs/prime_by_source.csv`) and totals against analysis 20's agent runs. Agreement: shares within agent noise (a point or two at 100k), totals within analysis 20's tolerance.
3. Checks that nothing else changed: `models/optima/v4/acceptance.py` (ODE acceptance) and the existing tests.
