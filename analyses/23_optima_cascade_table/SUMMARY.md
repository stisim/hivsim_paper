# 23: Optima's cascade as a table — summary

**Answer:** yes. Optima's own cascade is now a `common/cascade` table (`common/cascade/tables/malawi_optima_params.csv`, 35 rows). It uses two new pathway types, rates and levels, each in agent form (through `CascadeHost`) and compartmental form. **As an ODE**, Optima′ running the table reproduces Optima′ with its own cascade **to rounding** (largest difference 8e-13 % on any output 2000–2020, 6e-13 points on the Figure 3A shares). **As agents** (100k agents, 32 seeds), Optima′ running the table through `OptimaHost` agrees with its own cascade to within about 1%. On ART, suppressed, diagnosed and HIV deaths are within 0.4% (mean over 2000–2020). But PLHIV (−0.9% against the ODE) and new infections (−1.2%) run slightly low, at z ≈ −2.4 over 32 seeds; the own cascade's 48 seeds sit at +0.15% and +0.5%. This is not explained (see Open). EMOD′, HIVsim, the EMOD table on Optima′ and Optima′'s own cascade (ODE acceptance, and as agents) are unchanged, bit for bit.

## What's added

**Table columns** (`params.py`, optional, so the EMOD tables still read as before):
- `group`: `sex`, `age`, `cd4` band, `pop` (an Optima population), e.g. `cd4:gt50|lt50; pop:FSW|M15-19`;
- `places`: `src|src>dst` over undx, dx, care, lost, svl, usvl;
- `per`: for a rate given as a number a year, the places it is spread over.

**Pathway types:**
- `rate`: a rate per person-year by group and year, for people in given places, to another place. Or, with `per`, people a year spread over the people in those places (viral-load monitoring).
- `target` (`on_art`): the number on ART by year, filled from those in care. CD4 bands go lowest first before `lowest_cd4_before`; `p_suppressive` start suppressed. People are moved back to care when the number is over target.
- `hold` (`diagnosed`, `in_care`, `on_art`, `suppressed`): shares held from a year at their value then.
- `pmtct`/`prep`/`pep`: coverages a model applies itself.

**Agent form:** `common/cascade/pathways.py`, `TableCascade(params, host)`. Each step:
- People in the rate rows' places move as competing risks.
- Then the target and holds are met as each person's chance, using the arithmetic of `agents.CascadeTargets`.
- Moves go through the host: `test`, `link_to_care`, `lose_from_care`, `start_art`, `stop_art`, `switch_arm`.

**Compartmental form:** `rates.py`, `CompartmentalCascade.path_rates` / `path_numbers` / `targets` / `coverage()`. The condensed-cascade build now runs only if the table has `hct` rows. Optima′ (`optima_prime.make_sim(cascade=Params(...))`) handles a table of pathways as follows:
- It takes the rates with its CD4 moves, as Optima does (`table_probs`, `PRODUCT`).
- It reads regain from the `per` row and the target and holds from the table (`CascadeTargets(ref, cx=)`).
- It reads PMTCT, PrEP and PEP into its births and force of infection (`table_inputs`).
- Rates for other moves are added as transitions.

**`CascadeHost` additions** (`interface.py`). Defaults are no-ops or derived, and the existing cascade calls none of them:
- `care_place()`: default None, meaning the cascade keeps its own record of care and lost;
- `cd4_band()`: default from `cd4()`;
- `population()`: default raises;
- `link_to_care()`, `lose_from_care()`, `switch_arm()`: default no-ops;
- `rate_scale()`: default 1.

Per host:
- `OptimaHost` implements all of them on its states. `rate_scale()` is the survival of Optima's 0.2-year step (see the check).
- `StisimHost` and `EmodHost` implement only `population()` (from `poplabel`, at call time).
- HIVsim's `components.make_C` gets a `cascade='table'` variant.
- `agents.make_sim(cascade=Params)` builds Optima′ with the table's coverages and without its own cascade moves or rules.

## The table's rows (`models/optima/v4/extract_cascade.py`, from `simpars` at dt = 0.2)

Rates are Optima′'s: Optima's probability per 0.2-year step divided by 0.2. Series are kept at Optima's time points (points lying on a straight line are dropped). Where Optima combines rates for people with AIDS (CD4 < 200), the combined rate is its own row.

| Row | Places | Groups | Value |
|---|---|---|---|
| hivtest | undx>dx | non-AIDS bands; 11 sets of populations | hivtest |
| aidstest | undx>dx | AIDS bands | max(aidstest, hivtest) |
| linktocare / aidslinktocare | dx>care | non-AIDS / AIDS | (1 − exp(−0.2/years))/0.2; AIDS: the larger of the two |
| leavecare / aidsleavecare | care\|svl\|usvl>lost | non-AIDS / AIDS | leavecare; AIDS: min(aidsleavecare, leavecare) (Optima's loss from ART is leavecare) |
| returntocare / aidsreturntocare | lost>care | non-AIDS / AIDS | returntocare; AIDS: (1 − (1 − aidstest·0.2)(1 − returntocare·0.2))/0.2 |
| treatfail | svl>usvl | | 0.02 until fixpropsupp (2020), then 0 |
| regainvs | usvl>svl, per svl\|usvl | | regainvs × numvlmon a year, until 2020 |
| target on_art, lowest_cd4_before, p_suppressive | | | numtx; 2018 (allcd4eligibletx); 1 − exp(−0.2/treatvs) = 0.632 |
| hold diagnosed, in_care, on_art, suppressed | | | fixpropdx 2100, fixpropcare 2100, fixproptx 2020, fixpropsupp 2020 |
| pmtct number, efficacy; prep / pep coverage (by population), efficacy | | | numpmtct, effpmtct 0.9; prep, effprep 0.95; pep, effpep 0.73 |

## Check (`run.py`; `outputs/compare_100000.txt`, `outputs/fig3a_100000.txt`)

**ODE** (dt 0.02, 1990–2021): every output matches to rounding over 2000–2020.

| Output | Max difference |
|---|---|
| PLHIV, diagnosed | 1.1e-13 % |
| On ART | 2.8e-13 % |
| Suppressed | 2.4e-13 % |
| New infections | 8.2e-13 % |
| HIV deaths | 1.0e-13 % |
| Figure 3A shares | 6e-13 points |

In 2020 both give 42.13 / 39.35 / 12.73 / 5.79 (undiagnosed / diagnosed never treated / on ART / lost). Run time is the same (49 s).

**Agents** (100k agents; 32 seeds against analysis 20's own-cascade runs with the same seeds; the seeds are uncorrelated between the two). Mean over 2000–2020:

| Output | Table vs own | Table vs ODE | Own (analysis 20, 48 seeds) vs ODE |
|---|---|---|---|
| PLHIV | −0.75% (mean z −1.3) | −0.87% | +0.15% |
| Diagnosed | −0.39% (z −0.7) | −0.40% | +0.19% |
| On ART | −0.04% (z +0.5) | −0.11% | −0.06% |
| Suppressed | −0.30% (z +0.4) | −0.26% | −0.07% |
| New infections | −1.34% (z −1.0) | −1.16% | +0.50% |
| HIV deaths | −0.25% (z −0.2) | −0.36% | +0.37% |

- **First round, without `rate_scale`** (16 seeds, `outputs/round1/`): diagnosed were +1.4–2.0% in 2000–2008 (z ≈ 3). The cause is that Optima′ takes cascade moves only among those who survive its 0.2-year step (survival (1 − HIV death)(1 − background death), 0.93 at CD4 < 50). The host now gives that survival as `rate_scale()`, and the gap is gone.
- **Figure 3A sources as agents** (16 seeds with `Outputs(OptimaHost)`, against analysis 22's 16): mean |table − own| is 0.75 points, with 15% of category-years at |z| > 2.
  - In 2020 the table gives 39.8 / 43.2 / 11.4 / 5.6, against 43.3 / 39.4 / 12.2 / 5.1 for its own cascade and 42.1 / 39.3 / 12.7 / 5.8 for the ODE (z ±3.3 on undiagnosed and diagnosed never treated).
  - 2020 is the year in which analysis 20 found the agents' ART-target noise largest (care stock random walk, proportions frozen).
- **Run time:** 316 s against 244 s. `care_place()` loops over 36 states each call.

**Unchanged:**
- `tests/` and `common/swaps/tests`: 22 passed (`logs/tests.log`, with STIsim `/tmp/sti-613` and Starsim `/tmp/ss-nov09` on PYTHONPATH).
- Optima′ ODE acceptance with its own cascade: identical to analysis 22's (`logs/acceptance_after.log`).
- The EMOD table on the Optima′ ODE (analysis 08's path): bit-identical to starsim-rc3, to 2012.
- Optima′ agents with their own cascade, seed 1: identical to analysis 20's (`logs/own_seed1_identity.log`).
- New `tests/test_cascade_table.py` (2 tests): the ODE table equals its own cascade to 2010; as agents, the number on ART is met.

**Smoke run on HIVsim** (`run.py smoke`, 10k agents, `outputs/smoke_hivsim.txt`): it runs (10 s). In 2021 there are 890k on ART against the target of 864k; care 144k, lost 25k. Not interpreted.

## Open

The agents' small deficit in PLHIV and new infections (about 1%, z ≈ −2.4 over 32 seeds, from about 1995), and 2020's Figure 3A split, are not explained. Checked and ruled out:
- populations mapped wrongly (PLHIV and diagnosed by population look like the own cascade's);
- agents dying in the step (none are flagged).

Candidates:
- the order within a step: the table's moves come after the disease step, while the own cascade's moves are inside it, with the CD4 moves;
- Optima's products, such as failure and regain conditional on not being lost ((1 − loss) factors) and the moves among survivors in other forms;
- noise: the first 16 seeds alone showed no gap.

To settle it: 64 seeds, or 400k agents; or log each move by state and step against the own cascade's transition results (`/tmp/hp-optcx-scratch/dbg.py` has the start of this).

## Provisional decisions (for Robyn)

1. **Rates are Optima′'s per-step convention:** Optima's probability per 0.2 years divided by 0.2. They are not converted to continuous-time rates (−ln(1 − q)/0.2), so that the ODE is exact. On another host they are read as rates per year.
2. **Optima's AIDS combinations** (max, min, 1 − (1 − a)(1 − b)) are rows of their own for the AIDS bands, not a rule for combining rows. Rows on the same move add.
3. **Loss from ART is Optima's `leavecare`** (from care, svl and usvl to lost). Optima's "lost" covers both those lost from pre-ART care and those who stopped ART.
4. **Optima's "in care" is pre-ART care.** It is the place `link_to_care` leads to, and the number on ART is filled from it. EMOD's pre-ART, staging and linking, and HIVsim (which has none), aren't the same thing. On EMOD′ and HIVsim, care and lost exist only in the table cascade's own record, and the host doesn't see them (`link_to_care` and `lose_from_care` are no-ops there).
5. **`switch_arm` (failure, regain) is a no-op on EMOD′ and HIVsim:** they keep the arm drawn at the start.
6. **Over-target ART is moved to care** (Optima's rule), as `stop_art` then `link_to_care`. A held `diagnosed` share below its level would move dx back to undx; no host does this (it raises). Optima's never applies (2100).
7. **PMTCT, PrEP and PEP are coverages that the model applies itself.** Optima′ puts them in its births and force of infection, read from the table. `TableCascade` skips them on agent hosts; drawing a per-person status from the coverage is the alternative.
8. **Optima has no mothers:** there are no 12-week or infant events in Optima′, so the table has none.
9. **`rate_scale()`, survival of a model's longer step:** a host hook to reproduce Optima's order (death first, then moves among survivors). The default is 1. Alternatively, put the survival into the rates.
10. **Groups by population on other hosts:** `population()` labels at call time (`poplabel`, not the end-of-step `PopLabel` analyzer). The key-population mismatches of analysis 22 §1 apply.
11. **The number on ART is in people** (Optima's), scaled by `sim.current_scale`. On HIVsim this is Optima's Malawi numbers on HIVsim's population (about the same size).

## Starsim friction (rc3)

1. **States are inferred to be counts when nothing links them to where people start** (`starsim/internals/graph.py:503`, `check_counts`). Without the cascade's declared transitions (moves done by an intervention with `ss.move(uids=)`), Optima′'s care states fail this check. Workaround: `ss.states(..., people=True)` when the moves are external (`optima_prime.make_disease`). Proposed fix: count `ss.move` as a link, or a disease-level `people=True`.
2. **Still no per-agent "which state"** (analysis 22's #2). `care_place()` and `cd4_band()` loop over 36 `BoolState`s, about 25% of the agent run time. Proposed: `disease.state_of(uids)`.
3. **`ss.move(uids=)` moves from one source** (analysis 22's #3). Proposed: `ss.move({src: dst}, uids=)`.

**Files:**
- New: `common/cascade/pathways.py`; `common/cascade/tables/malawi_optima_params.csv`; `models/optima/v4/extract_cascade.py`; `tests/test_cascade_table.py`; this analysis (`run.py`, `launch_abm.sh`, `outputs/`, `logs/`).
- Changed: `common/cascade/{params,rates,interface}.py` and `README.md`; `models/optima/v4/{optima_prime,agents,host}.py`; `models/hivsim/{stisim_host,components}.py`; `models/emod/v4/cascade_host.py`.
