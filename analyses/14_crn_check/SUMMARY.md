# 14: Common random numbers on step 1 (C) — summary

**Answer:** pairing (the same seed for a model and its one-slot swap) helps stocks a lot and flows a little. It does nothing for ART coverage or transmission shares. Paired runs are identical until the cascades first act differently (1990–91), then decorrelate within a few years. Part of that is a leak in Starsim's agent engine (below), which should be fixed upstream and this check rerun.

**Setup:**
- HIVsim (all slots HIVsim) against HIVsim with EMOD's condensed table as C (`hivsim_cx.py`, which is what `swap.make_sim('condensed')` builds).
- Control: C from `hivsim_b`, HIVsim's own factory file under another name.
- 60 seeds per arm, through the runner's pool, one sim per process, 10,000 agents, 1985–2025. 20 more seeds run agent by agent to see where pairs diverge.
- EMOD's own `make_C` can't run on HIVsim's V yet: it has no ART holder and builds its table per node at EMOD's scale. That is step 11's adapter work.
- HIVsim over-transmits on this STIsim (analysis 11; since diagnosed as condoms being ignored by declared STIsim v4), so treat the numbers as orders of magnitude.

**Control:** Δ ≡ 0 (0 of 31,980 values differ). The divergence runs also match the runner's stored runs exactly. The runner tests pass (6/6).

## Variance reduction factor (var unpaired / var paired, 60 seeds)

| Output | 1995 | 2005 | 2010 | 2020 | 2024 (95% CI) |
|---|---|---|---|---|---|
| PLHIV | 34 | 5.5 | 5.1 | 3.1 | 2.6 (1.9–3.7) |
| New infections | 4.6 | 2.6 | 1.6 | 1.4 | 1.3 (1.0–1.8) |
| Incidence 15–49, women / men | 5.2 / 5.6 | 1.9 / 2.4 | 1.9 / 1.3 | 1.3 / 1.1 | 1.2 / 1.0 |
| HIV deaths | 1.6 | 1.3 | 1.6 | 1.0 | 1.4 (1.0–2.1) |
| ART coverage | – | 1.1 | 1.0 | 0.9 | 1.0 (0.8–1.3) |
| Share of deaths on ART | – | 0.9 | 1.0 | 1.0 | 1.0 |
| Transmission shares by stage (6) | – | 0.9–1.2 | 0.9–1.0 | 0.9–1.1 | 1.0–1.2 |

A dash means no ART yet in either arm. The empirical unpaired SD (each swap run against the next seed's base run) agrees with var(A)+var(B). `outputs/vrf.csv`, `figures/vrf.png`.

## Seeds for a 95% interval of ±10% of |Δ| (2024)

| Output | Paired | Unpaired |
|---|---|---|
| PLHIV (Δ −25%) | 7 | 18 |
| ART coverage (+0.086) | 16 | 15 |
| HIV deaths | 17 | 24 |
| Share of deaths on ART | 64 | 64 |
| Transmissions while on ART but not suppressed | 33 | 36 |
| Transmissions from undiagnosed (2020) | 119 | 111 |
| New infections (Δ −5%) | ~1,900 | ~2,550 |

`vrf.csv` also gives seeds for ±5% of the base mean.

## Where CRN breaks

- **First divergence:** the first differing starting agent appears between 1990.1 and 1991.5 (median 1991.1), when HIVsim's own ART starts (`n_art.csv`: 0 in 1990, 769 in 1991); the condensed cascade's ART starts only in 2004.
- **The first difference is not someone treated.** In seed 1 at ti 72, the same 5 people are infected in both arms, but 4 of them draw different acute durations. HIV's own distributions have the same call counts in both arms.
- **Cause: a Starsim agent-engine leak.** Engine draws are keyed by (step, number of the round of moves within the step). Every move made from code, such as ART starting someone, adds an untagged round: the base arm has 7 rounds per step from 1991, the swap 5–6. So one intervention's moves change the disease's later draws in that step for every agent, not just the people it touched (`outputs/engine_rounds_seed1.csv`).
- **How it spreads:** through infections, deaths and `StructuredSexual` rematching. The median share of starting agents differing in any state is 0.5% (1995), 2.3% (2000), 8% (2005), 14% (2010) and 19% (2024). By 2024, infection status differs for 13% and partner count for 11% (`outputs/divergence_by_year.csv`, `figures/divergence.png`, `figures/paired_trajectories.png`).
- **Why the VRF falls:** stocks keep the shared pre-1991 history, so PLHIV goes from 34 to 2.6. Flows lose it within about 10 years. Outputs defined by the cascade never shared any history, so they sit at VRF ≈ 1 from the start.

## Recommendation for the swap runs

1. Keep the shared seed per cell. It costs nothing, cuts seeds for stocks 2–3×, gives up to ~1.5× for flows, and the trivial-swap identity relies on it.
2. Size seeds by the unpaired variance for coverage, stage shares, deaths on ART and incidence. The planned 30 seeds (10 pairs × 3) is enough for the large Δs (PLHIV, coverage), marginal for deaths on ART and the smaller stage shares, and nowhere near enough for incidence Δs of a few percent. Report intervals, not point Δs.
3. Fix the engine leak in Starsim, then rerun this analysis (`run.py runs`, `diverge`, `rounds`, `analyse`). Some divergence is real and will remain. The EMOD′ steps have weaker CRN anyway, so they will need unpaired sizing.

## Upstream friction

- **Starsim engine leak** (`starsim/solvers/agents_moving.py`: `Moving.new_round()` ~714, `draw_ind()` ~734, `round_seed()` ~739). Untagged rounds come from `agents_stepping.py:35, 394, 430, 566, 633` and `agents_events.py:165, 333`. **Proposed fix:** give the engine's own passes fixed tags (e.g. `'step_state'`, `'transmission'`, `'transfers'`), and tag moves from code by the calling module plus its own counter, so one module's moves can't shift another module's draws. **Test:** two sims that differ only in an intervention that moves agents should give identical draws for agents it never touches.
- **Starsim `InfectionLog.to_df(ti=)`** (`starsim/diseases.py` ~1700) drops the `network` column in a step whose entries have no route. This crashed `common/outputs.py` `sexual_pairs()` on 15 of 60 seeds per hybrid (fixed on our side). **Proposed:** keep the column as NaN whenever the disease records routes.
- **Starsim `connectors.py:110` (`to_ti`)**, as in analysis 11: the condensed cascade's 2001 start blocks short runs, so the diagnostics run to 2005 (fixed since in starsim#1573).
- **Runner** (`common/swaps/runner.py:141`, `hybrid_id`): a config's model file outside `models/` is hashed by its path, not its contents, so editing it doesn't trigger a rerun (proposed: hash its contents). `build()` loads models in `set()` order, which changes between processes; `hivsim_cx.py` failed at random until it put HIVsim's folder on `sys.path` itself (proposed: load models in sorted order).

## What changed in code

- `common/outputs.py`: `yearly()` also gives year-end `pop`, `plhiv`, `diagnosed`, `on_art`, `suppressed` and `art_coverage`; `sexual_pairs()` handles a step's log with no `network` column. Both change hybrid IDs, so runs stored before these commits will rerun.

## Open

- EMOD's own `make_C` on HIVsim's V waits on step 11 (adapters for the path's hybrids).
- With its own ART, HIVsim has about 15 HIV deaths a year in 2020 and 37 in 2024 at 10,000 agents, 25–28% of them on ART; with the condensed cascade about 128 and 173, ~1% on ART (see analysis 13 on HIVsim's mortality on ART).
- Rerun this check after the Starsim engine fix.
