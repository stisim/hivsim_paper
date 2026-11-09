# 11: HIVsim on rc3, slot factories, EMOD′'s 1985 burn-in, and HIVsim's scale — summary

**Answers:**
1. **HIVsim on rc3** runs, short and full, with its own care and with the condensed cascade. The only code change was `ss.MaternalNet` → `ss.PrenatalNet`. Its outputs are **not** the same in distribution as before. With experiment 03's best set, prevalence 15–49 keeps rising, to 19% in 2020 against 8% before (×1.4–2.4 over 2000–2020). Incidence is ×2.2 and the share of infected on ART ×0.58 (the same numbers on ART over more people infected). The population is 4% lower, which is mostly the one-step shift in when rc3 records a stock (`hivsim_rc3_vs_before_by_measure.csv`, `figures/hivsim_outputs.png`). This is a change in STIsim, not in the migration (below). The set needs recalibrating, or STIsim's port needs checking, before HIVsim's numbers are used.
2. **Slot factories:** identical. `harness.make_sim`, `malawi.make_sim` and `swap.make_sim`, rebuilt on `components.py`, give the same runs as before (at 5df2972), same seed: 0 arrays differ out of 55–417 per case (EMOD′ network, HIV, and the cascade interpreted, compiled and condensed to 2010; HIVsim Malawi, own and condensed). The swap runner's `build(all HIVsim)` also equals `swap.make_sim('own')` (0 of 296), and its own tests pass on the real factories (`common/swaps/tests`, 6/6). `pytest tests` gives 10/10. A negative control (seed 2 against 3) finds 196 differing arrays.
3. **EMOD′'s burn-in:** resumed equals continuous **exactly** (0 of 240 arrays differ, for each of 3 seeds, to 2025), once `EmodNetwork` saves Numba's generator state with the sim. Before that fix, a resume in the same process was exact, but one in a new process differed from the first resumed step (ti 294, 1985.0) in 212–217 of 240 arrays. Starsim's save and resume needed no fix: the state that was lost was EMOD′'s own, Numba's global stream (the risk in §5). `tests/test_resume.py` checks it, in a new process.
4. **HIVsim's scale:** 1% (76,205 agents) agrees with 10,000 agents within the spread over seeds. Medians are within 3.5% for every measure and year, and 146 of 168 points (87%) have |z| ≤ 2. That is more beyond 2 than the 5% the README set. The cause is small systematic differences that 20 + 10 seeds resolve: prevalence 2–3.4% lower at 1% in 2004–2010 (z down to −3.0), and the population 0.5% lower throughout. Both are well inside the 10,000-agent runs' own 5–95% range, which is about ±6% for prevalence (`hivsim_1pct_vs_10k_by_measure.csv`). The number on ART is set by the target, so its z is not meaningful (ratio 1.00–1.03). So agent number moving with D is fine for the swaps. A path that needs differences below about 3% in prevalence should make the agent number its own step (§5 risks).

**Run times** (`burnin_times.csv`, 3 seeds, 3 at once; `hivsim_runtime.csv`, 20 at once):

| Run | Time |
|---|---|
| EMOD′ (condensed cascade), 1960.5–2025 continuous | 72–74 s |
| EMOD′ burn-in, 1960.5–1985 | 10.1–10.2 s |
| Save at 1985 (13–14 MB) / load | 3.2 s / 0.3–0.4 s |
| EMOD′ resumed, 1985–2025 | 58–60 s |
| HIVsim, 10,000 agents, 1985–2022: before / rc3 | 6.4 s / 8.7 s |
| HIVsim, 1%, 1985–2022 | 53 s |

The burn-in is cheap (14% of a run), since HIV is seeded only in 1983 and the cascade starts later. Saving it saves 10 s per hybrid. EMOD′ from 1985 costs about as much as HIVsim at 1%.

**What changed in code:**
- `models/emod/v4/components.py`, `models/hivsim/components.py`: one factory per slot of §1.1 (D, R, N, T, H, S, C, M), and the runner's coarse RN (R + N) and V (T + H + S), in the runner's interface (`HOST`, `SCALE`; `pars(pid)`; `make_<slot>(pars, ctx)` returning sim keywords). HIVsim also has `SIM = sti.Sim`, since its D gives a location for STIsim to build. The runner now makes the sim class of D's model (`common/swaps/runner.py`, one line). Where one module fills several slots, the first slot builds it and the others give nothing: T and S are inside `EmodHIV` and `sti.HIV`, and HIVsim's R is inside `StructuredSexual`. EMOD's C adds `Snapshots` and `EmodReport`, and HIVsim's V adds `ss.infection_log()`, which doesn't change the run. HIVsim's ART is made in one place (`components.make_art()`), for the move to the separate ART module.
- `models/emod/v4/cascade.py` → `emod_cascade.py`, and `models/hivsim/cascade_host.py` → `stisim_host.py`, so that no two modules on the runner's path share a name.
- `models/emod/v4/network.py`: `EmodNetwork.__getstate__/__setstate__` carry Numba's generator state.
- `common/identity.py`: every result array and every analyzer's or logging module's table, for exact comparisons.

**Why the HIVsim shift is STIsim's.** Run at the same seeds on rc3, the STIsim commits that work with rc3 (c40c58c, f395c9d and c1127e3: ART as its own module, ART's own random numbers, and the host's methods) give identical medians. So the shift isn't in the ART work. It is in the base the port was made on (v4-declared, #610, with the rc3 updates), compared with the port the "before" env has (rs-art-split, c815b0a, on rc1). The parameter values are the same (both modules' `pars`, diffed). The network has about as many edges and acts. Turning off the stage multipliers (`rel_trans_acute` and `rel_trans_falling` = 1) leaves rc3 at 10.5% in 2000 against 7.2%, so it isn't stage infectiousness either. New infections per person infected are about 25–30% higher from 1990, before ART, and that compounds. The next thing to check is transmission per act along an edge in declared `sti.HIV` (per-act `beta` against acts per step), with a fixed network on both versions.

**Upstream friction (Starsim rc3):**
- `starsim/connectors.py:110` (`to_ti`) raises if an intervention's `start` or round is after the sim's `stop`. For example, HIVsim's condensed cascade (`cx_anc_test`, 2001) can't be run to 2000, and EMOD′'s (`L17_yes`, 2004; `cx_hct_init`, 1990) can't be run to 2003. A short check run then has to go to 2010. **Proposed fix:** a start or round after the stop is allowed, and the module simply never acts (with a note at `verbose`), while a start before the sim's start still raises.
- Not a Starsim bug, but worth a line in `sim.save()`'s docstring: state that lives outside the sim (Numba's or NumPy's global generators) isn't saved, so a checkpoint resumed in a new process can differ.

**Open:**
- HIVsim on rs-v4-art needs either recalibration (experiment 03's set gives 19% prevalence) or a check of STIsim's transmission (above). The swaps can proceed (identity doesn't depend on fit), but HIVsim's endpoint numbers can't be used until then.
- HIVsim's `pars(pid)` gives experiment 03's best set for every pid until `models/hivsim/reference/` has the sets (§5 step 2).
- The runner using D's sim class (`sti.Sim` for HIVsim's D) has to be checked when an EMOD′ V or C runs on HIVsim's D (adapters, steps 6–7). `sti.Sim` also processes diseases and connectors.
- Numba's stream is still global: two EMOD′ sims stepped in one process share it. The adapters' `own_stream` (analysis 12) avoids that, and the runner runs one sim per process.
