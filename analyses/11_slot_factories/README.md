# 11: HIVsim on rc3, slot factories, EMOD′'s 1985 burn-in, and HIVsim's scale

Steps 3 and 4 of the first milestone (`docs/swap_design.md` §5), with HIVsim moved to Starsim 4.0.0rc3 first.

**Questions:**
1. **HIVsim on rc3.** With Starsim 4.0.0rc3 and STIsim's v4 port (branch `rs-v4-art`: declared `sti.HIV`, coded `sti.ART`, the `CascadeHost` methods), does HIVsim Malawi (`models/hivsim/`) run, and do its headline outputs agree in distribution with the version before (env `hivsim_paper`: Starsim 4.0.0rc1, its STIsim; `~/hivsim_paper` at main)? Exactness isn't expected: STIsim's HIV is now declared, and rc3's draws are keyed differently. The migration itself is one rename (`ss.MaternalNet` → `ss.PrenatalNet`, which is what MaternalNet was: prenatal only).
2. **Slot factories.** `models/emod/v4/components.py` and `models/hivsim/components.py` give one factory per slot of §1.1 (D, R, N, T+H, S, C, M; T, H and S are one module in both models, so one factory builds them and the others say so). Each model's `make_sim` is rebuilt on them. Test: `make_sim` before and after gives identical results (every array of `sim.results`, and the analyzers' tables), same seed, `tests/test_components.py`.
3. **EMOD′'s burn-in.** EMOD′ (full model, condensed cascade) from 1960.5 to 1985, saved, loaded and run on to the end, against one run from 1960.5 to the end, same seed: every result array must be identical (§1.4). Run times: the burn-in, the run on from 1985, the continuous run, saving and loading.
4. **HIVsim's scale.** HIVsim at 10,000 agents against HIVsim at 1% scale (Malawi's 1985 population × 0.01, about 80,000 agents), several seeds each: do outputs that don't depend on the number of agents (prevalence, incidence, share on ART, and the scaled counts) agree in distribution (§1.4, Decided 9 Nov)?

**What's run** (`run.py`, env `hp_rc3` with `PYTHONPATH=<STIsim rs-v4-art>`):
- `run.py hivsim --tag rc3` and, in env `hivsim_paper` with `--repo ~/hivsim_paper`, `--tag before`: `swap.make_sim('own')` (experiment 03's best set), 10,000 agents, seeds 1–20, to 2022; yearly outputs to `outputs/hivsim_<tag>.csv`.
- `run.py burnin`: EMOD′ set 0, seeds 1–3, to 2025: continuous; and burn-in to 1985, `sim.save()`, `ss.load()`, run on. Identity per result array to `outputs/burnin_identity.csv`, times to `outputs/burnin_times.csv`.
- `run.py scale`: HIVsim (`own`), seeds 1–20 at 10,000 agents (the rc3 runs above) and 1–10 at 1%, to `outputs/hivsim_scale1pct.csv`.
- `run.py compare` writes the tables; `run.py --plot` the figures.
- The identity tests: `pytest tests/test_components.py`.

**How results compare:** questions 2 and 3 exactly (any difference is a failure). Questions 1 and 4 in distribution, per output and year (2000–2020): z = (median₁ − median₂)/SE, SE from the two sets of runs (as analysis 06). Agreement: |z| > 2 in no more than about 5% of points and no systematic shift. For question 1 a shift is possible and would be STIsim's change rather than the migration's; it is reported, not fixed here.
