# HIVsim reference

HIVsim Malawi is calibrated in [hivsim_malawi](https://github.com/stisim/hivsim_malawi). When it has a finalized parameter set, the model and that parameter set are copied here (with the hivsim_malawi commit they came from), and are not tuned further in this repo. HIVsim is already a Starsim v4 model, so it has no separate `v4/` folder.

| File | What it is |
|---|---|
| `pars.json` | **Provisional** (9 Nov 2026): hivsim_malawi experiment 03's best run (commit 803b2f7), with its source and caveats. Read by `swap.BEST` and `components.pars()`. It needs STIsim with the condom fix (stisim#612 at 9c1dc30 or later); with it, prevalence 15–49 is ~8% in 2020. Chosen so the swap path can proceed (endpoint identity doesn't depend on fit); calibration resumes after on-ART mortality (stisim#613), and this set will be replaced |
