# 15: The five hybrids of the first path, through the runner — summary

**Answer:** yes. All five hybrids of HIVsim → EMOD′ (C, V, RN, D) build from the slot factories and run 1985→2025 through the runner: 25 of 25 runs (5 seeds) finished with no errors. The outputs analyzer gives every §2.2 output for each, and its counts agree with the disease module's own results. This is plumbing only. HIVsim over-transmits on rs-v4-art (analysis 11; diagnosed since as condoms ignored by declared STIsim v4) and has no reference sets, so no Δ is interpreted.

**Run times** (`outputs/run_times.csv`, 25 at once): HHHH 10.6–11.0 s; HHHE 15.2–15.8; HHEE 8.1–9.6; HEEE 10.1–11.3 (all 10,000 agents); EEEE 77.0–78.6 (1%, first run with the burn-in saved). Hybrids are named D, RN, V, C (H = HIVsim, E = EMOD′).

**Checks** (`run.py check`, seeds 1–2 per hybrid rerun in-process; `outputs/check.csv`, `counts_ratio.csv`, `counts_vs_disease.csv`):
- **Reproduced:** every stored table is reproduced exactly, including EEEE resumed from its stored 1985 burn-in.
- **No NaNs:** no output is NaN where its denominator is positive. The NaNs that occur are 0/0: `art_coverage` before any infection, `share_deaths_on_art` in years with no HIV deaths.
- **Counts against the disease's own results,** in agents, by year:
  - New infections and HIV deaths: ratio exactly 1 in every year, for every hybrid with an epidemic. EmodHIV's `new_infections` is sexual only, so for it the comparison is with the log's sexual transmissions (`tx_n`).
  - People infected at year end: median 1.00 (0.99–1.26). The outliers are early years with few infected, from the stock being recorded one step apart.
- **Take-off** (`outputs/epidemic_2020.csv`, PLHIV > 0 in 2020): HHHH 5/5, HHHE 5/5, HHEE 2/5, HEEE 3/5, EEEE 5/5. EMOD′'s seeding (2% of Risk HIGH) gives 0–5 seeds at 10,000 agents (analysis 12). Δs on those steps need many seeds, conditioning on take-off, or seeding as its own step (S).

**Bridges:** see `common/swaps/README.md`, "Hybrids of the first path".
- **New:** `adapters.EmodHostOnSti` (EMOD's cascade host on STIsim's demography); `EmodPeople.new_debuts`; HIVsim `components.for_cascade()` (ART holder, circumcision held at its pre-2010 level); HIVsim's D as modules when a foreign RN or V is present; HIVsim's RN without `ss.PrenatalNet` when V is foreign; EMOD's C summing flows over nodes and scaling by the sim on a foreign V or D. Factories detect a foreign slot with `own(ctx, slot)`.
- **Existing:** `EmodPeople`, `EmodHIVOnSti`, `EmodNetworkAt`, `EmodRiskAt`.
- **Runner:** `start()` resumes from `<hybrid>/<cell>.burnin.sim` when D's model has `BURNIN`. `hybrid_id` hashes config model files by contents, and `build()` loads models in sorted order (analysis 14's proposals).
- **Not built** (other paths): HIVsim's slots on EMOD′'s D, and EMOD's V on EMOD's D with HIVsim's RN (raises).

**Identity:**
- The all-HIVsim and all-EMOD′ builds still equal their `make_sim`.
- `build(C = EMOD′, rest HIVsim)` equals `swap.make_sim('condensed')` exactly (0 of 428 arrays), so step 1 is analysis 14's hybrid.
- 25 tests pass (`tests/`, `common/swaps/tests`, including `test_toy_burnin` and `test_path_hybrids`).

**Upstream friction:**
- HIVsim's demography passed as modules is not bit-identical to `sti.Sim`'s location form: pars and initial states are equal, but `ss.Pregnancy`'s distributions are registered twice, under `pars_demographics_pregnancy_*` traces (182 against 174 dists), so its draws differ from step 3. This weakens CRN between HHHE and HHEE. Proposed fix: register each distribution once, by the module's own trace; or let STIsim accept `demographics=['malawi', module]` (also proposed in analysis 12).

**Open:**
- Sharing one burn-in across hybrids needs slots swapped in at 1985 (§1.4); the EMOD′ → HIVsim path will need it.
- Low take-off with EMOD's V on HIVsim's D.
- `harness.py` imports `components` flat and picks up HIVsim's when HIVsim's folder is first on `sys.path` (harmless now, fragile): rename one of the two `components.py` or import by path.
