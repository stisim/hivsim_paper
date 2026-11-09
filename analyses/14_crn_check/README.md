# 14: Common random numbers on step 1 (C)

Step 10 of the first milestone (`docs/swap_design.md` §5; §3, "CRN check"). The runner gives every hybrid in a cell the same seed (§4). This checks whether that pairing pays off.

**Question:** does sharing the seed between a model and its one-slot swap lower the variance of the swap's effect Δ enough to matter? And where does pairing stop working?

**Swap:** step 1 of the path, C (care cascade). HIVsim (`D, RN, V, C = hivsim`) against HIVsim with EMOD's condensed cascade as C (`C = hivsim_cx`, `hivsim_cx.py`). That is `swap.make_sim('condensed')`: EMOD's Malawi table set 0 summed over nodes, on `StisimHost`, with `sti.ART` holding ART's states and circumcision held at its pre-2010 level. EMOD's own `make_C` (`models/emod/v4/components.py`) isn't usable on HIVsim's V yet. It builds the table from EMOD's inputs per node at EMOD's scale, with no ART holder, which is step 11's adapter work. So the C swap here is the condensed variant that HIVsim already has.

**Control:** a trivial swap, C from `hivsim_b`, which is HIVsim's own factory file under another name. With shared seeds Δ must be exactly 0 for every output, year and seed.

**Caveat:** HIVsim on STIsim rs-v4-art over-transmits (analysis 11; under investigation elsewhere). Absolute numbers and the size of Δ are not meaningful. How the variance behaves is.

**What's run** (`run.py`, env `hp_rc3`, `PYTHONPATH=/tmp/sti-rs-v4-art`):
- `run.py runs`: the three hybrids through the runner (`common/swaps`), seeds 1–60 (one pair, 60 seeds), 10,000 agents, 1985 to 2025. Stored in `raw_results/swaps/`, one sim per process.
  - `common.outputs.yearly()` now also gives the year-end stocks (PLHIV, diagnosed, on ART, suppressed, ART coverage). The tables before it had flows and shares only.
- `run.py diverge`: seeds 1–20, the three hybrids run one after another in one process per seed. An extra analyzer records, at every step and for each of the 10,000 starting agents, whether they are alive, infected, diagnosed and on ART, and how many sexual partners they have. Paired runs are compared agent by agent: at each step, the share of starting agents whose state differs between base and swap, and between base and control. The analyzer draws no random numbers. As a check, the base and control runs here must equal the runner's stored runs exactly.
- `run.py analyse` writes the tables and figures.

**Outputs and years:** ART coverage, PLHIV, new infections (all ages), HIV deaths, share of HIV deaths on ART, incidence 15–49 by sex, and the shares of sexual transmissions by the source's cascade stage (undiagnosed, diagnosed and never on ART, ART under 6 months, ART suppressed, ART not suppressed, interrupted). Years 1995, 2000, 2005, 2010, 2015, 2020 and 2024.

**Measures:** per output and year, over the 60 seeds:
- Paired: Δᵢ = Y_swap(seed i) − Y_base(seed i), and var_p = var(Δᵢ).
- Unpaired: var_u = var(Y_swap) + var(Y_base), the variance of the difference of independent runs. It is also checked empirically, with Δᵢ = Y_swap(seed i) − Y_base(seed i+1).
- Variance reduction factor VRF = var_u / var_p, with a 95% bootstrap interval over seeds.
- Seeds needed for a 95% interval on mean Δ of half-width h, n = var·(1.96/h)², paired and unpaired, for two choices of h: 10% of |mean Δ|, and 5% of the base mean.

**What would count as mattering:** VRF ≥ 2 halves the seeds needed. VRF near 1 means pairing buys nothing, so the swap runs should rely on more seeds (§3, consequences). The control must give Δ ≡ 0. The cascades start in 1990 (EMOD's table: HCT uptake, symptomatic testing) and 2001 or later (HIVsim's own: ART from its numbers on ART, testing from `test_start`). So under perfect CRN, divergence should start no earlier than the first step where C acts differently. Circumcision differs from 1985, since the condensed variant holds it at its pre-2010 level, so that is a possible earlier source.
