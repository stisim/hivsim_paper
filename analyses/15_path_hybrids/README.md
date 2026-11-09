# 15: The five hybrids of the first path, through the runner

Plumbing for step 11 of the first milestone (`docs/swap_design.md` §5): "the five hybrids of the path run, and the outputs analyzer gives every §2.2 output".

**Question:** do all five hybrids of HIVsim → EMOD′ (coarse slots in the order C, V, RN, D) build from the slot factories and run 1985→2025 through the runner (`common/swaps`), and does `common.outputs.Outputs` give sane values for each: no NaNs where the denominator isn't zero, and counts that agree with the disease module's own results?

**Hybrids** (runner names; the bridges each one uses are in `common/swaps/README.md`, "Hybrids of the first path"):

| # | D | RN | V | C | Scale |
|---|---|---|---|---|---|
| 0 | hivsim | hivsim | hivsim | hivsim | 10,000 agents |
| 1 | hivsim | hivsim | hivsim | emod | 10,000 agents |
| 2 | hivsim | hivsim | emod | emod | 10,000 agents |
| 3 | hivsim | emod | emod | emod | 10,000 agents |
| 4 | emod | emod | emod | emod | 1% (76,000 agents in 1985), resumed from its 1960.5–1985 burn-in |

**Not interpreted:** HIVsim on STIsim rs-v4-art over-transmits (analysis 11; diagnosed elsewhere) and has no reference parameter sets yet (every pid is hivsim_malawi experiment 03's best run). EMOD′ is parameter set 0. Δs between hybrids are not read scientifically here.

**What's run** (`run.py`, env `hp_rc3`, `PYTHONPATH=/tmp/sti-rs-v4-art`):
- `run.py runs`: the path through the runner, one pair (hivsim 0, emod 0), 5 seeds, stop 2025, into `raw_results/swaps/`. Run times from the manifest.
- `run.py check`: for seeds 1–2 of each hybrid, the run again in-process (from the stored burn-in for hybrid 4). It checks that the stored table is reproduced exactly, that the yearly outputs have no NaN where the denominator is positive, and compares the outputs' counts (new infections, HIV deaths, people infected at the end of each year) with the disease's own results (`sti.HIV`: `new_infections`, `new_deaths`, `infected`; `EmodHIV`: `new_infections`, `hiv_deaths`, `n_infected`), scaled by people per agent.

**What would count as done:** every run finishes without error; the check finds no unexplained NaN and the counts agree (ratio 1, up to a stock's one-step recording shift).
