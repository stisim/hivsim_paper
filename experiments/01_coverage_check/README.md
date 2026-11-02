# Exp 01 — Can HIVsim Malawi reach the targets?

**Question.** Over plausible ranges for its uncertain parameters, does HIVsim Malawi (`models/hivsim/malawi.py`) produce runs that bracket the agreed targets (`results/targets.csv`, `docs/targets.md`)? This is a check of the model and the ranges before any fitting. If some target falls outside every run, the ranges are wrong or the model is missing something, and fitting won't fix that. An exact match isn't the aim yet.

**Plan.**
- 200 parameter sets from a Latin hypercube over the ranges below, one run each, 10,000 agents, 1985-2021. Also 10 runs (different seeds) of one central set, to show how much of the spread is random noise, and the run time.
- Ranges (uniform unless marked log; the defaults are STIsim's):

  | Parameter | Range | Default | What it does |
  |---|---|---|---|
  | `hiv.beta_m2f` | 0.002-0.05 (log) | 0.05 | per-act transmission, male to female (female to male is 0.5 of it) |
  | `hiv.rel_init_prev` | 0.25-4 (log) | 1 | size of the epidemic in 1985 |
  | `nw.prop_f0` | 0.6-0.95 | 0.85 | share of women in the lowest risk group |
  | `nw.prop_m0` | 0.5-0.9 | 0.8 | share of men in the lowest risk group |
  | `nw.f1_conc` | 0-0.2 | 0.01 | concurrency, women in the middle risk group |
  | `nw.m1_conc` | 0.05-0.4 | 0.2 | concurrency, men in the middle risk group |
  | `nw.fsw_shares` | 0.01-0.08 | 0.05 | share of women who ever sell sex |
  | `nw.client_shares` | 0.05-0.25 | 0.12 | share of men who are ever clients |
  | `test_scale` | 0.25-4 (log) | 1 | multiplies the probabilities of general and low-CD4 testing |

- For each run, save the model's value for every target and check in `results/targets.csv`, and the share of the ART input reached, to `outputs/results.jsonl`.
- Figures: each target's data with its interval against the spread of the 200 runs (all runs, and the 5-95% range), and the 10 seeds of the central set.

**Success criteria.** Every primary target (prevalence by sex and age, five surveys) and the UNAIDS series lie inside the range of the runs, preferably inside the 5-95% band, with the shape across ages reached by some runs. The ART input should be met in most runs. If a target is outside every run, the result says which, and what it implies (a range to widen, a missing input such as condom use or circumcision, or a structural limit). Noise: the seed-to-seed spread should be small compared with the spread across parameter sets; if it isn't, 10,000 agents is too few.
