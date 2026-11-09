# 21: Seeding options on the path HIVsim → EMOD′-fixed

**Question:** on the path HIVsim → `emod_fixed` (analysis 16's addendum: order C, V, RN, D, 200 cells), the hybrids with EMOD′'s HIV on HIVsim's demography fail to take off in some cells (HHEE 144/200, HEEE 164/200). EMOD's seeding (2% of Risk HIGH at the first step on or after 1983) gives 0–5 seeds at 10,000 agents (analysis 12). How do take-off and the 2020 headline Δs (new infections; Figure 3A shares: undiagnosed, diagnosed never treated, on ART, interrupted) change under other seeding rules? This measures options for `swap_design` D11; it doesn't choose one (that is Robyn's decision).

**Options** (each a seeding option for EmodHIV, off by default, so every existing hybrid is bit-identical):

| Option | Rule | Hybrids run |
|---|---|---|
| 0 (addendum) | EMOD's: each HIGH agent seeded with probability 2% (`seed_coverage`) at the first step | the 5 of the path, rerun at this commit (bit-identity with analysis 16's stored tables checked) |
| (a) scaled | The same expected number of seeds per capita as EMOD′-fixed at its own scale: r = E[seeds]/alive at EMOD′-fixed's seeding step (measured here over 50 seeds of its burn-in), and each HIGH agent seeded with probability r·alive/HIGH (at most 1) | HHEE-a, HEEE-a |
| (b) S slot | Seeding as its own slot S. S = HIVsim: HIVsim's rule (`init_prev_data` × `rel_init_prev`, by risk group, sex and sex work, among the sexually active) applied to EmodHIV at its seeding step, in place of EMOD's. On HIVsim's network (HHEE) it reads StructuredSexual directly; on EMOD′'s network (HEEE) it is HIVsim's mean seeding probability by EMOD's Risk class and sex (measured on HIVsim's 1985 population, with Risk read as the adapters do: HIGH = FSW or client, LOW = risk group 0, MEDIUM otherwise), among those past debut. The seeds' time since infection is EmodHIV's (fast-forwarded as EMOD's), since the natural history is V's | HHEE-S, HEEE-S |
| (c) more seeds | 5× EMOD's share (10% of HIGH) | HHEE-c, HEEE-c, and EEEE-c |

With (b), the path gets a seeding step: C, V, S, RN, D (S right after V; HHEE-S between HHHE and HHEE), or C, V, RN, S, D (HEEE-S between HHEE-S and HEEE). "S with D" is the second without its S step: HEEE-S → EEEE as one step. "S with V" is the addendum's path.

With (a) and (c) the path is C, V, RN, D with HHEE and HEEE replaced; EEEE stays EMOD′-fixed, so under (c) the D step also takes the seeding back to 1×. EEEE-c (EMOD′-fixed seeding 10% of HIGH in 1983) shows how much that matters at EMOD′'s scale.

(a) and (c) use the same random draws as EMOD's rule (a uniform per HIGH agent, compared with a different threshold), so their seeds are a superset or subset of option 0's in each cell.

**Cells:** the same 200 (`p0_s1`–`p0_s200`), pair (hivsim 0, emod_fixed 0), stop 2025, tier 1 (uncompensated).

**Outputs:** per option, take-off (prevalence all ages ≥ 0.5% at the end of 2020) per hybrid with Wilson intervals, Δp per step, and the 2020 headline Δs per step: (all) all cells where defined, (both) both arms took off, and (path) cells where every hybrid of that path took off, so that the step Δs add up over a common set of cells. The same three summaries for the addendum's results (option 0), next to the take-off probabilities. Compared with the addendum's table.

**Done when:** all runs finish without error; option 0 reproduces analysis 16's stored tables exactly; tests in `tests/` and `common/swaps/tests` pass; take-off and the headline table for every option; a SUMMARY saying which conclusions hold under every option.
