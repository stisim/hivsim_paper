# 08: The condensed cascade in Optima′

**Question:** the condensed cascade (`common/cascade/`) acts on agents; Optima′ is compartmental. `common/cascade/rates.py` turns the cascade's table into rates between Optima′'s care states (undiagnosed, diagnosed, in care, on ART suppressive or not, lost), by CD4 stage, population and year (the mapping is in `common/cascade/README.md`). Two questions:

1. **Is the rate approximation consistent with the cascade it stands for?** In EMOD′ running the condensed cascade with the same table, measure the rates at which people with HIV move between those care states each year, and compare them with what `rates.py` gives for the same people.
2. **What does EMOD's cascade do to Optima′?** Optima′ with its own cascade against Optima′ with EMOD's cascade as rates, on the common targets and on the cascade.

**What's run:**
- (i) `harness.make_sim(0, seed, cascade=True, condensed=True)`: EMOD′, parameter set 0 (whose table is `common/cascade/tables/malawi_emod_set000_*.csv`), 3 seeds (seeded from its run number, as analysis 06), to 2025, with Starsim v4 from `~/starsim` (commit ebd52c1a, with PR #1538) on `PYTHONPATH`. An analyzer records, each month, each person with HIV's care state (mapped as below), CD4 stage (Optima's bands, from EMOD's CD4; acute by EMOD's stage), and group (sex and 5-year age band); person-time and moves between care states; entries into CD4 < 200 by the undiagnosed and diagnosed (where `rates.py` puts first symptoms); pregnancies at 12 weeks per woman, by group (the pregnancy rate `rates.py` is given); and the hidden shares `rates.py` approximates: the share in the testing loop of everyone past debut, of the undiagnosed and of the diagnosed, and the share of those dropped out of ART who are out rather than lost for good.
  - Care states from the condensed cascade's places: on ART → svl or usvl by arm; staging, pre-ART, linking → care; out, lost → lost; else diagnosed → dx; else undx.
  - Expected moves: person-time in each state, stage and group, times the rate of `rates.py` for that month, plus entries into CD4 < 200 times its probabilities on first symptoms. Compared with the moves seen, by year and by pair of states.
- (ii) Optima′ (all components, dt = 0.02, to 2025): with its own cascade; with EMOD's (set 0's table, with its VMMC and PrEP flows); and with EMOD's without the flows. Measured with `models/optima/reference/optima_targets.measure()` on the common targets, and the cascade over time (diagnosed, in care, on ART, suppressed, lost, as shares of PLHIV).
- `run.py` runs and writes `outputs/`; `run.py --plot` draws `figures/`. Optima′ runs in the `hivsim_paper` env (Starsim 4.0.0rc1); EMOD′ needs Starsim head.

**Expected differences in (i)**, by construction of the approximation:
- Linkage and the start of ART after staging are taken as at once in `rates.py` (EMOD: 1–120 days in linking), so undx → svl/usvl is larger, and undx → care and care → svl/usvl smaller, than seen; the totals out of undx should agree better than the pairs.
- First symptoms are tied to entering CD4 < 200; in EMOD they come a Weibull time before death.
- The pre-ART check (every 182 days) is a rate; the testing loop and the dropout trap are shares from small ODEs, not people.

**Expected differences in (ii):** Optima's ART numbers are forced (`numtx`), EMOD's cascade is not, so the number on ART becomes an output; EMOD's dropout trap holds many people with HIV out of care; EMOD's arms fix suppression at ART start (no failure or viral load monitoring); and VMMC and PrEP are added, which Optima doesn't have.
