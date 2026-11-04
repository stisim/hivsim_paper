# The condensed HIV care cascade

One description of an HIV care cascade (testing, staging and linkage, pre-ART, ART, dropout and return, PMTCT, VMMC and PrEP) as a table of parameters, which runs in any model that implements its interface, so the cascade can be swapped between models.

| File | What it is |
|---|---|
| `params.py` | The table (one row per pathway parameter: constants, sigmoids, steps, linear series, delays) and the flows table (VMMC and PrEP numbers by group and period); `Params` reads them |
| `tables/` | Tables for Malawi: `malawi_emod_set000_{params,flows}.csv`, EMOD's campaign for parameter set 0 (`models/emod/v4/extract_cascade.py`) |
| `interface.py` | `CascadeHost`: what the cascade needs from an agent-based model (people, the events of a step, HIV) |
| `cascade.py` | The cascade for agents, as Starsim deliveries acting through a `CascadeHost`; one place in care per person |
| `rates.py` | The cascade for compartmental models: the table as rates between care states, by CD4 band, group and year |

EMOD′ runs `cascade.py` (`models/emod/v4/cascade_host.py`); Optima′ runs `rates.py` (`optima_prime.make_sim(cascade=Params(...))`).

## The compartmental mapping (`rates.py`)

A compartmental model has counts by care state, CD4 band and group, not people with a place in care, a test schedule and a record of their lowest CD4. `CompartmentalCascade(params, groups, years)` gives:

- `rates[(src, dst)]`: rates per year between care states, bands × groups × years, within a CD4 band (moves between bands stay the host's natural history);
- `on_progression[(src, dst)]`: probabilities on entering the symptomatic band (first symptoms);
- `births[state]`: where an infected birth goes (the test at 6 weeks);
- `mtct_averted`: the share of mother-to-child transmission averted, for mothers not in care, by ANC testing and PMTCT;
- `flow_rows('vmmc' | 'prep', year)`: the numbers a year, by sex and ages.

Care states: **undx**; **dx**, diagnosed but not in care (the cascade's none, uptake or testing); **care** (staging, pre-ART); **svl** and **usvl**, on ART in the suppressive or non-suppressive arm; **lost**, dropped out of ART (the cascade's out and lost). These are Optima's states. The host gives: groups (sex, ages), the CD4 bands and which one starts symptoms, pregnancies per woman by group and year, and optionally the share of each group past debut (default: everyone from 15).

### Exact (given the care states)

- Testing in the routine loop (exponential intervals): a rate of 1/`test_interval` for those in the loop; entering the loop from uptake at `p_uptake`/`uptake_interval`.
- What follows a positive test: immediate ART (by sex), staging lost, linkage to ART if eligible, pre-ART if not, else back to dx; ANC tests staged with `p_stage`. Products of the table's probabilities.
- ART arms: suppressive with `p_suppressive` by sex and age; the arm is kept (no failure or regain, as in the cascade).
- Dropout from ART (exponential `dropout_delay`): lost for good, out, restart (with the arm drawn again), or back to dx.
- The dropout trap: with `p_restart` 0, those out never restart; they return only to dx, at (1 − `p_lost`)(1 − `p_out`)/`dropout_delay`.
- The Option B block (`blocked_if_staged`) and the PMTCT options: as expected efficacies by year.
- Event-triggered testing at an event the host has as a transition: first symptoms split the progression into the symptomatic band; the infant test splits infected births.

### Approximate

- **Hidden places in care.** The loop share (who in undx or dx is in the testing loop rather than waiting in uptake) and the share of lost who are out rather than lost for good are not care states. `rates.py` follows them with small ODEs: everyone past debut (into the loop from uptake, and at debut with `p_debut`; turnover `debut_turnover`), the undiagnosed (entering with everyone's share, leaving by tests and at `other_exit`, the pool steady in size), the diagnosed (entering in uptake), and the dropped out (from the first year anyone can start ART). Per person this is exact; the shares are of the whole pool, not of each band and group.
- **Delays as rates or skipped.** The pre-ART check (every 182 days) is a rate; staging (1 day) and the ART start delay (1–120 days in linking) are taken as at once, so linkage goes straight to ART and `care` holds only pre-ART.
- **Eligibility by CD4 band.** The rule (WHO stage, CD4, pregnancy, children's ages) is applied to a band as the share of it eligible, with people spread evenly in time over the band and WHO stages from EMOD's natural history (√CD4 linear over the prognosis, WHO stages Weibull shares of it: `Eligibility`). Pregnancy is a share of women (pregnancies × 280 days); ages are spread evenly over each group. Staging uses the current band: the lowest CD4 a person has had (which matters after ART dropout) has no compartmental counterpart.
- **First symptoms** are the entry into the symptomatic band (Optima's CD4 < 200), once per infection; in EMOD they come a Weibull time before death, and again after stopping ART.
- **Flows.** VMMC numbers are spread evenly over the men of the ages and HIV-negative (regions added up); PrEP numbers become coverage (starts × days protected over the women of the ages). Risk groups (EMOD's `risk:HIGH|MEDIUM`) can't be told apart.
- **PMTCT** applies to mothers in undx or dx (those the ANC test reaches); mothers in care or lost get none, as in the cascade; those on ART get the host's ART effect.

### Not representable

- Per-person schedules and memory: the next test date, the lowest CD4, the PMTCT window of each mother, repeated first symptoms after ART.
- Extra event-triggered testing pathways (`testing=` in `cascade.py`) whose trigger isn't a transition of the host.
- Cascade states that are levels, not rates (the number on ART, proportions held fixed): the table has none; Optima's own (`CascadeTargets`) are off when a table is used.

### Other cascades

STIsim/HIVsim's cascade is already in rates or probabilities per agent: `HIVTest`'s annual testing probability p is a rate −ln(1 − p); `HIVCare`'s `rate_loss`, `rate_loss_art` and `rate_return` are rates between dx/care, svl/usvl and lost; `ART`'s `p_effective_art` is the arm. Its coverage targets (`ART(coverage=...)`) are levels, like Optima's `numtx`, and map to a rule on levels (Optima′'s `CascadeTargets`), not to rates.

Analysis 08 (`analyses/08_optima_condensed_cascade/`) checks these rates against EMOD′ running the same table, and runs Optima′ on them.
