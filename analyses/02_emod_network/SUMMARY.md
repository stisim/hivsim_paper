# 02: EMOD′ sexual network against EMOD — summary

**Answer:** the first EMOD′ network reproduces EMOD's closely. Run on the same parameter sets without HIV, nearly every statistic compared agrees with EMOD to within a few percent from 1961 to 2025: the share in each relationship type by sex and age, concurrency, formation rates, durations, coital acts, age differences, Risk mixing, and contacts per FSW and per client. The largest systematic difference is commercial sex: about 5% more people are in a commercial relationship at any time, and about 7% more commercial relationships are formed. Against the published model (with HIV), the two diverge after about 1990, by as much as HIV deaths change the population; most of all for marital relationships, many of which HIV ends.

## What matches

v4 against EMOD without HIV, on the same parameter sets (0, 100, 200). Figures are in `figures/`, numbers in `outputs/`. In the figures, lines are medians over the sets each model was run on (10 for v4 and EMOD with HIV, 3 for EMOD without HIV), so the dashed line isn't strictly comparable with the others; the numbers below use matched sets.

- **Population** (`population.png`): v4 is 0.3–1.2% larger than EMOD in 1980 and 0.8–2.1% larger in 2024. Numbers of FSW and clients (Risk HIGH) agree within run-to-run noise (−2% to +8% in single years). With HIV, EMOD's population is 12% smaller by 2024.
- **In a relationship, by type and sex** (`in_relationship.png`), as a share of those past debut: the mean v4 − EMOD difference over sets and sexes is between −1.2% and +1.7% for transitory, informal and marital, in 1980, 2000 and 2020. This holds across a wide range over parameter sets (transitory 3–7%, informal 4–30%). EMOD's marital peak (28% around 1978) and decline to 25% by 2005 are reproduced.
- **By age** (`by_age.png`, 1990 and 2020): shapes match, including marital peaking at 55–60 for men and 50–55 for women, and commercial sex confined to women aged 15–30. The 2020 gap at older ages is against EMOD *with* HIV.
- **Concurrency** (`concurrency.png`): the share with two or more partners and the mean number of partners are 1.8–3.6% above EMOD, mostly the commercial excess below.
- **Formation** (`formation.png`): relationships formed in 1970–2024 are +1.9% (transitory), +3.5% (informal) and +1.0% (marital), mostly v4's slightly larger population. The early spike, the marital decline to about 10 per 1000 a year by 1985, and the burn-in shape match. Rates are frozen from 1975.5, as in EMOD.
- **Queues** (`queues.png`): queue lengths match for every type, including EMOD's runaway commercial queue: after the burn-in, men's commercial entries outpace women's and the queue grows to about 1,800 entries per 1000 people past debut by 1995, all men. This is EMOD behaviour (frozen rates), reproduced.
- **Durations and acts** (`durations.png`): drawn durations are identical (mean 329, 727, 7300 and 7.0 days in both). Actual durations of ended relationships and coital acts per relationship match, including one act per commercial relationship.
- **Age differences** (`age_difference.png`): mean 3.4 / 4.6 / 4.8 / 6.7 years in v4 against 3.3 / 4.6 / 4.9 / 6.7 in EMOD without HIV (transitory / informal / marital / commercial); the distributions overlie each other.
- **Commercial contacts** (`commercial.png`), 1990–2024: 101.6 a year per FSW in v4 against 99.4 in EMOD without HIV, and 24.3 against 23.8 per client.
- **Risk mixing** (`outputs/mixing.csv.gz`): shares of each Risk pair are within about 1.5 percentage points for every type.

## What doesn't, and likely causes

1. **Commercial sex runs about 5% high**: +4% to +6% in a commercial relationship at every year compared, and +7.3% commercial formations in 1970–2024 (+4.6% to +12.7% over the 3 sets), against about +1.5% for population. Women's commercial queue-entry rate is frozen at its 1975.5 value, which depends on the FSW-to-client ratio then. v4 sees debut one step late (the risk module runs after the network), which shifts sex-work uptake around the freeze. Likeliest cause; not tested.
2. **Population is 1–2% larger** than EMOD's without HIV by 2024, growing. Candidates: pregnancy timing (EMOD's 280-day gestation starting within the step) and where in the step the yearly snapshot falls. Shares per person past debut are unaffected; counts are.
3. **With HIV, EMOD diverges from about 1990** (marital −20% and concurrency −8% by 2020, relative to v4), as expected with no HIV deaths in the harness. Up to the late 1980s v4 is within about 1% of EMOD with HIV for every type.
4. **Small simplifications** not visible in these statistics: when a woman queued twice for one type is matched, the matched entry is removed rather than EMOD's most recent one; EMOD's 50,000-candidate cap isn't applied (women's queues stay short, so it never binds). See `models/emod/v4/README.md`.

## Caveats

- One v4 run per set against one EMOD run per set (3 without HIV). Single-year differences of a few percent are within noise; the systematic ones (commercial sex, population) appear in all 3 sets.
- The harness demography is minimal and copies EMOD's 1000-times-too-low background mortality. It is not EMOD′'s demography component.
- Snapshots are taken in the step of EMOD's report times (x.92), but not at the same point within the step.
- v4 uses its own random numbers (EMOD's within-step queues need sequential draws per person), so it isn't safe for common random numbers.

## Next steps

- Find the commercial excess: compare EMOD's frozen rates with `net.rates`, and detect debut within the network step.
- Add HIV (EMOD′ transmission from `edges.acts` and `edges.condoms`, and natural history), then rerun against EMOD with HIV.
- Replace the harness demography with the EMOD′ demography component.
