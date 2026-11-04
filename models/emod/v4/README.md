# EMOD′ (EMOD-HIV re-expressed in Starsim v4)

The whole published model: demography, the sexual network, Risk and sex work, HIV transmission and natural history, and the campaign's care cascade (testing, staging, pre-ART, ART, dropout), VMMC, PrEP and PMTCT.

| File | What it is |
|---|---|
| `inputs.py` | `load(pid)`: the parameters for one of Daniel's parameter sets, read from the EMOD inputs that `reference/run_emod.py` builds for it, so his quirks come with them (e.g. the Risk overlay's node order) |
| `network.py` | `EmodNetwork`, an `ss.SexualNet`: EMOD's pair-formation algorithm and its four relationship types (transitory, informal, marital, commercial): queues with per-person entry rates set during the 15-year burn-in and then frozen, matching on the joint age tables and Risk weights, the extra-relational flag chain and randomly rounded maxima, Weibull durations in whole steps, and Poisson coital acts with dilution and a forced first act. The module docstring lists the steps and where each is in the EMOD source |
| `risk.py` | `EmodRisk`, an intervention: `Risk` drawn by node at first sight; the campaign's sex work (at debut from 1961, will ever be FSW or client by node; uniform delay to Risk HIGH, then Weibull or uniform time to MEDIUM); and the STI co-infection flag, given at debut by Risk |
| `hiv.py` | `EmodHIV`, a disease written as code: per-act transmission over the network's acts with EMOD's multipliers (stage, STI flag, condoms, circumcision, male-to-female by age), acts taken in random order per person; transmission at birth; the age-dependent prognosis, acute and AIDS stages and deterministic CD4; traditional circumcision from 1975; and the seeding at SeedYr. The docstring gives EMOD's order within a step |
| `demography.py` | `EmodDemography`: EMOD's initial population (node sizes at 1% scale, ages from its CDF, pregnant at the start), pregnancies by age and year with 280-day gestation, background deaths by sex, age and year with the mortality tables read 1000 times too low as EMOD reads them (`fix_mortality=True` reads them as labelled), 3 nodes, no migration, and ages in single-precision days as EMOD keeps them |
| `cascade.py` | `EmodCascade`: the rest of Daniel's campaign, read from the campaign JSON and run with EMOD's semantics (listeners, `Disqualifying_Properties` checked when given and at each update, `New_Property_Value` at the first update, chains and short delays resolved within the step, muxers, NChooser targets scaled and truncated). EMOD's quirks come out of the file as they do in EMOD: the dropout trap, the ignored ANC uptake table, and the PMTCT combination being disqualified by ART staging. The docstring lists the semantics |
| `compile_campaign.py` | The same campaign compiled into native Starsim deliveries (`ss.Delivery` with conditions, `coverage=`, `delay=`, `at=`, `deliver=`), with small coded products for what has no native form; `harness.make_sim(..., cascade=True, compiled=True)`. Needs Starsim v4 after 4.0.0rc1 (PR #1538: cycles through `delay=`). The docstring gives the mapping; `analyses/05_emod_compiled_cascade/` compares it with the literal interpreter |
| `cascade_host.py` | `EmodHost`: EMOD′'s side of the condensed cascade (`common/cascade/`), the `CascadeHost` interface on EmodHIV, the demography and Risk; ART arms are the campaign's ARTMortalityTable configs |
| `extract_cascade.py` | The condensed cascade's parameter and flow tables for one set, read from Daniel's campaign, each row naming its source event; `python extract_cascade.py --pid 0` writes `common/cascade/tables/malawi_emod_set000_{params,flows}.csv`. `harness.make_sim(..., cascade=True, condensed=True)` runs it; `analyses/06_condensed_cascade/` compares it with the compiled campaign |
| `harness.py` | `make_sim(pid, seed, network=True, hiv=False, cascade=False, compiled=False, condensed=False)`: the modules for one parameter set, from 1960.5 in monthly steps, with `Snapshots`, an analyzer giving the same counts as EMOD's `ReportHIVByAgeAndGender` (with relationships, and with HIV, infections by stage), and with `cascade`, `EmodReport`, EMOD's ReportHIVByAgeAndGender table (by Risk, ARTStatus and InterventionStatus) |

People properties shared between the modules: `people.node` (0-2, EMOD's nodes 1-3), `people.risk` (0 LOW, 1 MEDIUM, 2 HIGH), `people.emod_age` (EMOD's age in days, after the step's ageing), `people.sti`, `people.circumcised`. `EmodNetwork.debut` holds debut ages. Edges have `p1` the man, `p2` the woman, `type` (0-3, EMOD's order), `acts` this step, `condoms` (the share of acts with a condom, from each type's sigmoid by node), `total_acts`.

## Running

In the `hivsim_paper` (or `starsim`) environment, from this folder or with it on the path:

```python
import harness
sim = harness.make_sim(pid=0, seed=1)  # about 45 s to 2025
sim.run()
sim.networks.emodnet.starts            # every relationship started: step, type, node, ages, Risk of both, drawn duration
sim.analyzers.snapshots.df             # yearly counts by node, sex, Risk and age

sim = harness.make_sim(pid=0, seed=1, hiv=True)  # about a minute to 2025
sim.run()
sim.diseases.hiv.infections            # every infection: year, relationship type (-1 seeded, -2 at birth), node, sex, age, stage of the source
sim.diseases.hiv.deaths                # every death from HIV: year, age, years since infection, CD4
```

`harness.make_sim(pid, seed, cascade=True)` is the full model (about a minute to 2025); `sim.analyzers.emodreport.df` is its report. Some of Daniel's 250 sets repeat others' parameters (e.g. 0 and 125), which EMOD runs apart only by Run_Number, so seed EMOD′ from `inputs.load(pid).run_number`.

The comparisons with EMOD are `analyses/02_emod_network/` (demography and network) and `analyses/03_emod_hiv_natural_history/` (HIV without treatment) and `analyses/04_emod_full/` (the published model); each has `run.py`, then `run.py --plot`, and needs the EMOD runs it names.

## Known differences from EMOD

- When a woman queued more than once in one type's queue is matched, EMOD removes her most recent entry, which may be in another age bin; here the entry matched is removed.
- EMOD's 50,000-candidate cap on partner selection in one age bin is not applied. Only women are chosen from a bin, and women's queues stay short (the commercial queue that grows is men's).
- Random numbers can't match EMOD's, so runs agree only in distribution.

See the two analyses for how close it is.
