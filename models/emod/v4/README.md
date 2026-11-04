# EMOD′ (EMOD-HIV re-expressed in Starsim v4)

So far: the sexual network, with Risk and sex work, and a harness to run it on its own. No HIV yet.

| File | What it is |
|---|---|
| `inputs.py` | `load(pid)`: the parameters for one of Daniel's parameter sets, read from the EMOD inputs that `reference/run_emod.py` builds for it, so his quirks come with them (e.g. the Risk overlay's node order) |
| `network.py` | `EmodNetwork`, an `ss.SexualNet`: EMOD's pair-formation algorithm and its four relationship types (transitory, informal, marital, commercial): queues with per-person entry rates set during the 15-year burn-in and then frozen, matching on the joint age tables and Risk weights, the extra-relational flag chain and randomly rounded maxima, Weibull durations in whole steps, and Poisson coital acts with dilution and a forced first act. The module docstring lists the steps and where each is in the EMOD source |
| `risk.py` | `EmodRisk`, an intervention: `Risk` drawn by node at first sight, and the campaign's sex work (at debut from 1961, will ever be FSW or client by node; uniform delay to Risk HIGH, then Weibull or uniform time to MEDIUM) |
| `harness.py` | `make_sim(pid, seed)`: the network and risk modules on a minimal EMOD-like population (`EmodDemography`: EMOD's initial ages, node sizes at 1% scale, pregnancy and mortality tables, the latter 1000 times too low as in EMOD), from 1960.5 in monthly steps, with `Snapshots`, an analyzer giving the same counts as EMOD's `ReportHIVByAgeAndGender` with relationships. The demography is only what the network needs; it is not the EMOD′ demography component |

People properties shared between the modules: `people.node` (0-2, EMOD's nodes 1-3), `people.risk` (0 LOW, 1 MEDIUM, 2 HIGH). `EmodNetwork.debut` holds debut ages. Edges have `p1` the man, `p2` the woman, `type` (0-3, EMOD's order), `acts` this step, `condoms` (the share of acts with a condom, from each type's sigmoid by node), `total_acts`.

## Running

In the `hivsim_paper` (or `starsim`) environment, from this folder or with it on the path:

```python
import harness
sim = harness.make_sim(pid=0, seed=1)  # about 45 s to 2025
sim.run()
sim.networks.emodnet.starts            # every relationship started: step, type, node, ages, Risk of both, drawn duration
sim.analyzers.snapshots.df             # yearly counts by node, sex, Risk and age
```

The comparison with EMOD is `analyses/02_emod_network/` (`run.py`, then `run.py --plot`; it needs the EMOD runs it names).

## Known differences from EMOD

- When a woman queued more than once in one type's queue is matched, EMOD removes her most recent entry, which may be in another age bin; here the entry matched is removed.
- EMOD's 50,000-candidate cap on partner selection in one age bin is not applied. Only women are chosen from a bin, and women's queues stay short (the commercial queue that grows is men's).
- Sexual debut is seen at the step after it (here the risk module runs after the network; in EMOD the debut event fires as people age within the step). Sex-work delays are months to years, so this shifts them by a month.
- Random numbers can't match EMOD's, so runs agree only in distribution.

See `analyses/02_emod_network/SUMMARY.md` for how close it is.
