# 02: EMOD′ sexual network against EMOD

**Question:** how closely does the first Starsim v4 translation of EMOD-HIV's sexual network (`models/emod/v4/`: the pair-formation algorithm, the four relationship types, Risk and sex work) reproduce the network in Daniel's EMOD Malawi model, and where does it differ?

**What's run:**
- **EMOD:** the same 10 parameter sets as analysis 01 (ids 0, 25, ..., 225), at 1% population, to 2025, with relationship reports switched on: `models/emod/reference/run_emod.py --sets 0 25 50 75 100 125 150 175 200 225 --cores 10 --out net10 --network`. This is the published model, with HIV.
- **EMOD without HIV:** sets 0, 100 and 200 with the epidemic seeding removed (`--sets 0 100 200 --out net3_nohiv --network --no-hiv`), since the v4 harness has no HIV and HIV deaths change the population from the late 1980s.
- **EMOD′ network:** `models/emod/v4/harness.py`, the network and risk modules on EMOD's population (initial ages, node sizes, pregnancy and mortality tables, 3 nodes, 1% scale, from 1960.5, monthly steps), no HIV, for the same 10 sets, one run each.
- `run.py` reduces both to the same summaries (`outputs/`), and `run.py --plot` draws the figures (`figures/`).

**Compared, over time:** people in relationships by type, sex and age; concurrency; relationship durations by type; partner age differences; commercial contacts per FSW and per client; formation rates by type; also queue sizes, population and the number of FSW and clients, which explain the others.

**Not compared:** transmission, or anything that depends on HIV. Runs are compared as distributions over parameter sets and seeds; EMOD's random number streams can't be matched.
