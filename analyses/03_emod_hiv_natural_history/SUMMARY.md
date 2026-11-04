# 03: EMOD′ HIV transmission and natural history against EMOD — summary

**Answer:** without testing or treatment, EMOD′ gives EMOD's epidemic within the spread of its own seeds, for all 3 parameter sets: prevalence by age, sex and node, people with HIV, sexual infections, deaths from HIV, time from infection to death, and who infects whom (by relationship type, stage and years since infection). The one notable difference, set 0's lower totals, comes from the seeding being stochastic in node 1.

## What matches

EMOD (one run per set) against EMOD′ (median [range] of 5 seeds), with the campaign reduced to sex work, STI co-infection, traditional circumcision and seeding (`run_emod.py --natural-history`).

- **Prevalence 15–49** (`totals.png`, `prevalence_by_age.png`, `prevalence_by_node.png`):

  | | Set 0 | Set 100 | Set 200 |
  |---|---|---|---|
  | Women 1995 | 0.141 / 0.145 | 0.151 / 0.143 | 0.154 / 0.149 |
  | Women 2020 | 0.135 / 0.120 [0.115–0.136] | 0.128 / 0.125 [0.115–0.129] | 0.156 / 0.158 [0.137–0.166] |
  | Men 2020 | 0.089 / 0.083 | 0.091 / 0.087 | 0.120 / 0.119 |

  The shape by age (peak at 30–45, rising with time) matches.
- **People with HIV, 2020**: 24.5k / 22.0k, 24.0k / 23.0k, 29.0k / 29.4k (agents, 1% scale).
- **Sexual infections, 2010–24**: 31.1k / 28.3k [27.7–32.2k], 30.8k / 30.0k, 36.0k / 37.3k.
- **HIV deaths, 2010–24**: 30.2k / 27.6k [25.6–30.5k], 28.7k / 28.4k, 34.9k / 35.5k.
- **Time from infection to HIV death** (`natural_history.png`): the distributions overlie each other; mean 8.51 whole years in both (deaths 2000–24).
- **By relationship type** (`transmissions.png`), share of sexual infections 2000–24: transitory 15.4% / 15.3%, informal 39.9% / 39.4%, marital 32.0% / 32.6%, commercial 12.7% / 12.7%. Each set's own path over time is followed (transitory falling after 2000, marital rising).
- **Stage of the partner who transmitted**, 2000–24: acute 23% / 21%, latent 59% / 61%, AIDS 18% / 18%. The distribution of that partner's years since infection overlies EMOD's.

## What doesn't, and likely causes

1. **Set 0 runs low** (women's prevalence 0.120 against 0.135 in 2020) because node 1's epidemic never starts in 4 of 5 seeds. EMOD′ seeds 0.6 people there on average (1.4 and 2.4 in sets 100 and 200), and there's no migration between nodes. EMOD's epidemic started in node 1 in all 3 runs, which at these rates would happen about 13% of the time: probably chance, not shown. Next: count EMOD's seeds over more runs.
2. **The acute share of transmissions** is 2.5 points higher, and **the AIDS share of people with HIV** is 6.1% against 6.7% (2005–24). Probably reporting timing rather than mechanism: EMOD's TransmissionReport takes the source's stage after its own progression in that step, and its population report likely includes people dying that step (who are in AIDS). Not verified.
3. **Transmission at birth** isn't reported by EMOD. Prevalence at ages 0–4 matches for set 200 (0.024 / 0.024) and is about 13% lower in sets 0 and 100, following their lower adult prevalence.

## Caveats

- This is not the published model: without ART the epidemic is larger and keeps growing, so this checks the mechanisms, not the fit.
- 3 sets, one EMOD run each.

## Next steps

- Testing, care and ART (including the dropout trap), then VMMC, PrEP and PMTCT, and compare with the full published EMOD.
