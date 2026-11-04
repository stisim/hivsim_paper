# 08: The condensed cascade in Optima′ — summary

**Answer:** `common/cascade/rates.py` turns the condensed cascade's parameter table into rates between Optima′'s care states (`make_sim(cascade=...)`; Optima′ without a table is unchanged, and `acceptance.py` gives the same results as before).
- **Consistent with the cascade it stands for.** In EMOD′ running the same table (set 0, 3 seeds), the rates out of each care state seen among people with HIV match what `rates.py` gives for the same person-time: within about 5% for the undiagnosed (2008–2020), on ART and lost, and within 10–26% for the diagnosed after 2012. The hidden shares match closely (testing-loop share of the undiagnosed 0.40 seen against 0.38; out-rather-than-lost share of dropouts 0.91 against 0.89, 2019).
- **In Optima′, EMOD's cascade behaves as it does in EMOD′ and HIVsim.** ART coverage in 2020 is 40% (EMOD′ 41%, HIVsim with the same cascade 44%, Optima′'s own cascade 78%, data 84%), 38% of PLHIV are lost, and new infections are 91k against 34k with Optima′'s own cascade.

## The mapping (details in `common/cascade/README.md`)

Care states are Optima's: undiagnosed, diagnosed not in care, in care (pre-ART), lost (dropped out of ART), on ART suppressed and not suppressed.
- **Exact given those states:** the routine testing loop (exponential intervals), the outcomes of a positive test (immediate ART by sex, lost at staging, linked if eligible, pre-ART if not), ANC staging, the ART arm by sex and age, dropout and its outcomes including the dropout trap, the PMTCT options and Option-B block (as an expected efficacy by year), testing at first symptoms and the infant test on infected births.
- **Approximate:** two hidden shares that aren't Optima states (in the testing loop vs waiting in uptake; out vs lost for good) are followed by a small ODE per pool; the 182-day pre-ART check becomes a rate, and staging and the ART start delay are taken as instant; eligibility is the share of each CD4 band eligible (time spread evenly over a band, WHO stage from EMOD's natural history, pregnancy from the birth rate); first symptoms as entry into CD4 < 200; VMMC spread over men of eligible ages with the regions summed; PrEP as coverage among women of eligible ages (EMOD's Risk targeting lost); PMTCT only for undiagnosed or diagnosed mothers not in care.
- **Not representable without extra states:** per-person memory and schedules (lowest CD4 ever, next test date, each mother's PMTCT window, symptoms recurring after stopping ART) and extra event-triggered testing pathways whose trigger isn't a compartment transition.
- **HIVsim/STIsim-style cascades** map directly (described, not implemented): an annual testing probability p becomes −ln(1−p); `HIVCare`'s loss and return rates are already rates; `ART(coverage=...)` is a target on levels, like Optima's `numtx`.

## (i) Rates against EMOD′ (`emod_rates.png`, `emod_hidden_shares.png`)

Starsim ebd52c1a. Total rate out of each state, seen in EMOD′ against `rates.py`, per person-year:

| Year | Undiagnosed | Diagnosed | On ART, suppressive | On ART, non-suppressive | Lost |
|---|---|---|---|---|---|
| 2010 | 0.280 / 0.288 | 0.117 / 0.144 | 0.173 / 0.178 | 0.168 / 0.178 | 0.022 / 0.021 |
| 2015 | 0.360 / 0.378 | 0.315 / 0.282 | 0.175 / 0.178 | 0.175 / 0.178 | 0.021 / 0.023 |
| 2020 | 0.450 / 0.429 | 0.370 / 0.273 | 0.173 / 0.178 | 0.173 / 0.178 | 0.0075 / 0.0076 |

Where pairs of states differ, by construction: EMOD passes people through linking (1–120 days, counted as in care) before ART, which `rates.py` skips, so undiagnosed → ART is 0.32 from `rates.py` against 0.24 seen and undiagnosed → care about 0 against 0.14; and `rates.py` applies pre-ART loss to everyone in care, while after 2016 EMOD's care is all linking. Exits from the diagnosed are about 26% low by 2020, not verified: likely EMOD firing first symptoms again after ART dropout (one event per infection in `rates.py`) and staging by the lowest CD4 ever rather than the current band.

## (ii) Optima′: its own cascade against EMOD's as rates (`optima_cascade.png`)

Starsim 4.0.0rc1, dt 0.02, to 2025. Score on the common targets: own cascade 6.35 (the Optima reference in analysis 01 scores 6.50); EMOD's cascade 59.8 (60.3 without VMMC and PrEP).

| 2020 | Own cascade | EMOD's cascade |
|---|---|---|
| On ART | 855k | 486k |
| ART coverage | 78% | 40% |
| Diagnosed (men / women) | 89% / 91% | 83% / 88% |
| On ART among diagnosed | 87% | 48% |
| Suppressed among on ART | 93% | 79% |
| Lost, share of PLHIV | 2% | 38% |
| New infections | 34k | 91k |
| HIV deaths | 9.8k | 37k |
| People with HIV | 1.11M | 1.19M |

Relative to each model's own cascade, EMOD's cascade changes the shares by nearly the same factors in Optima′ as in EMOD′ (analysis 06): ART coverage 0.618 against 0.610, diagnosed 0.946 against 0.949, suppressed among on ART 0.85 against 0.87, suppressed among PLHIV 0.45 against 0.46.

**What differs, and why:**
- **Number vs rate.** Optima starts ART by a forced number (`numtx`, lowest CD4 first); EMOD's cascade starts it by rates, so the number on ART becomes an output and falls far short of the data.
- **The dropout trap.** Optima's lost compartment returns people to care; EMOD's trap holds them out: 38–44% of PLHIV lost by 2020–2025, against 2%.
- **Suppression.** Optima's cascade has 63% suppressed at ART start plus failure and viral-load monitoring; EMOD's arms are fixed at ART start.
- **Pre-ART waiting.** Optima's own cascade holds 30–48% of PLHIV in care waiting for ART in 2000–2010, because `numtx` caps starts; with EMOD's rates almost nobody waits, but many diagnosed people wait outside care while ineligible (up to 38% of PLHIV around 2008).
- **VMMC and PrEP.** 0.78M men circumcised by 2020 cut 2020 infections by 2.7% (5% by 2025); PrEP is negligible.

## Next steps

- Decide whether the swaps take the cascade with or without EMOD's quirks (the dropout trap, the Option-B block): swap design D10.
- Test whether repeated symptoms and lowest-CD4 staging explain the gap in exits from the diagnosed.
- Run with the other parameter sets' tables.
- Optionally split Optima′'s lost and undiagnosed compartments to make the hidden shares exact.
