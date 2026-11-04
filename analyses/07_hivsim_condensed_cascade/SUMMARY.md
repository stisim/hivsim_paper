# 07: HIVsim with EMOD's condensed cascade — summary

**Answer:** the condensed cascade drives HIVsim (STIsim's `sti.HIV`) through the same `CascadeHost` interface as EMOD′, after five small additions to STIsim (WHO stage and first symptoms from CD4, PMTCT for mothers not on ART, ART in an arm the caller chooses, PrEP courses ended by HIV; branch `rs-cascade-host` on STIsim PR #607). The host behaves. Swapped in for HIVsim's own care, EMOD's cascade halves ART coverage (2020: 44% against 81%), so HIV deaths (54k against 2.3k a year; UNAIDS 16k) and new infections (62k against 27k; UNAIDS 21k) are much higher, with people with HIV about the same (about 1.0M), because more die and more are infected.

Versions: Starsim v4 ebd52c1a, STIsim `rs-cascade-host` 393e3c6, HIVsim Malawi from hivsim_malawi 803b2f7 (best run of experiment 03, provisional; copied to `models/hivsim/malawi/`), 10,000 agents, 50 seeds per variant, about 8 s per run.

## The host (condensed variant)

- Events match HIVsim's: 12-week pregnancies equal conceptions (744.5 against 744 a year in 2010), 6-week-olds within 1% of births, 42–78 first symptoms a year; ART starts and stops match the cascade's deliveries.
- No one is on ART outside care on ART (0 exceptions in any run); 80% of those on ART are suppressed, as the table's arms give.
- PMTCT reaches HIVsim: child infections over 2001–2021 fall from 328 to 258 agents per run (−21%) with it.
- One EMOD-faithful leak: EMOD's path for symptomatic people who didn't present goes straight to staging without a test (`symptomatic.p_present_again`), leaving 45 person-years on ART without a diagnosis over all runs.

## The swap, 2020 medians (`cascade_over_time.png`)

| | HIVsim own care | HIVsim, condensed cascade | EMOD′, condensed (06) | Data |
|---|---|---|---|---|
| ART coverage | 80.7% | 43.7% | 40.3% | |
| Off ART after starting it, share of PLHIV | 0% | 34.2% | 42.3% | |
| Diagnosed | 81% | 80% | 84.5% | |
| Suppressed among on ART | 96% | 80% | 80% | |
| HIV deaths a year | 2.3k | 54k | | UNAIDS 16k |
| New infections a year | 27k | 62k | | UNAIDS 21k |
| Infections of children a year | 1.5k | 7.6k | | |
| People with HIV | ≈1.0M | ≈1.0M | | |

On the common targets, condensed over own (medians): ART coverage 0.63, on ART 0.59, VLS among PLHIV 0.46, on ART among diagnosed 0.54, AIDS deaths 1.06 and new infections 1.03 over all years, PLHIV 1.00. Median target score: own 5.9, condensed 49.5.

## Why they differ

- **ART by numbers vs by mechanism.** HIVsim's ART follows the programme's numbers on ART by construction (everyone diagnosed ends up on ART, nobody is off it after starting), with suppression held at MPHIA's 96%. EMOD's cascade starts people through staging and linkage and stops them at its dropout rate, and with the dropout trap (from 2009, 85–95% of those who drop out are never offered ART again) 34% of PLHIV are off ART after starting it by 2020 (EMOD′: 42%). Its arms give 80% suppression on ART.
- **The cascade carries over from EMOD′.** In HIVsim it gives nearly the same shares as in EMOD′ (coverage 44% against 40%, diagnosed 80% against 85%, suppressed on ART 80% against 80%): the cascade's behaviour is mostly its own, not the disease's.
- **Neither combination is recalibrated.** EMOD's cascade was fitted with EMOD's disease, and HIVsim's parameters with its own care, so the swapped model fits the targets worse. Note also that HIVsim's own care gives far fewer HIV deaths than UNAIDS (2.3k against 16k in 2020); its calibration is provisional.

## Gaps

WHO stage and first symptoms come from CD4 in HIVsim; staging reads STIsim's CD4 nadir rather than the lowest count measured; PrEP is not renewed for those already on it; circumcision before the programmes is held at HIVsim's pre-2010 survey level (EMOD keeps traditional circumcision outside the cascade); STIsim's tests fail on Starsim v4 head independently of these changes (head requires `ss.Pregnancy` for `MaternalNet`), though HIVsim Malawi runs on it.

## Next steps

- Recalibrate (HIVsim with the condensed cascade, or the cascade's linkage and dropout with HIVsim).
- The same swap for Optima′ (analysis 08).
- Decide whether the dropout trap belongs in the paper's shared cascade or only in EMOD's (swap design D10).
