# EMOD rerun: status and background-mortality problem (3 Nov 2026)

## What we did

- Reran Daniel Citron's EMOD-HIV Malawi model (`models/emod/original/`) with the published `emod-hiv` 2.37.0 binary instead of his NYU cluster build (setup in `models/emod/README.md`).
- `models/emod/run_emod.py` fills in his resampled parameter sets the same way dtk-tools did. Runs stop at 2025 to match the other models (his setup runs to 2100).
- Test: 10 of his 250 sets (ids 0, 25, ..., 225, spanning his fit scores), 1% population, about 11 minutes per run.
- `models/emod/compare_emod.py` compares national outputs with his calibration targets: `results/emod/sub10_2025_vs_targets.png`.

## Results

- **Match:** prevalence aged 15-49 (women and men) is inside his target error bars, and the number on ART is close to his targets.
- **Don't match:** population is about twice the census (about 25M against 13.1M in 2008), PLHIV is about double, and ART coverage comes out around 40% instead of 80-90%.

## Cause: background mortality about 1000 times too low

- `MortalityDistributionMale/Female` in `Demographics_Malawi.json` hold annual deaths **per person** (infant mortality 0.18 in 1960).
- Their `ResultScaleFactor` is 1/(365*1000), the factor for deaths **per 1000**, so every rate is divided by an extra 1000.
- In the output, non-HIV deaths are about 0.03 per 1000 per year (Malawi: about 10-15), so the population grows about 4% a year.
- His runs use the same input files, so they almost certainly have the same problem. It didn't show up because his calibration rescaled population and on-ART counts to 8.26M aged 15-49 in 2018 and gave population only 2.5% of the weight. That hides the overall size but not the growth or the age structure.

## Open decision

- **Fix** the scale factor (one line): his 250 parameter sets were fitted without background mortality, so they probably no longer fit his targets and would need recalibrating.
- **Reproduce as published:** EMOD enters the comparison with a population that's far too large.
- Either way, raise it with Daniel.

## Figures against the shared targets (`results/emod/sub10/figures/`)

`models/emod/emod_targets.py` measures EMOD runs against `results/targets.csv` and draws the same three figures as the HIVsim experiments, using `models/hivsim/targetfit.py`. Our runs now report from 1990 in 5-year age bins; this is a reporting change only. Same 10 sets, run to 2025.

- **Prevalence by age and sex** matches DHS and MPHIA well. Prevalence was most of his fit (80% of the weight).
- **UNAIDS numbers are off from about 1995:** PLHIV reaches about 2M by 2020 (data about 1M), new infections stay around 150k a year (data falls to about 20k), and AIDS deaths peak too high and rise again after 2015.
- **ART coverage levels off at about 40% from 2014** (data reaches 86% by 2021), even though about 90% of PLHIV are diagnosed, matching MPHIA. So in this model most diagnosed people never start ART or don't stay on it. The number on ART matches the programme data, so the low coverage comes from PLHIV being about double. ART in his campaign is triggered by events, not handed out as fixed numbers, so the link to the mortality bug is not mechanical. I haven't found the cause yet.
- **Not shown:** viral suppression (EMOD has no viral load) and FSW prevalence (the model has risk groups but no sex-work flag).
- **Check of the 2.37 build:** across the 10 sets, our scores agree only weakly with his likelihoods (Spearman 0.43 in the right direction). That's too few sets to say much.
