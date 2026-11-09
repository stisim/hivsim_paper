# EMOD

- `reference/source/`: Daniel Citron's EMOD-HIV Malawi model as used for MIHPSA (inputs, 250 resampled parameter sets, calibration form, dtk-tools scripts for the NYU cluster), unchanged.
- `reference/run_emod.py` reruns his parameter sets here, without dtk-tools or the cluster. It maps each parameter set to model inputs as `source/optim_script.py` does and runs the EMOD binary directly. Runs stop at 2025 by default (his setup runs to 2100; `--years` changes it). The report is widened to start in 1990, use 5-year age bins and split by `Risk` and `ARTStatus`; this changes the output only. Outputs go to `raw_results/emod/<out>/setNNN/output/`.
- `reference/emod_targets.py --out <out>` measures the runs against `common/targets.csv` and draws the shared figures in `analyses/01_reference_fits/`.
- `v4/`: EMOD′, the Starsim v4 translation; so far the sexual network, Risk and sex work, with a harness (see `v4/README.md`). `reference/run_emod.py --network` adds the relationship reports it is compared against, and `--no-hiv` drops the seeding, and `--natural-history` keeps only the campaign events EMOD′ has so far (sex work, STI co-infection, traditional circumcision, seeding) and reports transmissions and HIV deaths. `--transmission` turns on Daniel's transmission reports (`source/InputFiles/Templates/config_transmission.json`: each transmission, diagnosis and ART start and stop events, from 1960) for `common/transmission.py`, and `--seeds` runs each set with several seeds.

## Binary

Daniel ran `EradicationMalaria-Ongoing_2023_01_05`. We use the published build in `emod-hiv` 2.37.0 (EMOD main, Aug 2026). `to_v237()` in `reference/run_emod.py` lists the input changes this needs.

One-time setup:

```
conda create -n emod -c conda-forge python=3.13 mpich numpy pandas openpyxl pip
conda run -n emod python -m pip install emod-hiv==2.37.0
conda run -n emod python -c "import emod_hiv.bootstrap as b; b.setup('raw_results/emod_bin')"   # from the repo root
mkdir raw_results/emod_bin/lib && ln -s ~/miniconda/envs/emod/lib/libmpi.so.12 raw_results/emod_bin/lib/libmpich.so.12
chmod +x raw_results/emod_bin/Eradication
```

The binary links to `libpython3.13` and `libmpich.so.12`. On Azure VMs MPI fails to start unless `UCX_TLS=self,sm,tcp UCX_NET_DEVICES=lo` (set in `reference/run_emod.py`).

## Known quirks of the original setup (reproduced as is)

- In the calibration form, `Pr Ex Infmrl Fem LOW` maps to the male parameter (`...INFORMAL...LOW.Prob_Extra_Relationship_Male`), so it overwrites `Pr Ex Infmrl Male LOW`, and the female value is never set.
- At 1% population, the 2011 VMMC target in the Northern region (54 men) rounds to zero.
- Background mortality is about 1000 times too low. `MortalityDistributionMale/Female` in `Demographics_Malawi.json` hold annual deaths per person (infant mortality 0.18 in 1960), but `ResultScaleFactor` is 1/(365*1000), the factor for deaths per 1000. Non-HIV deaths come out at about 0.03 per 1000 per year, and the population reaches about twice the census by 2008. His calibration rescaled population and on-ART counts to 8.26M aged 15-49 in 2018 and gave population 2.5% of the weight, which hides the size but not the growth or the age structure.
- People who drop out of ART can't restart (from 2009): `InterventionStatus` stays `OnEffectiveART`/`OnNonSuppressiveART` after dropout, and the ART-start events are blocked for those values. See `analyses/01_reference_fits/SUMMARY.md` and the other problems listed in `docs/emod_inventory.md`.
