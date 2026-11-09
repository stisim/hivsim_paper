Responding to known quirks of original setup, in README.md

* In the calibration form, `Pr Ex Infmrl Fem LOW` maps to the male parameter (`...INFORMAL...LOW.Prob_Extra_Relationship_Male`), so it overwrites `Pr Ex Infmrl Male LOW`, and the female value is never set.
  * Acknowledged, it is unlikely this would affect the calibration in a significant way
* At 1% population, the 2011 VMMC target in the Northern region (54 men) rounds to zero.
  * Acknowledged, it is unlikely this would have a
* Background mortality is about 1000 times too low. `MortalityDistributionMale/Female` in `Demographics_Malawi.json` hold annual deaths per person (infant mortality 0.18 in 1960), but `ResultScaleFactor` is 1/(365*1000), the factor for deaths per 1000. Non-HIV deaths come out at about 0.03 per 1000 per year, and the population reaches about twice the census by 2008. His calibration rescaled population and on-ART counts to 8.26M aged 15-49 in 2018 and gave population 2.5% of the weight, which hides the size but not the growth or the age structure.
  * This changes the population (and PLHIV), as the current ersion has background mortality that is too low.
  * This does not affect the calibration by much - all simulation results are rescaled by a population scaling factor to match census population data in 2018 - so prevalence and ART coverage (rates) remain the same.
  * This does affect the future trajectory of the epidemic, where the PLHIV population and incidence are lower when we use the correct background mortality.
  * I'm guessing that the original version I calibrated had the correct scaling factor and at some point in the shuffle that got overwritten. I checked our configuration files for our other models and they do not have this same problem.
* People who drop out of ART can't restart (from 2009): `InterventionStatus` stays `OnEffectiveART`/`OnNonSuppressiveART` after dropout, and the ART-start events are blocked for those values. See `analyses/01_reference_fits/SUMMARY.md` and the other problems listed in `docs/emod_inventory.md`.
  * This error was likely introduced when we converted the model to use the ARTMortalityTable.
  * It can be fixed by switching InterventionStatus back to OffART once someone drops out.
