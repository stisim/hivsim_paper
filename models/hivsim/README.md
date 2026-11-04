# HIVsim (STIsim's HIV, Malawi)

| File | What it is |
|---|---|
| `malawi/` | HIVsim Malawi (`malawi.py`, `targetfit.py`, `data/`) copied unchanged from hivsim_malawi at commit 803b2f7 (3 Nov 2026), while its calibration is still under way. `targetfit.py` measures the common targets in a HIVsim run (pass it `common/targets.csv`, which is the same file) |
| `cascade_host.py` | `StisimHost`: HIVsim's side of the condensed cascade (`common/cascade/`), the `CascadeHost` interface on `sti.HIV`, `sti.ART`, `ss.Pregnancy` and `sti.StructuredSexual`; and `art_for_cascade()`, an `sti.ART` that only holds ART's states |
| `swap.py` | `make_sim(variant, seed)`: HIVsim Malawi with its own care (`own`), or with the condensed cascade and EMOD's Malawi table (`condensed`); `BEST` is hivsim_malawi experiment 03's best run |
| `reference/` | Where the finalized HIVsim Malawi will go |

Needs STIsim with the host's methods (`HIV.who_stage()`, `HIV.new_symptoms()`, `HIV.start_pmtct()`, `ART.start_art(effective=)`, PrEP courses ended by HIV): branch `rs-cascade-host` (from `rs-art-split`, PR #607), and Starsim v4 after 4.0.0rc1 for the condensed cascade:

```bash
PYTHONPATH=<stisim rs-cascade-host>:<starsim v4> python -c "import swap; swap.make_sim('condensed').run()"
```

How the interface maps onto HIVsim, and what is approximate:
- WHO stage and first symptoms come from the CD4 count (stage 4 and symptoms below 200), not from a clinical course; first symptoms come once per person.
- Staging reads STIsim's CD4 nadir (lowest off ART), not the lowest count measured.
- The ART arms are STIsim's suppressed and non-suppressed ART (its mortality and infectiousness on each), chosen by the cascade.
- PMTCT for mothers not on ART reduces their children's susceptibility in pregnancy and breastfeeding, like ART's `pmtct_efficacy`.
- Flows (VMMC, PrEP) are summed over EMOD's nodes; PrEP's restriction to risk MEDIUM/HIGH is STIsim's risk groups 1-2.
- Circumcision before the programmes is held at HIVsim's pre-2010 survey level (`swap.traditional_vmmc()`), as EMOD keeps traditional circumcision outside the cascade.
