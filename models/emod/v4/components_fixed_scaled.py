"""
EMOD′-fixed (components_fixed.py) with its seeding scaled to the hybrid's population (analysis 21, option a): the same
expected seeds per person alive as EMOD′-fixed's own seeding at its scale (`seed_coverage` of Risk HIGH in 1980.85,
set 0: 1.97e-4 per person, analyses/21_seeding/outputs/rates.json), drawn among Risk HIGH (adapters.seed_per_capita).
Otherwise EMOD′-fixed's factories, so a hybrid treats its slots as EMOD′'s.
"""
from components_fixed import * # noqa: F401,F403
import components_fixed
from adapters import seed_per_capita

SEEDS_PER_CAPITA = 1.974e-4


def pars(pid):
    inp = components_fixed.pars(pid)
    inp.seeding = seed_per_capita(SEEDS_PER_CAPITA)
    return inp
