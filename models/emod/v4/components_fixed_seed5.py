"""
EMOD′-fixed (components_fixed.py) seeding 5× its share of Risk HIGH (analysis 21, option c: set 0's 1.09% → 5.45%).
Otherwise EMOD′-fixed's factories, so a hybrid treats its slots as EMOD′'s.
"""
from components_fixed import * # noqa: F401,F403
import components_fixed

SEED_MULT = 5


def pars(pid):
    inp = components_fixed.pars(pid)
    inp.campaign.seed_coverage *= SEED_MULT
    return inp
