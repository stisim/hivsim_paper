"""
EMOD′ with Daniel's fixes (hivsim_paper PR #4; swap_design D10): background mortality read as deaths per person
(`mortality`) and InterventionStatus released at ART dropout so people can start ART again (`reenrol`); see
inputs.FIXES. Otherwise EMOD′'s own factories (components.py), so a hybrid treats its slots as EMOD′'s. In a config,
`models: {'emod_fixed': 'models/emod/v4/components_fixed.py'}`. The other fixes in inputs.FIXES are off here.
"""
from components import * # noqa: F401,F403 (the slot factories, HOST, SCALE, BURNIN)
import components

FIXES = ('mortality', 'reenrol')


def pars(pid):
    return components.pars(pid, FIXES)
