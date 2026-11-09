""" EMOD′ with fix art_male (models/emod/v4/inputs.py FIXES) """
import _variant # noqa: F401 (puts EMOD′ on the path)
from components import * # noqa: F401,F403
import components
import inputs

FIXES = ('art_male',)


def pars(pid):
    return components.pars(pid, FIXES)
