""" Optima's cascade as models/optima/v4/components.py (no arm switching), with its PMTCT number drawn per person (see optima_c_switch.py) """
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('optima_c_switch', Path(__file__).resolve().parent / 'optima_c_switch.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.SWITCH, base.PMTCT = False, True
pars, make_C = base.pars, base.make_C
