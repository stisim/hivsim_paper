""" As optima_c_switch.py, with EMOD′'s own VMMC and PrEP flows (see there) """
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('optima_c_switch', Path(__file__).resolve().parent / 'optima_c_switch.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.PREVENTION = True
pars, make_C = base.pars, base.make_C
