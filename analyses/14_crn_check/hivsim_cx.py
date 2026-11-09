"""
HIVsim's slot factories with EMOD's condensed cascade as C (swap.make_sim('condensed'): EMOD's Malawi table set 0, on
StisimHost, with sti.ART holding ART's states and circumcision held at its pre-2010 level). Every other slot is
HIVsim's own (models/hivsim/components.py), so as a C from this model, it is step 1 of the path (C) for HIVsim.
"""
import sys
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'models' / 'hivsim')) # HIVsim's factories, whichever model the runner loads first
from components import HOST, SIM, pars, make_D, make_R, make_N, make_RN, make_T, make_H, make_S, make_V, make_M # noqa: F401
import components


def make_C(pars, ctx):
    return components.make_C(pars, SimpleNamespace(**vars(ctx), cascade='condensed'))
