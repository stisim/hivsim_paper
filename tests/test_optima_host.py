"""
Optima's population label and Optima′'s CascadeHost (models/optima/v4/poplabel.py, host.py), as agents.

    pytest tests/test_optima_host.py    # env hp_rc3; about a minute
"""
import sys
from pathlib import Path
import numpy as np
import pytest
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'models' / 'optima' / 'v4'))
import poplabel as pl # noqa: E402


def test_label():
    """ Key populations only at 15-49 and of the right sex; everyone else by sex and band """
    age = np.array([16, 16, 49.9, 50, 3, 90, 30, 30])
    female = np.array([1, 1, 1, 1, 0, 0, 0, 0], bool)
    fsw = np.array([1, 0, 1, 1, 0, 0, 0, 0], bool)
    client = np.array([0, 0, 0, 0, 1, 0, 1, 0], bool)
    got = [pl.POPS[i] for i in pl.label(age, female, fsw, client)]
    assert got == ['FSW', 'F15-19', 'FSW', 'F50-54', 'M0-4', 'M85+', 'Clients', 'M30-34']


@pytest.fixture(scope='module')
def sim():
    import agents as ag, host as H
    from common.outputs import Outputs
    sim = ag.make_sim(n_agents=5000, end=1995, interventions=[H.OptimaTracker()], analyzers=[Outputs(H.OptimaHost)])
    sim.run()
    return sim


def test_outputs_match_disease(sim):
    """ Every sexual infection is logged with a source with HIV, and the logged number is the disease's """
    out = sim.analyzers.outputs
    r = sim.results.optima
    n = sum(np.asarray(r[k]).sum() for k in ('infection_susreg', 'infection_progcirc'))
    assert np.isclose(len(out.df)*out.scale, n)
    assert all((st != '').all() for _, st, _ in out.host.sources) # Each in a state with HIV when drawn


def test_moves(sim):
    """ Testing diagnoses; starting ART marks ever treated; stopping puts people in lost, as interrupted """
    h = sim.analyzers.outputs.host
    au = sim.people.auids
    undx = au[h.state_in(au, ['undx'])][:20]
    h.test(undx)
    assert h.cascade_stage(undx)['diagnosed'].all()
    h.start_art(undx[:10], 'suppressive')
    st = h.cascade_stage(undx[:10])
    assert st['on_art'].all() and st['suppressive'].all() and st['ever_treated'].all()
    h.stop_art(undx[:10])
    assert h.state_in(undx[:10], ['lost']).all()
    assert (h.cascade_stage(undx[:10])['category'] == 'interrupted').all()
    assert (h.cascade_stage(undx[10:])['category'] == 'diagnosed_never_treated').all()
