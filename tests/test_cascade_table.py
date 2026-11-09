"""
Optima's cascade as a table (common/cascade/tables/malawi_optima_params.csv; pathways.py, rates.py): Optima′ running it.

    pytest tests/test_cascade_table.py    # env hp_rc3; about a minute
"""
import sys
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'models' / 'optima' / 'v4'))
TABLE = REPO / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv'


def test_ode_table_is_own():
    """ As an ODE, the table gives Optima′'s own cascade (catches a row or rule read differently from Optima's inputs) """
    import optima_prime as opp, compare as cmp
    from common.cascade.params import Params
    runs = [opp.make_sim(components=cmp.COMPONENTS['E'], end=2010, **kw).run().results.optima for kw in ({}, dict(cascade=Params(TABLE)))]
    for k in ('undx_gt200', 'dx_gt50', 'care_gt350', 'lost_gt50', 'svl_gt500', 'usvl_gt200', 'infection_susreg'):
        assert np.allclose(np.asarray(runs[0][k]), np.asarray(runs[1][k]), rtol=1e-9, atol=1e-6), k


def test_agents_table_hits_target():
    """ As agents through OptimaHost, the table's number on ART is met (catches moves that don't reach the host's states) """
    import agents as ag, host as H, optima_prime as opp
    from common.cascade.params import Params
    from common.cascade.pathways import TableCascade
    P = Params(TABLE)
    sim = ag.make_sim(n_agents=20_000, end=2008, cascade=P, interventions=[TableCascade(P, H.OptimaHost), H.OptimaTracker()]).run()
    tc = sim.interventions.table_cascade
    au = sim.people.auids
    place = tc.places(au)
    on_art = np.isin(place, [5, 6]).sum()*sim.current_scale
    target = P.value('target', 'on_art', sim.t.now('year'))
    assert abs(on_art/target - 1) < 0.05, (on_art, target)
    assert (place == 3).sum() > 0 and (place == 4).sum() > 0 # Some in care and lost: linkage and loss move people
