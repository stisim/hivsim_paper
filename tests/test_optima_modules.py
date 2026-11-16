"""
Optima′ as agent modules (models/optima/v4/modules.py, components.py; analysis 26): the arithmetic each module took
from Optima′ declared gives the same numbers on the same people, and the runner's hybrid with every slot Optima′'s
is the modules' own sim.

    pytest tests/test_optima_modules.py    # env hp_rc3; about a minute
"""
import sys
from pathlib import Path
import numpy as np
import pytest
REPO = Path(__file__).resolve().parents[1]
OPT = REPO / 'models' / 'optima' / 'v4'
sys.path[:0] = [str(REPO), str(OPT)]


@pytest.fixture(scope='module')
def sim():
    import components as oc
    sim = oc.make_sim(seed=2, stop=2008, n_agents=20_000)
    sim.run(until=2007.5) # With people diagnosed, in care, lost and on ART in each arm
    return sim


def declared_counts(sim):
    """ The modules' people as Optima′'s state counts (states × populations, people) """
    import optima_prime as opp, modules as om
    names = [s.replace('-', '_') for s in opp.load().settings['statelabels']]
    ppl, h = sim.people, sim.diseases.hiv
    au = ppl.auids
    pop, inf = ppl.optima_pop.raw[au], h.infected.raw[au]
    place, band = h.place.raw[au], h.band.raw[au]
    circ = ppl.circumcised.raw[au]
    y = np.zeros((len(names), om.NP))
    for i, n in enumerate(names):
        if n == 'susreg': m = ~inf & ~circ
        elif n == 'progcirc': m = ~inf & circ
        else:
            c, k = n.split('_')
            m = inf & (place == om.PLACES.index(c)) & (band == opp.CD4.index(k))
        y[i] = np.bincount(pop[m], minlength=om.NP)
    return names, y*sim.current_scale


def test_force_of_infection(sim):
    """ OptimaT's force of infection (by source: place × band, from H's infectiousness and effdx) equals optima_prime.Transmission's on the same people """
    import optima_prime as opp
    names, y = declared_counts(sim)
    sim.diseases.optima_t.start_step()
    tr = opp.Transmission(opp.load())
    tr.people = lambda s: y
    assert np.allclose(sim.diseases.optima_t.rate, tr.rates(sim), rtol=1e-10, atol=0)
    assert (tr.rates(sim) > 0).any()


def test_births_and_forcing(sim):
    """ OptimaDR's births (with MTCT and PMTCT) and population forcing equal optima_prime.Births' and agents.ForcePopsize's on the same people """
    import optima_prime as opp, agents as ag, modules as om
    names, y = declared_counts(sim)
    d = sim.demographics.optima_dr
    au = sim.people.auids
    cls, dx = d.places(au)
    pop, year, dt = d.optima_pop.raw[au], d.year(), d.t.dt_year
    orig = opp.state_counts, ag.counts
    try:
        opp.state_counts = ag.counts = lambda s, n, fresh=False: np.array([y[names.index(k)] for k in n])
        b = opp.Births(opp.load()).inflows(sim)
        fp = ag.ForcePopsize(opp.load())
        arrive_s, leave_s = fp.flows(sim)
    finally:
        opp.state_counts, ag.counts = orig
    mine = d.births(pop, cls, dx, year, sim.current_scale)
    for got, key in zip(mine, ['birth', 'birth_undx', 'birth_dx']):
        assert np.allclose(got, b[key], rtol=1e-10, atol=1e-9)
    y0 = np.bincount(cls*om.NP + pop, minlength=3*om.NP).reshape(3, om.NP)*sim.current_scale
    arrive, leave = d.forcing_flows(y0, year, dt)
    sus = np.array([n in ('susreg', 'progcirc') for n in names])
    hiv = ~sus & ~fp.art
    assert np.allclose(arrive[0], arrive_s[sus].sum(0)*dt) and np.allclose(arrive[1], arrive_s[hiv].sum(0)*dt)
    for cl, m in ((0, sus), (1, hiv)): # Departures: the same rate for each state of a class
        has = y[m] > 0
        assert np.allclose(np.broadcast_to(leave[cl], y[m].shape)[has], leave_s[m][has])


def test_runner_endpoint():
    """ The runner's hybrid with every slot Optima′'s is components.make_sim(), exactly """
    import components as oc
    from common import swaps
    cfg = dict(stop=1993, pairs=dict(optima=[None]), models=dict(optima='models/optima/v4/components.py'))
    cell = swaps.cells(cfg)[0]
    a = swaps.build(dict(D='optima', T='optima', H='optima', S='optima', C='optima'), cell, cfg).run()
    b = oc.make_sim(seed=cell['seed'], stop=1993).run()
    assert swaps.differences(a, b) == []
    assert len(a.analyzers.outputs.df) > 0
