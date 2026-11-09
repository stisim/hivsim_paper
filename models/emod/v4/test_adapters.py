"""
Tests of the EMOD′ adapters on HIVsim's slots (adapters.py). From this folder, with STIsim v4 on the path:

    pytest test_adapters.py
"""
import sys
from pathlib import Path
import numpy as np
import numba as nb
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adapters # noqa: E402

N = 2000


@nb.njit
def other_numba_code(n):
    """ Numba code elsewhere that draws from the global stream """
    s = 0.0
    for _ in range(n):
        s += np.random.random()
    return s


def summary(sim):
    out = [sim.results.alive.values.copy()]
    if 'hiv' in sim.diseases:
        out.append(sim.diseases.hiv.results.new_infections.values.copy())
    if 'emodnet' in sim.networks:
        out += [sim.networks.emodnet.results[f'new_{k}'].values.copy() for k in 'TIMC']
    return np.concatenate(out)


@pytest.mark.parametrize('network', ['sti', 'emod'])
def test_runs(network):
    """ Each hybrid builds and runs 1985 to 2025, with people in every node, Risk set for adults, and HIV seeded and its records kept """
    sim = adapters.make_sim(pid=0, seed=2, n_agents=N, stop=2025, network=network)
    sim.run()
    ppl, hiv = sim.people, sim.diseases.hiv
    au = ppl.auids
    adults = au[ppl.age.raw[au] >= 15]
    assert len(np.unique(ppl.node.raw[au])) == 3
    assert (ppl.risk.raw[adults] >= 0).all()
    assert ppl.sti.raw[adults].any()
    assert len(hiv.results.new_infections) == len(sim.t)
    assert np.isclose(ppl.emod_age.raw[adults]/365, ppl.age.raw[adults], atol=0.1).all() # emod_age is the age at the end of the step (the people were aged since)
    if network == 'emod':
        net = sim.networks.emodnet
        assert len(net) > 0 and all(net.results[f'n_{k}'][1] > 0 for k in 'TIMC') # Relationships of every type after the first step (the warm-up is in it, after its stocks are recorded)
        assert net.burnin_left == 0
    if len(hiv.infections): # Seeding 2% of Risk HIGH may seed nobody at this size
        assert hiv.infections.type.isin([-2, -1, 0, 1, 2, 3]).all()


def test_births_and_nodes():
    """ Births seen by EMOD′ modules are the pregnancy's deliveries, and newborns are in their mother's node """
    sim = adapters.make_sim(pid=0, seed=3, n_agents=N, stop=1990, network='sti', hiv=False)
    sim.init()
    dem, n_born, n_check = sim.demographics.emoddemog, 0, 0
    for _ in range(60):
        sim.run_one_step()
        if len(dem.born):
            assert (sim.people.node.raw[dem.born] == sim.people.node.raw[dem.mothers]).all()
            assert (sim.people.parent.raw[dem.born] == dem.mothers).all()
            n_born += len(dem.born)
        n_check += sim.demographics.pregnancy.results.birth[sim.ti]
    assert n_born > 0 and np.isclose(n_born, n_check/sim.pop_scale) # The pregnancy's results are scaled


@pytest.mark.parametrize('network', ['sti', 'emod'])
def test_repeat(network):
    """ The same seed gives the same run """
    runs = [summary(adapters.make_sim(pid=0, seed=4, n_agents=N, stop=2000, network=network).run()) for _ in range(2)]
    assert np.array_equal(*runs)


def test_numba_stream():
    """ EmodNetwork's draws don't depend on other Numba code's use of the global stream, with own_stream (the default);
    without it they do, since EmodNetwork's base seeds that stream once, at init """
    def run(own, other):
        sim = adapters.make_sim(pid=0, seed=5, n_agents=N, stop=1995, network='emod', hiv=False)
        sim.init()
        sim.networks.emodnet.own_stream = own
        for _ in range(24):
            if other:
                other_numba_code(1000)
            sim.run_one_step()
        return summary(sim)
    assert np.array_equal(run(True, False), run(True, True))
    assert not np.array_equal(run(False, False), run(False, True))
