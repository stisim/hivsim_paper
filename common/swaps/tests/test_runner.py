"""
The runner's tests (docs/swap_design.md 3, endpoint identity). EMOD′ ones take a few minutes; `-k toy` for the quick ones.

    conda run -n hp_rc3 python -m pytest common/swaps/tests -v
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

here = Path(__file__).resolve().parent
repo = here.parents[2]
sys.path[:0] = [str(repo / 'models' / 'emod' / 'v4'), str(repo)]
import harness # noqa: E402,F401 (first: pytest puts common/ on the path, where `cascade` is common/cascade, not EMOD′'s)
from common import swaps # noqa: E402

TOY = dict(models=dict(a='common/swaps/tests/toy_a.py', b='common/swaps/tests/toy_b.py'), stop=2003)
REAL = repo / 'models' / 'emod' / 'v4' / 'components.py'
EMOD = 'models/emod/v4/components.py' if REAL.exists() else 'common/swaps/tests/emod_standin.py'
STOP = 2005 # Short horizon: the cascade's first deliveries start in 2001
ALL = dict(D='emod', RN='emod', V='emod', C='emod')


def toy(**kw):
    return dict(TOY, slots=['V', 'N', 'D'], **{'from': 'a', 'to': 'b'}, **kw)


def test_toy_hybrids_and_cells():
    assert [swaps.name(h) for h in swaps.hybrids(toy())] == ['D=a,N=a,V=a', 'D=a,N=a,V=b', 'D=a,N=b,V=b', 'D=b,N=b,V=b']
    assert len(swaps.hybrids(toy(lattice=True))) == 8
    c = swaps.cells(dict(pairs=dict(a=[0, 1], b=[5, 6]), seeds=2))
    assert [(x['cell'], x['pids'], x['seed']) for x in c][-1] == ('p1_s2', dict(a=1, b=6), 102)


def test_toy_swap_identity():
    """ Swapping a slot for the same model's (as another name) changes nothing; swapping for the other model does """
    cfg = dict(TOY, models=dict(TOY['models'], a2=TOY['models']['a']))
    cell = swaps.cells(dict(pairs=dict(a=[1], a2=[1], b=[1])))[0]
    run = lambda choice: swaps.build(choice, cell, cfg).run()
    a = run(dict(D='a', N='a', V='a'))
    assert swaps.differences(a, run(dict(D='a', N='a2', V='a'))) == []
    assert swaps.differences(a, run(dict(D='a', N='b', V='a')))


def test_toy_run_store_resume(tmp_path):
    """ The pool runs every hybrid x cell once, stores them and the manifest, and a second call skips them all """
    cfg = toy(pairs=dict(a=[0, 1], b=[0, 1]), seeds=2)
    rows = swaps.run(cfg, out=tmp_path, cores=4)
    assert len(rows) == 4*4 and not any(r['error'] for r in rows)
    man = pd.read_csv(tmp_path / 'manifest.csv')
    assert len(man) == 16 and set(man.columns) >= {'hybrid', 'name', 'cell', 'seed', 'commit', 'starsim', 'seconds'}
    assert swaps.run(cfg, out=tmp_path) == []
    df = swaps.load(cfg, out=tmp_path)
    assert df.groupby(['name', 'cell']).ngroups == 16
    # Rerunning one cell by hand gives the stored table
    choice, cell = swaps.hybrids(cfg)[1], swaps.cells(cfg)[3]
    sim = swaps.build(choice, cell, cfg).run()
    stored = df[(df.name == swaps.name(choice)) & (df.cell == cell['cell'])].reset_index(drop=True)
    np.testing.assert_array_equal(stored['sis_I'].to_numpy(), sim.to_df()['sis_I'].to_numpy())


@pytest.fixture(scope='module')
def emod_ref():
    """ EMOD′ from its own harness, with the outputs analyzer """
    import harness
    from cascade_host import EmodHost
    from common.outputs import Outputs
    sim = harness.make_sim(0, seed=1, stop=STOP, cascade=True, condensed=True, record=False, analyzers=[Outputs(EmodHost, scale=100)])
    return sim.run()


def emod_cfg(**models):
    return dict(models=dict(emod=EMOD, **{m: EMOD for m in models}), stop=STOP, pairs=dict(emod=[0], **{m: [0] for m in models}))


def test_emod_endpoint(emod_ref):
    """ build(all EMOD′) equals EMOD′'s make_sim: every result, the reports, the outputs' counts and the transmission log """
    cfg = emod_cfg()
    sim = swaps.build(ALL, swaps.cells(cfg)[0], cfg).run()
    assert list(sim.analyzers.keys()) == list(emod_ref.analyzers.keys())
    assert swaps.differences(sim, emod_ref) == []
    assert len(sim.analyzers.outputs.df) > 1000 # Transmissions logged


def test_emod_trivial_swap(emod_ref):
    """ C and RN from 'emod_b' (EMOD′'s factories under another name): still EMOD′ exactly """
    cfg = emod_cfg(emod_b=True)
    sim = swaps.build(dict(ALL, C='emod_b', RN='emod_b'), swaps.cells(cfg)[0], cfg).run()
    assert swaps.differences(sim, emod_ref) == []


def test_emod_run_and_outputs(tmp_path):
    """ The pool's stored outputs for the EMOD′ endpoint: yearly outputs and targets measured """
    cfg = dict(emod_cfg(), hybrids=[ALL])
    rows = swaps.run(cfg, out=tmp_path, cores=1)
    assert len(rows) == 1 and not rows[0]['error'], rows
    y = swaps.load(cfg, out=tmp_path)
    assert {'tx_undiagnosed', 'share_deaths_on_art', 'incidence_15_49_f', 'new_infections'} <= set(y.columns)
    assert y.loc[y.year == 2003, 'tx_n'].item() > 0
    folder = next(p for p in tmp_path.iterdir() if p.is_dir())
    t = json.loads((folder / 'p0_s1.targets.json').read_text())
    assert np.isfinite(t['88']) and np.isnan(t['212']) # plhiv 1990.5; diagnosed 2020.5 is after the run


def test_toy_burnin(tmp_path, monkeypatch):
    """ With D's model's BURNIN, a run is burnt in, saved, and resumed from it the next time; both equal a continuous run """
    cfg = toy(pairs=dict(a=[0], b=[0]), hybrids=[dict(D='a', N='a', V='b')])
    monkeypatch.setattr(swaps.load_model('a', cfg['models']['a']), 'BURNIN', 2001.5, raising=False)
    choice, cell = cfg['hybrids'][0], swaps.cells(cfg)[0]
    ref = swaps.build(choice, cell, cfg).run().to_df()['sis_I'].to_numpy()
    for _ in range(2): # Burnt in, then resumed
        swaps.run(cfg, out=tmp_path, cores=1)
        folder = next(p for p in tmp_path.iterdir() if p.is_dir())
        assert (folder / 'p0_s1.burnin.sim').exists()
        np.testing.assert_array_equal(swaps.load(cfg, out=tmp_path)['sis_I'].to_numpy(), ref)
        (folder / 'p0_s1.csv.gz').unlink()


PATH_CFG = {'from': 'hivsim', 'to': 'emod', 'slots': ['C', 'V', 'RN', 'D'], 'pairs': {'hivsim': [0], 'emod': [0]}, 'stop': STOP}


def test_path_hybrids(tmp_path):
    """ The five hybrids of the first path (HIVsim → EMOD′: C, V, RN, D) build from the factories, with their bridges, and
    run through the pool to a short horizon (EMOD′'s D resumed from its burn-in), with finite outputs """
    hs = swaps.hybrids(PATH_CFG)
    rows = swaps.run(PATH_CFG, out=tmp_path, cores=len(hs))
    assert len(rows) == 5 and not any(r['error'] for r in rows), [r['error'] for r in rows]
    y = swaps.load(PATH_CFG, out=tmp_path)
    assert y.groupby('name').ngroups == 5
    for col in ('new_infections', 'hiv_deaths', 'incidence_15_49_f', 'incidence_15_49_m'):
        assert np.isfinite(y[col]).all(), col
    assert len(list(tmp_path.glob('*/p0_s1.burnin.sim'))) == 1 # All EMOD′ only
