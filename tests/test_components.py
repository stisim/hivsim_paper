"""
Identity: each model's make_sim, rebuilt on its slot factories (models/*/components.py), gives exactly the runs it gave
before (at REF, the commit before the factories), same seed: every result array and every analyzer's table. And EMOD′'s
V as three modules (EmodS, EmodT, EmodH) gives exactly the runs of EmodHIV as one module (hiv.py at SPLIT), with every
cascade, the fixes, and on HIVsim's D and RN (EmodTOnSti), with EMOD's infection and death logs and the outputs' table.

    PYTHONPATH=<STIsim rs-v4-art-mortality>:<Starsim rs-nov16> pytest tests/test_components.py -n 8    # env hp_rc3; about 2 minutes on 8 cores

Each case runs in its own process (`python tests/test_components.py <case>`), since the two models have modules of the
same name (as they were at REF: cascade_host, and EMOD's cascade against common/cascade).
"""
import sys
import types
import subprocess
from pathlib import Path
import pytest

repo = Path(__file__).resolve().parents[1]
REF = '5df2972' # HIVsim migrated to rc3; harness.py, malawi.py and swap.py as they were before the factories
EMOD = repo / 'models' / 'emod' / 'v4'
HIVSIM = repo / 'models' / 'hivsim'
SPLIT = 'd7f950e' # EmodHIV as one module (transmission, natural history and seeding), before the split
EMOD_CASES = dict(network={}, hiv=dict(hiv=True), interpreted=dict(cascade=True), compiled=dict(cascade=True, compiled=True), condensed=dict(cascade=True, condensed=True))
HIVSIM_CASES = ['hivsim_malawi', 'hivsim_own', 'hivsim_condensed', 'hivsim_runner', 'hivsim_runner_emod_c']
SPLIT_CASES = dict(split_hiv=dict(hiv=True), split_interpreted=dict(cascade=True), split_compiled=dict(cascade=True, compiled=True),
                   split_condensed=dict(cascade=True, condensed=True), split_fixes=dict(cascade=True, condensed=True, fixes=('mortality', 'reenrol')),
                   split_hhee=None) # The last: the runner's D=hivsim,RN=hivsim,V=emod,C=emod, in a cell where the epidemic takes off


def at_ref(path, ref=REF, **replace):
    """ The module at `path` as it was at commit `ref` (with its own __file__, so its data paths hold), with some globals replaced """
    src = subprocess.check_output(['git', 'show', f'{ref}:{path.relative_to(repo)}'], cwd=repo, text=True)
    mod = types.ModuleType(path.stem + '_ref')
    mod.__file__ = str(path)
    exec(compile(src, str(path), 'exec'), mod.__dict__)
    mod.__dict__.update(replace)
    return mod


def art_mortality(sim):
    """ The pre-factory sim with mortality on ART by EMOD's table, which the factories have on by default since (models/hivsim/components.py, make_art) """
    import stisim as sti
    arts = [i for i in sim.pars.interventions if isinstance(i, sti.ART)]
    assert len(arts) == 1, arts
    arts[0].pars.art_death_table = sti.ART.emod_death_table
    return sim


def unsplit(comp, old):
    """ The factories `comp` with EmodHIV as one module (`old`, hiv.py at SPLIT), from make_H, as they had it at SPLIT """
    def make_H(inp, ctx):
        if not comp.with_hiv(ctx):
            return {}
        kw = dict(record=comp.opt(ctx, 'record', False), seeding=comp.seeding(inp, ctx))
        if comp.own(ctx, 'RN'):
            return dict(diseases=[old.EmodHIV(inp, **kw)])
        import adapters
        on_sti = type('EmodHIVOnSti', (old.EmodHIV,), dict(contacts=adapters.EmodTOnSti.contacts))
        return dict(diseases=[on_sti(inp, net='structuredsexual', **kw)], demographics=[adapters.EmodPeople(inp, risk_net='structuredsexual')])
    comp.make_S = comp.make_T = lambda inp, ctx: {}
    comp.make_H = make_H


def check(case):
    """ Run one case before and after; print the arrays that differ, and exit 1 if any do """
    sys.path[:0] = [str(repo), str(repo / 'common')]
    if case in SPLIT_CASES:
        sys.path.insert(0, str(EMOD))
        import harness
        import components
        old_hiv = at_ref(EMOD / 'hiv.py', ref=SPLIT)
        if case == 'split_hhee':
            from common import swaps
            choice, cfg = dict(D='hivsim', RN='hivsim', V='emod', C='emod'), dict(stop=2010, pairs=dict(hivsim=[0], emod=[0]), seeds=4)
            cell = swaps.cells(cfg)[3] # Seed 4: HIV takes off by 1995 (it doesn't with seeds 1-3)
            new = swaps.build(choice, cell, cfg)
            comp = swaps.runner.load_model('emod')
            unsplit(comp, old_hiv)
            old = swaps.build(choice, cell, cfg)
        else:
            kw = dict(seed=3, stop=2010) | SPLIT_CASES[case]
            new = harness.make_sim(0, **kw)
            unsplit(components, old_hiv)
            old = harness.make_sim(0, **kw)
    elif case in EMOD_CASES:
        sys.path.insert(0, str(EMOD))
        import harness
        import emod_cascade
        sys.modules['cascade'] = emod_cascade # Its name at REF
        ref = at_ref(EMOD / 'harness.py')
        kw = dict(seed=3, stop=2010) | EMOD_CASES[case]
        old, new = ref.make_sim(0, **kw), harness.make_sim(0, **kw)
    else:
        sys.path[:0] = [str(HIVSIM / 'malawi'), str(HIVSIM)]
        import malawi
        import swap
        import stisim_host
        sys.modules['cascade_host'] = stisim_host # Its name at REF
        ref_malawi = at_ref(HIVSIM / 'malawi' / 'malawi.py')
        if case == 'hivsim_runner': # The runner's endpoint: build(all HIVsim) is swap.make_sim('own'), with the outputs analyzer
            from common import swaps
            from common.outputs import Outputs
            cfg = dict(stop=2010, pairs=dict(hivsim=[0]))
            cell = swaps.cells(cfg)[0]
            old = swap.make_sim('own', seed=cell['seed'], stop=2010, analyzers=[Outputs(stisim_host.StisimHost)])
            new = swaps.build(dict(D='hivsim', RN='hivsim', V='hivsim', C='hivsim'), cell, cfg)
        elif case == 'hivsim_runner_emod_c': # The path's first step: EMOD′'s C on HIVsim is swap.make_sim('condensed') (EMOD's table set 0 on StisimHost)
            del sys.modules['cascade_host'] # EMOD′'s, from here
            from common import swaps
            from common.outputs import Outputs
            cfg = dict(stop=2010, pairs=dict(hivsim=[0], emod=[0]))
            cell = swaps.cells(cfg)[0]
            old = swap.make_sim('condensed', seed=cell['seed'], stop=2010, analyzers=[Outputs(stisim_host.StisimHost)])
            new = swaps.build(dict(D='hivsim', RN='hivsim', V='hivsim', C='emod'), cell, cfg)
        elif case == 'hivsim_malawi':
            old, new = art_mortality(ref_malawi.make_sim(seed=2, n_agents=3e3, verbose=0)), malawi.make_sim(seed=2, n_agents=3e3, verbose=0)
        else:
            ref = at_ref(HIVSIM / 'swap.py', malawi=ref_malawi)
            variant = case.removeprefix('hivsim_')
            old, new = art_mortality(ref.make_sim(variant, seed=2, n_agents=3e3)), swap.make_sim(variant, seed=2, n_agents=3e3)
    from common.identity import arrays, differences
    a, b = arrays(old.run()), arrays(new.run()) # Each initialised as it runs: EmodNetwork seeds Numba's global stream at its start
    if case in SPLIT_CASES:
        assert list(new.diseases.keys()) == ['emods', 'emodt', 'hiv'] and list(old.diseases.keys()) == ['hiv'], (new.diseases.keys(), old.diseases.keys())
    d = differences(a, b)
    print(f'{case}: {len(a)} arrays, {len(d)} differ')
    if len(d):
        print(d.head(20).to_string())
    sys.exit(int(len(d) > 0))


@pytest.mark.parametrize('case', list(EMOD_CASES) + HIVSIM_CASES + list(SPLIT_CASES))
def test_identity(case):
    r = subprocess.run([sys.executable, __file__, case], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    print(r.stdout.strip().splitlines()[-1])


if __name__ == '__main__':
    check(sys.argv[1])
