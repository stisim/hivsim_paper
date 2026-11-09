"""
EMOD′ burnt in to 1985, saved, and loaded and run on in a new process, equals the run from 1960.5 exactly (swap_design §1.4).

    pytest tests/test_resume.py    # env hp_rc3; about 2 minutes

The new process matters: EmodNetwork draws from Numba's generator, whose state is outside the sim (network.py saves it).
"""
import sys
import subprocess
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
EMOD = repo / 'models' / 'emod' / 'v4'
KW = dict(seed=4, stop=2010, cascade=True, condensed=True)


def setup():
    sys.path[:0] = [str(EMOD), str(repo)]
    import harness
    return harness


def burnin(fn):
    """ Run from 1960.5 and save its arrays; burn in to 1985 and save the sim """
    import numpy as np
    import sciris as sc
    harness = setup()
    from common.identity import arrays
    sc.save(f'{fn}.arrays', arrays(harness.make_sim(0, **KW).run()))
    sim = harness.make_sim(0, **KW)
    sim.init()
    years = np.asarray(sim.t.yearvec, dtype=float)
    sim.run(until=float(years[np.argmin(np.abs(years - 1985)) - 1]))
    sim.save(f'{fn}.sim')


def resume(fn):
    import sciris as sc
    import starsim as ss
    setup()
    from common.identity import arrays, differences
    d = differences(sc.load(f'{fn}.arrays'), arrays(ss.load(f'{fn}.sim').run()))
    print(f'{len(d)} arrays differ')
    if len(d):
        print(d.head(20).to_string())
    sys.exit(int(len(d) > 0))


def test_resume(tmp_path):
    fn = str(tmp_path / 'emod')
    for step in ('burnin', 'resume'): # Each in its own process
        r = subprocess.run([sys.executable, __file__, step, fn], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]


if __name__ == '__main__':
    dict(burnin=burnin, resume=resume)[sys.argv[1]](sys.argv[2])
