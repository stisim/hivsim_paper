"""
Exp 03: run HIVsim Malawi (with testing by sex and over time) over a Latin hypercube of the ranges in config.yaml, and 10 seeds of the central set;
write each run's target values to outputs/results.jsonl as it finishes; then plot (models/hivsim/targetfit.py)

    nice python run.py
"""
import sys
import json
import traceback
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import sciris as sc
from scipy.stats import qmc

here = Path(__file__).parent
sys.path.insert(0, str(here.parents[1] / 'models' / 'hivsim'))
import malawi
import targetfit

cfg = sc.loadyaml(here / 'config.yaml')
outfile = here / 'outputs' / 'results.jsonl'


def to_value(u, low, high, scale):
    return float(np.exp(np.log(low) + u*(np.log(high) - np.log(low)))) if scale == 'log' else float(low + u*(high - low))


def run_one(job):
    kind, i, seed, pars = job
    try:
        hiv_pars = {k[4:]: v for k, v in pars.items() if k.startswith('hiv_')}
        nw_pars = {k[3:]: v for k, v in pars.items() if k.startswith('nw_')}
        other = {k: v for k, v in pars.items() if not k.startswith(('hiv_', 'nw_'))} # test_scale, test_rel_m, test_start, condom_scale
        T = sc.timer()
        sim = malawi.make_sim(seed=seed, n_agents=cfg['n_agents'], stop=cfg['stop'], verbose=0, analyzers=targetfit.TargetRecorder(),
                              hiv_pars=hiv_pars, nw_pars=nw_pars, **other).run()
        values = targetfit.measure(sim)
        years = np.asarray(sim.t.yearvec)
        on_art = np.asarray(sim.results.hiv.on_art)
        n_art = malawi.read('n_art.csv', index_col='year').n_art
        art_reached = {int(y): float(on_art[np.argmin(abs(years - y))] / n_art[y]) for y in (2010, 2015, 2020)}
        return dict(kind=kind, i=i, seed=seed, pars=pars, values={int(k): v for k, v in values.items()}, art_reached=art_reached, time=T.toc(output=True))
    except Exception:
        return dict(kind=kind, i=i, seed=seed, pars=pars, error=traceback.format_exc())


def jobs():
    names = list(cfg['pars'])
    u = qmc.LatinHypercube(d=len(names), seed=cfg['seed']).random(cfg['n_samples'])
    for i, row in enumerate(u):
        yield 'sample', i, i + 1, {n: to_value(x, *cfg['pars'][n]) for n, x in zip(names, row)}
    central = {n: to_value(0.5, *cfg['pars'][n]) for n in names}
    for s in range(cfg['n_seeds_central']):
        yield 'central', s, 1000 + s, central


if __name__ == '__main__':
    outfile.parent.mkdir(exist_ok=True)
    outfile.write_text('')
    alljobs = list(jobs())
    with Pool(min(len(alljobs), 150)) as pool, outfile.open('a') as f:
        for k, res in enumerate(pool.imap_unordered(run_one, alljobs)):
            f.write(json.dumps(res) + '\n')
            f.flush()
            if 'error' in res:
                print(f'Run {res["kind"]} {res["i"]} failed:\n{res["error"]}')
            print(f'{k+1}/{len(alljobs)} done', flush=True)
    targetfit.plot(here)
