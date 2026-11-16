"""
Analysis 26: Optima′ as agent modules against Optima′ declared as agents, and against the ODE (see README.md).

    python run.py one ARM SEED [N] [STOP]   # One run: ARM 'declared' (agents.make_sim) or 'modules' (components.make_sim)
    python run.py sweep CORES [SEEDS] [N]   # Every run not yet stored, both arms, seeds 1..SEEDS (64), in a pool
    python run.py ode                       # The ODE (optima_prime.make_sim with Optima's table), in the same measures
    python run.py runner                    # The runner's hybrid with every slot Optima′'s against components.make_sim (to 1995)
    python run.py tables [N]                # Tables: outputs/compare_N.txt, outputs/by_year_N.csv
"""
import os
for k in ('RAYON_NUM_THREADS', 'NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import sys
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OPT = REPO / 'models' / 'optima' / 'v4'
sys.path[:0] = [str(REPO), str(OPT)]
import starsim as ss # noqa: E402
import optima_prime as opp # noqa: E402

OUT, LOGS, RUNS = HERE / 'outputs', HERE / 'logs', HERE / 'outputs' / 'runs'
N_AGENTS, STOP, SEEDS = 100_000, 2021, 64
YEARS = np.arange(1990, STOP)
COUNTS = ['pop', 'plhiv', 'diagnosed', 'on_art', 'suppressed', 'new_infections', 'hiv_deaths', 'sexual']
FIG3A = ['undiagnosed', 'diagnosed_naive', 'on_art', 'interrupted'] # D8; on_art = art_lt6m + art_supp + art_nonsupp
PLACES = ['undx', 'dx', 'care', 'lost', 'svl', 'usvl']


class BandPlace(ss.Analyzer):
    """ People with HIV by place in care × CD4 band, at the step nearest each whole year (agents × people per agent) """
    def __init__(self, host, name='bandplace', **kwargs):
        super().__init__(name=name, **kwargs)
        self.host_class = host
        self.rows = {}
        return

    def init_post(self):
        super().init_post()
        self.host = self.host_class(self.sim)
        return

    def step(self):
        year = self.t.now('year')
        y = int(np.round(year))
        if abs(year - y) < 0.5*self.t.dt_year - 1e-9 and y not in self.rows:
            ppl = self.sim.people
            au = ppl.auids
            au = au[(ppl.ti_dead.raw[au] != self.ti) & (ppl.ti_removed.raw[au] != self.ti)]
            au = au[np.asarray(self.host.infected(au), bool)]
            p, b = np.asarray(self.host.care_place(au)), np.asarray(self.host.cd4_band(au))
            self.rows[y] = np.bincount((p - 1)*6 + b, minlength=36).reshape(6, 6)*self.sim.current_scale
        return


def declared(seed, n=N_AGENTS, stop=STOP):
    """ Optima′ declared, as agents, with Optima's table (as analysis 24's `optima`; PMTCT as components.PMTCT, which adds nothing here) """
    import agents as ag, host as H, components as oc
    from common.cascade.pathways import TableCascade
    from common.outputs import Outputs
    P = oc.pars()
    intvs = [TableCascade(P, H.OptimaHost, pmtct=oc.PMTCT), H.OptimaTracker()]
    return ag.make_sim(n_agents=n, rand_seed=seed, end=stop, cascade=P, interventions=intvs, analyzers=[BandPlace(H.OptimaHost), Outputs(H.OptimaHost)])


def modules(seed, n=N_AGENTS, stop=STOP):
    """ Optima′ as modules, from its slot factories (as the runner builds it) """
    import components as oc
    from modules import ModulesHost
    return oc.make_sim(seed=seed, stop=stop, n_agents=n, analyzers=[BandPlace(ModulesHost)])


def extract(sim):
    from common.outputs import yearly
    o = sim.analyzers.outputs
    y = yearly(o).set_index('year')
    out = pd.DataFrame(index=y.index)
    for c in COUNTS[:-1]:
        out[c] = y[c]*o.scale
    out['sexual'] = y.tx_n*o.scale
    out['undiagnosed'] = 100*y.tx_undiagnosed
    out['diagnosed_naive'] = 100*y.tx_diagnosed_naive
    out['on_art_tx'] = 100*(y.tx_art_lt6m + y.tx_art_supp + y.tx_art_nonsupp)
    out['interrupted'] = 100*y.tx_interrupted
    bp = sim.analyzers.bandplace.rows
    return sc.objdict(yearly=out, bandplace=np.array([bp[k] for k in sorted(bp)]), bp_years=np.array(sorted(bp)))


def run_one(arm, seed, n=N_AGENTS, stop=STOP):
    RUNS.mkdir(parents=True, exist_ok=True)
    T = sc.timer()
    sim = (declared if arm == 'declared' else modules)(seed, n, stop)
    sim.run()
    out = extract(sim)
    out.update(arm=arm, seed=seed, n_agents=n, time=T.toc(output=True), final_agents=len(sim.people.auids))
    sc.save(RUNS / f'{arm}_{n}_{seed}.obj', out)
    print(f'{arm} seed {seed}, {n} agents: {out.time:.0f} s, {out.final_agents} agents at the end', flush=True)
    return out


def ode_measures(sim):
    names = [s.replace('-', '_') for s in opp.load().settings['statelabels']]
    r = sim.results.optima
    t = np.array([d.years for d in r['susreg'].timevec])
    stock = lambda ns: sum(np.asarray(r[n]) for n in ns)
    at_end = lambda x: np.array([x[np.argmin(abs(t - (y + 1)))] for y in YEARS])
    hiv = [n for n in names if n not in ('susreg', 'progcirc')]
    out = pd.DataFrame(index=pd.Index(YEARS, name='year'))
    out['pop'] = at_end(stock(names))
    out['plhiv'] = at_end(stock(hiv))
    out['diagnosed'] = at_end(stock([n for n in hiv if not n.startswith('undx')]))
    out['on_art'] = at_end(stock([n for n in hiv if n.startswith(tuple(opp.TREATED))]))
    out['suppressed'] = at_end(stock([n for n in hiv if n.startswith('svl')]))
    yr = np.floor(t + 1e-9).astype(int)
    flow = lambda keys: np.array([sum(np.asarray(r[k]) for k in keys)[yr == y].sum() for y in YEARS])
    out['new_infections'] = flow([k for k in r.keys() if k.startswith('infection_') or k.startswith('birth_')])
    out['hiv_deaths'] = flow([k for k in r.keys() if k.startswith('hivdeath_')])
    out['sexual'] = flow([k for k in r.keys() if k.startswith('infection_')])
    by = opp.infections_by_source(sim)
    byr = np.floor(by.index.to_numpy() + 1e-9).astype(int)
    grp = by.T.groupby(by.columns.str.split('_').str[0]).sum().T.groupby(byr).sum().reindex(index=YEARS, columns=PLACES, fill_value=0)
    tot = grp.sum(axis=1)
    out['undiagnosed'] = 100*grp['undx']/tot
    out['diagnosed_naive'] = 100*(grp['dx'] + grp['care'])/tot
    out['on_art_tx'] = 100*(grp['usvl'] + grp['svl'])/tot
    out['interrupted'] = 100*grp['lost']/tot
    return out


def ode():
    import compare as cmp, components as oc
    OUT.mkdir(exist_ok=True)
    T = sc.timer()
    sim = opp.make_sim(components=cmp.COMPONENTS['E'], end=STOP + 1, cascade=oc.pars()).run()
    res = sc.objdict(yearly=ode_measures(sim), time=T.toc(output=True))
    sc.save(OUT / 'ode.obj', res)
    print(res.yearly.round(1).to_string())
    return res


def runner():
    """ The runner's hybrid with every slot Optima′'s against components.make_sim(), to 1995 at 100,000 agents """
    import components as oc
    from common import swaps
    cfg = dict(stop=1995, pairs=dict(optima=[None]), models=dict(optima='models/optima/v4/components.py'))
    cell = swaps.cells(cfg)[0]
    a = swaps.build(dict(D='optima', T='optima', H='optima', S='optima', C='optima'), cell, cfg).run()
    b = oc.make_sim(seed=cell['seed'], stop=1995).run()
    d = swaps.differences(a, b)
    text = f'runner build(D, T, H, S, C = optima) against components.make_sim(), seed {cell["seed"]}, 100,000 agents, 1990-1995: {len(d)} results differ {d}; sexual transmissions logged: {len(a.analyzers.outputs.df)} and {len(b.analyzers.outputs.df)}'
    print(text)
    OUT.mkdir(exist_ok=True)
    (OUT / 'runner_identity.txt').write_text(text + '\n')


def _job(job):
    arm, seed, n = job
    if (RUNS / f'{arm}_{n}_{seed}.obj').exists():
        return arm, seed, 'stored'
    LOGS.mkdir(exist_ok=True)
    with open(LOGS / f'{arm}_{n}_{seed}.log', 'w') as f:
        r = subprocess.run([sys.executable, __file__, 'one', arm, str(seed), str(n)], stdout=f, stderr=subprocess.STDOUT)
    return arm, seed, r.returncode


def sweep(cores, seeds=SEEDS, n=N_AGENTS):
    from multiprocessing import Pool
    jobs = [(a, s, n) for s in range(1, seeds + 1) for a in ('declared', 'modules')]
    with Pool(cores) as pool:
        for arm, seed, rc in pool.imap_unordered(_job, jobs):
            print(arm, seed, rc, flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'one':
        a = sys.argv[2:]
        run_one(a[0], int(a[1]), int(a[2]) if len(a) > 2 else N_AGENTS, int(a[3]) if len(a) > 3 else STOP)
    elif cmd == 'sweep':
        a = sys.argv[2:]
        sweep(int(a[0]), int(a[1]) if len(a) > 1 else SEEDS, int(a[2]) if len(a) > 2 else N_AGENTS)
    elif cmd == 'ode':
        ode()
    elif cmd == 'runner':
        runner()
    elif cmd == 'tables':
        import tables
        tables.main(int(sys.argv[2]) if len(sys.argv) > 2 else N_AGENTS)
