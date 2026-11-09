"""
Analysis 23: Optima's cascade as a table (see README.md)

    python run.py ode              # Optima′ as an ODE: its own cascade, and the table (common/cascade/tables/malawi_optima_params.csv)
    python run.py abm SEED [N]     # Optima′ as agents running the table through OptimaHost (common/cascade/pathways.py)
    python run.py compare [N]      # Tables: ODE table vs own; agents table vs analysis 20's agents with the own cascade
"""
import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd, sciris as sc
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'models' / 'optima' / 'v4'))
import optima_prime as opp, compare as cmp, agents as ag # noqa: E402
from common.cascade.params import Params # noqa: E402

spec = importlib.util.spec_from_file_location('a20', REPO / 'analyses' / '20_optima_agents' / 'run.py')
a20 = importlib.util.module_from_spec(spec); spec.loader.exec_module(a20)
OUT, A20 = HERE / 'outputs', REPO / 'analyses' / '20_optima_agents' / 'outputs'
TABLE = REPO / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv'
END = a20.END
YRS = np.arange(2000, 2021)
STATES = ['undx', 'dx', 'care', 'lost', 'usvl', 'svl']
FIG3A = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted']
KEYS = dict(plhiv='PLHIV', diag='Diagnosed', treat='On ART', supp='Suppressed', inci_year='New infections (year)', deaths_year='HIV deaths (year)')


def by_state(sim):
    """ Sexual infections per year by the source's cascade state (analysis 19's attribution, ODE) """
    by = opp.infections_by_source(sim)
    yr = np.floor(by.index.to_numpy() + 1e-9).astype(int)
    return by.T.groupby(by.columns.str.split('_').str[0]).sum().T.groupby(yr).sum().reindex(columns=STATES, fill_value=0)*sim.pars.dt.years


def run_ode():
    OUT.mkdir(exist_ok=True)
    for label, kw in [('ode_own', {}), ('ode_table', dict(cascade=Params(TABLE)))]:
        T = sc.timer()
        sim = opp.make_sim(components=cmp.COMPONENTS['E'], end=END, group_results=a20.group_results(), **kw).run()
        out = a20.extract(sim)
        out.by_state, out.time = by_state(sim), T.toc(output=True)
        sc.save(OUT / f'{label}.obj', out)
        print(f'{label}: {out.time:.1f} s')


def run_abm(seed, n=100_000):
    import host as H
    from common.cascade.pathways import TableCascade
    OUT.mkdir(exist_ok=True)
    T = sc.timer()
    sim = ag.make_sim(n_agents=n, rand_seed=seed, end=END, group_results=a20.group_results(), cascade=Params(TABLE),
                      interventions=[TableCascade(Params(TABLE), H.OptimaHost), H.OptimaTracker()])
    sim.run()
    out = a20.extract(sim)
    out.time, out.n_agents = T.toc(output=True), n
    sc.save(OUT / f'abm_{n}_{seed}.obj', out)
    print(f'abm seed {seed}, {n} agents: {out.time:.0f} s')


def shares(st):
    s = 100*st.div(st.sum(axis=1), axis=0)
    return pd.DataFrame({'undiagnosed': s.undx, 'diagnosed_never_treated': s.dx + s.care, 'on_art': s.usvl + s.svl, 'interrupted': s.lost})


def compare(n=100_000):
    Y = list(a20.YEARS)
    idx = [Y.index(y) for y in YRS]
    tot = lambda out, k: np.asarray(out[k]).sum(axis=0)[idx] if np.ndim(out[k]) == 2 else np.asarray(out[k])[idx]
    own, tab = sc.load(OUT / 'ode_own.obj'), sc.load(OUT / 'ode_table.obj')
    lines = ['== ODE: the table against Optima′\'s own cascade, 2000-2020 (table/own - 1) ==']
    for k, lab in KEYS.items():
        rel = tot(tab, k)/tot(own, k) - 1
        lines.append(f'{lab:24s} mean {100*rel.mean():+.2e}%, max |%| {100*np.abs(rel).max():.2e}; 2020: {tot(own, k)[-1]:,.0f} vs {tot(tab, k)[-1]:,.0f}')
    so, stb = shares(own.by_state.loc[YRS]), shares(tab.by_state.loc[YRS])
    lines.append(f'Figure 3A source shares (ODE mapping: care with dx, lost = interrupted): max |difference| {np.abs(so - stb).to_numpy().max():.2e} points; 2020 own ' +
                 ' / '.join(f'{x:.2f}' for x in so.loc[2020]) + '; table ' + ' / '.join(f'{x:.2f}' for x in stb.loc[2020]))
    files = sorted(OUT.glob(f'abm_{n}_*.obj'))
    if files:
        seeds = [int(f.stem.split('_')[-1]) for f in files]
        tb = [sc.load(f) for f in files]
        ow = [sc.load(A20 / f'abm_{n}_{s}.obj') for s in seeds]
        lines += ['', f'== Agents ({n:,}, seeds {seeds}): the table through OptimaHost against analysis 20\'s own cascade (same seeds); and each against the ODE ==',
                  'mean % = mean over years of (table mean / own mean - 1); z = difference / SE of the difference (independent runs)']
        for k, lab in KEYS.items():
            a, b = np.array([tot(x, k) for x in tb]), np.array([tot(x, k) for x in ow])
            se = np.sqrt(a.var(0, ddof=1)/len(a) + b.var(0, ddof=1)/len(b))
            z = (a.mean(0) - b.mean(0))/np.maximum(se, 1e-9)
            rel = a.mean(0)/b.mean(0) - 1
            o = tot(own, k)
            lines.append(f'{lab:24s} mean {100*rel.mean():+.2f}% (max |%| {100*np.abs(rel).max():.2f}, {YRS[np.abs(rel).argmax()]}); mean z {z.mean():+.2f}, |z| > 2 in {np.sum(np.abs(z) > 2)} of {len(YRS)} years; vs ODE: table {100*(a.mean(0)/o - 1).mean():+.2f}%, own {100*(b.mean(0)/o - 1).mean():+.2f}%')
        lines.append(f'Run time: table {np.mean([x.time for x in tb]):.0f} s, own {np.mean([x.time for x in ow]):.0f} s (analysis 20, 48 at once)')
    (OUT / f'compare_{n}.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ode':
        run_ode()
    elif cmd == 'abm':
        run_abm(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 100_000)
    elif cmd == 'compare':
        compare(int(sys.argv[2]) if len(sys.argv) > 2 else 100_000)
