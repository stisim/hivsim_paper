"""
Analysis 23: Optima's cascade as a table (see README.md)

    python run.py ode              # Optima′ as an ODE: its own cascade, and the table (common/cascade/tables/malawi_optima_params.csv)
    python run.py abm SEED [N]     # Optima′ as agents running the table through OptimaHost (common/cascade/pathways.py)
    python run.py outputs SEED [N] # As abm, with Outputs(OptimaHost): sources by Figure 3A category
    python run.py fig3a            # Their shares against analysis 22's (own cascade)
    python run.py compare [N]      # Tables: ODE table vs own; agents table vs analysis 20's agents with the own cascade
    PYTHONPATH=/tmp/sti-613:/tmp/ss-nov09 python run.py smoke    # The table on HIVsim's people (runs; not interpreted)
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


def run_outputs(seed, n=100_000):
    """ As run_abm, with Outputs(OptimaHost): each year's sexual transmissions by the source's Figure 3A category (as analysis 22's, with Optima's own cascade) """
    import host as H
    from common.cascade.pathways import TableCascade
    from common.outputs import Outputs
    T = sc.timer()
    P = Params(TABLE)
    sim = ag.make_sim(n_agents=n, rand_seed=seed, end=END, cascade=P, interventions=[TableCascade(P, H.OptimaHost), H.OptimaTracker()], analyzers=[Outputs(H.OptimaHost)])
    sim.run()
    srcs = sim.analyzers.outputs.host.sources
    yr = np.floor(np.concatenate([np.full(len(st), t) for t, st, _ in srcs]) + 1e-9).astype(int)
    cat = pd.crosstab(yr, np.concatenate([c for _, _, c in srcs])).reindex(columns=FIG3A, fill_value=0)
    sc.save(OUT / f'outputs_{n}_{seed}.obj', sc.objdict(by_category=cat, time=T.toc(output=True)))
    print(f'outputs seed {seed}: {T.toc(output=True):.0f} s')


def compare_fig3a(n=100_000):
    """ Figure 3A source shares (D8) as agents: the table against analysis 22's own cascade (16 seeds), and the ODE """
    sh = lambda d: 100*d.div(d.sum(axis=1), axis=0)
    tb = np.array([sh(sc.load(f).by_category.loc[YRS]).to_numpy() for f in sorted(OUT.glob(f'outputs_{n}_*.obj'))])
    ow = np.array([sh(sc.load(f).by_category.loc[YRS]).to_numpy() for f in sorted((REPO / 'analyses' / '22_optima_host' / 'outputs').glob(f'outputs_{n}_*.obj'))])
    ode = shares(sc.load(OUT / 'ode_own.obj').by_state.loc[YRS]).to_numpy()
    se = np.sqrt(tb.var(0, ddof=1)/len(tb) + ow.var(0, ddof=1)/len(ow))
    z = (tb.mean(0) - ow.mean(0))/np.maximum(se, 1e-9)
    lines = [f'== Figure 3A source shares (%), agents: table ({len(tb)} seeds) vs own cascade (analysis 22, {len(ow)} seeds), and the ODE (own cascade; lost = interrupted) ==']
    for y in (2000, 2005, 2010, 2015, 2020):
        i = list(YRS).index(y)
        lines.append(f'{y}: ' + '; '.join(f'{c}: {tb[:, i, j].mean():.1f} vs {ow[:, i, j].mean():.1f} (z {z[i, j]:+.1f}), ODE {ode[i, j]:.1f}' for j, c in enumerate(FIG3A)))
    lines.append(f'Over 2000-2020: mean |table - own| {np.abs(tb.mean(0) - ow.mean(0)).mean():.2f} points, max {np.abs(tb.mean(0) - ow.mean(0)).max():.2f}; |z| > 2 in {np.mean(np.abs(z) > 2)*100:.0f}% of category-years')
    (OUT / f'fig3a_{n}.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def smoke_hivsim(n=10_000):
    """ The table on HIVsim's people (StisimHost, which has no care or lost of its own and keeps the ART arm): it runs; not interpreted """
    sys.path.insert(0, str(REPO / 'models' / 'hivsim'))
    import components as hc, starsim as ss, stisim as sti
    T = sc.timer()
    sim = hc.make_sim(hc.pars(0), n_agents=n, stop=2021, cascade='table', params=Params(TABLE)).run()
    tc = sim.interventions.table_cascade
    au = sim.people.auids
    place = np.bincount(tc.places(au), minlength=7)*sim.current_scale
    lines = [f'HIVsim ({n:,} agents; Starsim {ss.__file__}, STIsim {sti.__file__}) with Optima\'s table, 2021: {T.toc(output=True):.0f} s',
             'people by place: ' + ', '.join(f'{p}: {x:,.0f}' for p, x in zip(['uninfected', 'undx', 'dx', 'care', 'lost', 'svl', 'usvl'], place)),
             f'target on ART {Params(TABLE).value("target", "on_art", sim.t.now("year")):,.0f}']
    (OUT / 'smoke_hivsim.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ode':
        run_ode()
    elif cmd == 'abm':
        run_abm(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 100_000)
    elif cmd == 'compare':
        compare(int(sys.argv[2]) if len(sys.argv) > 2 else 100_000)
    elif cmd == 'outputs':
        run_outputs(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 100_000)
    elif cmd == 'fig3a':
        compare_fig3a()
    elif cmd == 'smoke':
        smoke_hivsim()
