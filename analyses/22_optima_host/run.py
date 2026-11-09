"""
Analysis 22: Optima's population label, and CascadeHost on Optima′'s agents (see README.md)

    python run.py label emod|hivsim   # Derived population sizes on EMOD′'s / HIVsim's people (HIVsim: PYTHONPATH=/tmp/sti-rs-v4-art)
    python run.py sizes               # Their comparison with Optima's
    python run.py outputs SEED [N]    # Optima′ as agents with Outputs(OptimaHost)
    python run.py compare             # Against analyses 19 and 20
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, sciris as sc
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'models' / 'optima' / 'v4'))
import poplabel as pl # noqa: E402

OUT = HERE / 'outputs'
YEARS = np.arange(1990, 2021)


def run_label(model):
    OUT.mkdir(exist_ok=True)
    T = sc.timer()
    if model == 'emod':
        sys.path.insert(0, str(REPO / 'models' / 'emod' / 'v4'))
        import components as cm
        sim = cm.make_sim(cm.pars(0), stop=2021, analyzers=pl.PopLabel(pl.emod_label, scale=cm.SCALE))
    else:
        sys.path.insert(0, str(REPO / 'models' / 'hivsim'))
        import components as cm
        sim = cm.make_sim(cm.pars(0), stop=2021, n_agents=100e3, sim_kwargs=dict(analyzers=pl.PopLabel(pl.stisim_label))) # 10x HIVsim's 10,000, for less noise in the key populations
    sim.run()
    sizes = sim.analyzers.poplabel.sizes
    sizes.to_csv(OUT / f'sizes_{model}.csv')
    print(f'{model}: {T.toc(output=True):.0f} s; pop_scale {sim.pop_scale}')


def optima_sizes():
    """ Optima's population sizes (its forced popsize, which the ODE follows) and Optima′ agents' (analysis 20, 100k agents, mean of seeds), years × pops """
    import optima_prime as opp
    ref = opp.load()
    pl.check(ref)
    sp = ref.simpars
    ode = pd.DataFrame(np.array([[np.interp(y, sp['tvec'], row) for row in sp['popsize']] for y in YEARS]), index=YEARS, columns=pl.POPS)
    files = sorted((REPO / 'analyses' / '20_optima_agents' / 'outputs').glob('abm_100000_*.obj'))
    ag = np.mean([sc.load(f).popsize for f in files], axis=0) # Pops × years (1990-2021)
    yrs = sc.load(files[0]).years
    agents = pd.DataFrame(ag[:, [list(yrs).index(y) for y in YEARS]].T, index=YEARS, columns=pl.POPS)
    return ode, agents, len(files)


def sizes():
    ode, agents, nseeds = optima_sizes()
    lines = [f'Optima′ agents: mean of {nseeds} seeds at 100k agents (analysis 20); Optima: its forced popsize (= the ODE)']
    lines.append('Optima′ agents / Optima − 1, % (2000, 2010, 2020), largest |%| over populations: ' +
                 ', '.join(f'{100*np.abs(agents.loc[y]/ode.loc[y] - 1).max():.2f}' for y in (2000, 2010, 2020)))
    tabs = {}
    for model in ('emod', 'hivsim'):
        f = OUT / f'sizes_{model}.csv'
        if not f.exists():
            continue
        d = pd.read_csv(f, index_col='year').reindex(YEARS)
        tabs[model] = d
        for y in (1990, 2000, 2010, 2020):
            o, m = ode.loc[y], d.loc[y]
            fem = [p for p in pl.POPS if p.startswith('F')]
            mal = [p for p in pl.POPS if p.startswith('M') and p != 'MSM']
            tab = pd.DataFrame({'optima': o, model: m, 'ratio': m/o,
                                'optima_share_of_sex_%': 100*o/np.where([p in fem for p in pl.POPS], o[fem].sum(), o[mal + ['MSM']].sum()),
                                f'{model}_share_of_sex_%': 100*m/np.where([p in fem for p in pl.POPS], m[fem].sum(), m[mal + ['MSM']].sum())})
            lines += ['', f'== {model}, {y} (people) ==', tab.round(3).to_string()]
    rows = []
    for y in YEARS:
        r = dict(year=y)
        for name, d in [('optima', ode)] + list(tabs.items()):
            f1549 = d.loc[y, [f'F{a}-{a+4}' for a in range(15, 50, 5)]].sum() + d.loc[y, 'FSW']
            m1549 = d.loc[y, [f'M{a}-{a+4}' for a in range(15, 50, 5)]].sum() + d.loc[y, ['Clients', 'MSM']].sum()
            r[f'{name}_fsw_%f15_49'] = 100*d.loc[y, 'FSW']/f1549
            r[f'{name}_clients_%m15_49'] = 100*d.loc[y, 'Clients']/m1549
            r[f'{name}_msm_%m15_49'] = 100*d.loc[y, 'MSM']/m1549
            r[f'{name}_total'] = d.loc[y].sum()
        rows.append(r)
    kp = pd.DataFrame(rows).set_index('year')
    kp.to_csv(OUT / 'key_pops.csv')
    lines += ['', '== Key populations, % of their sex aged 15-49, and total population ==', kp.iloc[::2].round(2).to_string()]
    (OUT / 'sizes.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


STATES = ['undx', 'dx', 'care', 'lost', 'usvl', 'svl']
FIG3A = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted']
PAPER = dict(undiagnosed='undiagnosed', diagnosed_naive='diagnosed_never_treated', art_lt6m='on_art', art_supp='on_art', art_nonsupp='on_art', interrupted='interrupted')
END = 2021


def run_outputs(seed, n=100_000):
    """ Optima′ as agents with Outputs(OptimaHost): yearly outputs, and each year's sexual transmissions by the source's state and category """
    import agents as ag, host as H
    from common.outputs import Outputs, yearly
    OUT.mkdir(exist_ok=True)
    T = sc.timer()
    sim = ag.make_sim(n_agents=n, rand_seed=seed, end=END, interventions=[H.OptimaTracker()], analyzers=[Outputs(H.OptimaHost)])
    sim.run()
    out = sim.analyzers.outputs
    y = yearly(out)
    srcs = out.host.sources
    yr = np.floor(np.concatenate([np.full(len(st), t) for t, st, _ in srcs]) + 1e-9).astype(int)
    st = pd.crosstab(yr, np.concatenate([s for _, s, _ in srcs])).reindex(columns=STATES, fill_value=0)
    cat = pd.crosstab(yr, np.concatenate([c for _, _, c in srcs])).reindex(columns=FIG3A, fill_value=0)
    r = sim.results.optima
    t = np.array([d.years for d in r['infection_susreg'].timevec])
    sexual = sum(np.asarray(r[k]) for k in ('infection_susreg', 'infection_progcirc'))
    res = pd.Series({yy: sexual[np.floor(t + 1e-9).astype(int) == yy].sum() for yy in st.index})
    sc.save(OUT / f'outputs_{n}_{seed}.obj', sc.objdict(yearly=y, by_state=st, by_category=cat, sexual_results=res, scale=out.scale, time=T.toc(output=True)))
    print(f'seed {seed}: {T.toc(output=True):.0f} s')


def run_ode():
    """ The agent formulation as an ODE: sexual infections per year by the source's state (analysis 19's attribution, summed over each year's steps) """
    import agents as ag, optima_prime as opp
    sim = ag.make_sim(solver='ode', end=END).run()
    by = opp.infections_by_source(sim)
    dt = sim.pars.dt.years
    yr = np.floor(by.index.to_numpy() + 1e-9).astype(int)
    st = by.T.groupby(by.columns.str.split('_').str[0]).sum().T.groupby(yr).sum()*dt
    by.to_csv(OUT / 'ode_by_source_steps.csv')
    st.reindex(columns=STATES, fill_value=0).to_csv(OUT / 'ode_by_state.csv') # svl: none (U=U in this project)


def compare():
    """ Agents (Outputs on OptimaHost) against the ODE (states), Optima (D8, analysis 19) and analysis 20's totals """
    files = sorted(OUT.glob('outputs_100000_*.obj'))
    runs = [sc.load(f) for f in files]
    yrs = np.arange(2000, 2021)
    lines = [f'{len(runs)} seeds at 100,000 agents; years {yrs[0]}-{yrs[-1]}; sexual transmissions summed over each year', '']
    sh = lambda d: 100*d.div(d.sum(axis=1), axis=0)
    # 1. By the source's state, against the ODE
    ode = pd.read_csv(OUT / 'ode_by_state.csv', index_col=0).loc[yrs]
    ag_st = np.array([sh(r.by_state.loc[yrs]).to_numpy() for r in runs])
    m, se = ag_st.mean(0), ag_st.std(0, ddof=1)/np.sqrt(len(runs))
    o = sh(ode).to_numpy()
    lines.append('== Share (%) by the source\'s state: agents mean (SE) vs ODE; z = (agents - ODE)/SE ==')
    tab = pd.DataFrame({f'{c}': [f'{m[i, j]:.2f} ({se[i, j]:.2f}) vs {o[i, j]:.2f}, z {(m[i, j] - o[i, j])/max(se[i, j], 1e-9):+.1f}' for i in range(len(yrs))] for j, c in enumerate(STATES)}, index=yrs)
    lines += [tab.iloc[::5].to_string(), tab.loc[[2020]].to_string()]
    z = (m - o)/np.maximum(se, 1e-9)
    lines.append(f'|z| > 2: {np.mean(np.abs(z) > 2)*100:.0f}% of state-years; mean |agents - ODE| {np.abs(m - o).mean():.2f} points, max {np.abs(m - o).max():.2f}')
    # 2. Figure 3A's categories (D8), agents' two ways (the host's ever-treated flag at the time, and the log's ART history), against analysis 19
    lines.append('')
    host = np.array([sh(r.by_category.loc[yrs]).to_numpy() for r in runs])
    tx = []
    for r in runs:
        y = r.yearly.set_index('year').loc[yrs]
        tx.append(np.array([100*sum(y[f'tx_{k}'] for k, c in PAPER.items() if c == cat) for cat in FIG3A]).T)
    tx = np.array(tx)
    a19 = pd.read_csv(REPO / 'analyses' / '19_fig3a_optima' / 'outputs' / 'optima_by_source.csv')
    a19 = a19[(a19.counted == 'sexual') & (a19.annual == 'step')]
    ref = {mp: sh(a19[a19.mapping == mp].pivot(index='year', columns='category', values='number').rename(columns={'interrupted_art': 'interrupted'})[FIG3A]) for mp in ('history', 'states')}
    pst = sh(ode.assign(dnt=ode.dx + ode.care, art=ode.usvl + ode.svl)[['undx', 'dnt', 'art', 'lost']].set_axis(FIG3A, axis=1))
    lines.append('== Figure 3A categories (%), D8: agents (host flag; transmission log), Optima′ ODE with lost = interrupted, Optima (analysis 19, step at the whole year: D8; lost = interrupted) ==')
    for yy in (2000, 2005, 2010, 2015, 2020):
        i = list(yrs).index(yy)
        lines.append(f'{yy}: ' + '; '.join(f'{c}: {host[:, i, j].mean():.1f}±{host[:, i, j].std(ddof=1)/np.sqrt(len(runs)):.1f}, {tx[:, i, j].mean():.1f} | ODE {pst.loc[yy, c]:.1f} | Optima {ref["history"].loc[yy, c]:.1f}, {ref["states"].loc[yy, c]:.1f}' for j, c in enumerate(FIG3A)))
    lines.append(f'Host flag vs transmission log: max difference {np.abs(host.mean(0) - tx.mean(0)).max():.2f} points')
    # 3. Totals, against analysis 20 (agents at 100k, ODE)
    lines.append('')
    a20 = REPO / 'analyses' / '20_optima_agents' / 'outputs'
    a20ode = sc.load(a20 / 'ode_agentpieces.obj')
    a20ag = [sc.load(f) for f in sorted(a20.glob('abm_100000_*.obj'))]
    Y = list(a20ode.years)
    idx = [Y.index(yy) for yy in yrs]
    tot = np.array([r.yearly.set_index('year').loc[yrs, 'tx_n'].to_numpy()*r.scale for r in runs])
    res = np.array([r.sexual_results.loc[yrs].to_numpy() for r in runs])
    a20m = np.mean([x.sexinci_year[idx] for x in a20ag], axis=0)
    lines.append('== Sexual infections per year: logged sources × scale (agents here), the disease\'s results (here), analysis 20 agents (48 seeds), ODE ==')
    for yy in (2000, 2005, 2010, 2015, 2018, 2020):
        i = list(yrs).index(yy)
        lines.append(f'{yy}: {tot[:, i].mean():,.0f} ± {tot[:, i].std(ddof=1)/np.sqrt(len(runs)):,.0f}; {res[:, i].mean():,.0f}; {a20m[i]:,.0f}; ODE {a20ode.sexinci_year[idx][i]:,.0f}')
    same = [np.allclose(r.sexual_results.loc[yrs].to_numpy(), sc.load(a20 / f.name.replace('outputs_', 'abm_')).sexinci_year[idx]) for r, f in zip(runs, files)]
    lines.append(f'The same seeds of analysis 20 (no host or Outputs): sexual infections identical every year in {sum(same)} of {len(same)} seeds')
    rel = tot.mean(0)/a20ode.sexinci_year[idx] - 1
    lines.append(f'Agents here / ODE - 1, mean over years {100*rel.mean():+.2f}%, max |%| {100*np.abs(rel).max():.2f} ({yrs[np.abs(rel).argmax()]}); logged = results exactly: {np.allclose(tot, res)}')
    (OUT / 'compare.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'label':
        run_label(sys.argv[2])
    elif cmd == 'sizes':
        sizes()
    elif cmd == 'outputs':
        run_outputs(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 100_000)
    elif cmd == 'ode':
        run_ode()
    elif cmd == 'compare':
        compare()
