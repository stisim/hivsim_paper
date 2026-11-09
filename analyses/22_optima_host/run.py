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


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'label':
        run_label(sys.argv[2])
    elif cmd == 'sizes':
        sizes()
