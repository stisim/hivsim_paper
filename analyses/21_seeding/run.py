"""
Analysis 21: seeding options on the path HIVsim → EMOD′-fixed. See README.md.

    python run.py rates            # EMOD′-fixed's seeds per capita (50 seeds); HIVsim's seeding by EMOD's Risk class (20 seeds)
    python run.py runs [cores]     # option 0 (the 5 hybrids of the path) and the options' hybrids, 200 cells each
    python run.py bitcheck         # option 0's stored tables against analysis 16's (addendum), cell by cell
    python run.py analyse          # take-off and 2020 headline Δs per option and path
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))

OUT = here / 'outputs'
SEEDS = 200
E = 'emod_fixed'
V4 = 'models/emod/v4/'
MODELS = {E: V4 + 'components_fixed.py', 'emod_fixed_a': V4 + 'components_fixed_scaled.py', 'emod_fixed_c': V4 + 'components_fixed_seed5.py'}
TAKEOFF = 0.005 # Prevalence (all ages) at the end of 2020, as analysis 16
HEADLINE = ['new_infections', 'tx_undiagnosed', 'tx_diagnosed_naive', 'tx_on_art', 'tx_interrupted']


def h(D, RN, V, C, S=None):
    c = dict(D=D, RN=RN, V=V, C=C)
    if S:
        c['S'] = S
    return c


H = 'hivsim'
HYBRIDS = { # short name: choice
    'HHHH': h(H, H, H, H), 'HHHE': h(H, H, H, E), 'HHEE': h(H, H, E, E), 'HEEE': h(H, E, E, E), 'EEEE': h(E, E, E, E),
    'HHEE-a': h(H, H, 'emod_fixed_a', E), 'HEEE-a': h(H, E, 'emod_fixed_a', E),
    'HHEE-S': h(H, H, E, E, S=H), 'HEEE-S': h(H, E, E, E, S=H),
    'HHEE-c': h(H, H, 'emod_fixed_c', E), 'HEEE-c': h(H, E, 'emod_fixed_c', E), 'EEEE-c': h('emod_fixed_c', 'emod_fixed_c', 'emod_fixed_c', 'emod_fixed_c'),
}
CONFIG = {'hybrids': list(HYBRIDS.values()), 'pairs': {H: [0], E: [0], 'emod_fixed_a': [0], 'emod_fixed_c': [0]}, 'seeds': SEEDS, 'stop': 2025, 'models': MODELS}
PATHS = { # option: (steps, hybrids along the path)
    '0 (EMOD rule)': (['C', 'V', 'RN', 'D'], ['HHHH', 'HHHE', 'HHEE', 'HEEE', 'EEEE']),
    '(a) scaled': (['C', 'V', 'RN', 'D'], ['HHHH', 'HHHE', 'HHEE-a', 'HEEE-a', 'EEEE']),
    '(b) S after V': (['C', 'V', 'S', 'RN', 'D'], ['HHHH', 'HHHE', 'HHEE-S', 'HHEE', 'HEEE', 'EEEE']),
    '(b) S after RN': (['C', 'V', 'RN', 'S', 'D'], ['HHHH', 'HHHE', 'HHEE-S', 'HEEE-S', 'HEEE', 'EEEE']),
    '(b) S with D': (['C', 'V', 'RN', 'D+S'], ['HHHH', 'HHHE', 'HHEE-S', 'HEEE-S', 'EEEE']),
    '(c) 5x': (['C', 'V', 'RN', 'D'], ['HHHH', 'HHHE', 'HHEE-c', 'HEEE-c', 'EEEE']),
    '(c) 5x, EEEE 5x': (['C', 'V', 'RN', 'D'], ['HHHH', 'HHHE', 'HHEE-c', 'HEEE-c', 'EEEE-c']),
}


def a16():
    """ Analysis 16's run.py (its scales, intervals) """
    import importlib.util
    spec = importlib.util.spec_from_file_location('run16', here.parent / '16_hivsim_to_emod' / 'run.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def swaps():
    from common import swaps as s
    return s


# Rates
def emod_seeds(seed):
    """ EMOD′-fixed to its seeding step (1983): alive, HIGH and seeds there """
    sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
    import harness
    import hiv
    import components_fixed
    got = {}
    camp = hiv.EmodHIV.campaign
    def campaign(self):
        was = self.seeded
        camp(self)
        if self.seeded and not was:
            ppl = self.sim.people
            got.update(year=self.t.now('year'), alive=len(ppl.auids), high=int((ppl.risk.raw[ppl.auids] == 2).sum()), seeds=int(self.infected.raw[ppl.auids].sum()))
    hiv.EmodHIV.campaign = campaign
    harness.make_sim(0, seed=seed, stop=1984, network=True, hiv=True, record=False, fixes=components_fixed.FIXES).run()
    return dict(seed=seed, **got)


def hivsim_seeding(seed):
    """ HIVsim's 1985 population (HHHH's build, initialised): its seeding probability (init_prev_data, before rel_init_prev)
    among the active, by EMOD's Risk as the adapters read it, and sex; with HIGH and alive """
    s = swaps()
    hc = s.runner.load_model(H)
    sim = s.build(HYBRIDS['HHHH'], dict(cell=f'p0_s{seed}', pair=0, pids={H: 0}, seed=seed), dict(CONFIG, hybrids=[HYBRIDS['HHHH']]))
    sim.init()
    ppl, nw = sim.people, sim.networks.structuredsexual
    uids = ppl.auids
    data = hc.malawi.read('init_prev_hiv.csv').sort_values('sw')
    p = np.zeros(len(uids))
    active = nw.active(ppl)
    for r in data.itertuples():
        cond = ppl[r.sex] & (nw.risk_group == r.risk_group) & active
        if r.sw:
            cond = cond & (nw.fsw if r.sex == 'female' else nw.client)
        p[cond[uids]] = r.init_prev
    sw = nw.fsw.raw[uids] | nw.client.raw[uids]
    risk = np.where(sw, 2, np.where(nw.risk_group.raw[uids] == 0, 0, 1))
    act = active.raw[uids]
    rows = []
    for k in range(3):
        for f in (0, 1):
            m = act & (risk == k) & (ppl.female.raw[uids] == f)
            rows.append(dict(seed=seed, risk=k, female=f, n_active=int(m.sum()), mean_p=p[m].mean() if m.any() else np.nan))
    rel = hc.sub(hc.pars(0), 'hiv_').get('rel_init_prev', 1)
    extra = dict(seed=seed, alive=len(uids), high=int(sw.sum()), expected_hivsim=float((p*rel).sum()), expected_emod_rule=s.runner.load_model(E, MODELS[E]).pars(0).campaign.seed_coverage*int(sw.sum()), rel_init_prev=rel,
                 infected_at_init=int(sim.diseases.hiv.infected.raw[uids].sum()))
    return rows, extra


def rates():
    from multiprocessing import Pool
    OUT.mkdir(exist_ok=True)
    with Pool(50, maxtasksperchild=1) as pool:
        e = pd.DataFrame(pool.map(emod_seeds, range(1, 51)))
    with Pool(20, maxtasksperchild=1) as pool:
        hs = pool.map(hivsim_seeding, range(1, 21))
    sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
    import components_fixed
    cov = components_fixed.pars(0).campaign.seed_coverage
    e['per_capita'] = e.seeds/e.alive
    e.to_csv(OUT / 'emod_seeds_1983.csv', index=False)
    r = pd.DataFrame([x for rows, _ in hs for x in rows])
    x = pd.DataFrame([x for _, x in hs])
    r.to_csv(OUT / 'hivsim_seeding_by_risk.csv', index=False)
    x.to_csv(OUT / 'hivsim_seeding_1985.csv', index=False)
    # Pooled: E[seeds]/E[alive] for EMOD′; HIVsim's mean p by Risk and sex weighted by the active in each
    tab = r.assign(w=r.mean_p*r.n_active).groupby(['risk', 'female'])[['w', 'n_active']].sum()
    tab = (tab.w/tab.n_active).unstack('female')
    summ = dict(emod_seeds_mean=e.seeds.mean(), emod_alive_mean=e.alive.mean(), emod_high_mean=e.high.mean(), emod_per_capita=e.seeds.sum()/e.alive.sum(),
                emod_per_capita_expected=cov*e.high.sum()/e.alive.sum(), emod_year=e.year.iloc[0], seed_coverage=cov,
                hivsim_alive_mean=x.alive.mean(), hivsim_high_mean=x.high.mean(), hivsim_expected_seeds_emod_rule=x.expected_emod_rule.mean(),
                hivsim_expected_seeds_own_rule=x.expected_hivsim.mean(), hivsim_infected_at_init=x.infected_at_init.mean(),
                risk_prev=tab.to_numpy().round(6).tolist())
    summ['expected_seeds_scaled'] = summ['emod_per_capita_expected']*summ['hivsim_alive_mean']
    (OUT / 'rates.json').write_text(json.dumps(summ, indent=1))
    print(e.describe().round(5).to_string())
    print(json.dumps(summ, indent=1))


def runs():
    s = swaps()
    rows = s.run(CONFIG, cores=int(sys.argv[2]) if len(sys.argv) > 2 else 150)
    print(f'{len(rows)} runs, {sum(bool(r["error"]) for r in rows)} errors')
    errs = [r for r in rows if r['error']]
    for r in errs[:3]:
        print(r['name'], r['cell'], r['error'][-1500:])
    ids = {s.runner.hybrid_id(c, CONFIG, dict(s.runner.code_version(), **s.runner.pins())): k for k, c in HYBRIDS.items()}
    man = pd.read_csv(s.runner.OUT / 'manifest.csv')
    man = man[man.hybrid.isin(ids)].drop_duplicates(['hybrid', 'cell'], keep='last').assign(short=lambda d: d.hybrid.map(ids))
    t = man.groupby('short').agg(runs=('cell', 'size'), errors=('error', lambda e: int(e.fillna('').ne('').sum())), s_median=('seconds', 'median'),
                                  commit=('commit', 'last'), dirty=('dirty', 'last'), stisim=('stisim', 'last'), starsim=('starsim', 'last'))
    OUT.mkdir(exist_ok=True)
    t.to_csv(OUT / 'run_times.csv')
    print(t.to_string())


def short_of(df):
    key = {json.dumps(c, sort_keys=True): k for k, c in HYBRIDS.items()}
    slots = ['D', 'RN', 'V', 'S', 'C']
    return [key[json.dumps({k: v for k, v in zip(slots, r) if isinstance(v, str)}, sort_keys=True)] for r in df.reindex(columns=slots).itertuples(index=False)]


def table():
    s = swaps()
    df = s.load(CONFIG)
    df['short'] = short_of(df)
    sc_ = a16().scales()
    df['scale'] = np.where(df.D == H, sc_['hivsim'], 100.0)
    for c in ['plhiv', 'new_infections', 'hiv_deaths']:
        df[c + '_agents'] = df[c]
        df[c] = df[c]*df.scale
    df['tx_on_art'] = df.tx_art_lt6m + df.tx_art_supp + df.tx_art_nonsupp
    df['prev'] = df.plhiv_agents/df['pop']
    take = df[df.year == 2020].set_index(['short', 'cell']).prev.ge(TAKEOFF).rename('took_off')
    return df.join(take, on=['short', 'cell'])


def bitcheck():
    """ Option 0's stored tables against analysis 16's addendum (a /tmp worktree's raw_results), cell by cell """
    s = swaps()
    old = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('/tmp/hp-path/raw_results/swaps')
    olds = {}
    for f in old.glob('*/hybrid.json'):
        j = json.loads(f.read_text())
        if j.get('models', {}).get(E):
            olds[json.dumps(j['choice'], sort_keys=True)] = f.parent
    v = dict(s.runner.code_version(), **s.runner.pins())
    rows = []
    for k in ['HHHH', 'HHHE', 'HHEE', 'HEEE', 'EEEE']:
        c = HYBRIDS[k]
        new = s.runner.OUT / s.runner.hybrid_id(c, CONFIG, v)
        o = olds.get(json.dumps(c, sort_keys=True))
        same = n = 0
        for i in range(1, SEEDS + 1):
            a, b = new / f'p0_s{i}.csv.gz', (o / f'p0_s{i}.csv.gz') if o else None
            if a.exists() and b is not None and b.exists():
                n += 1
                x, y = pd.read_csv(a), pd.read_csv(b)
                same += int(x.shape == y.shape and list(x.columns) == list(y.columns) and np.array_equal(x.to_numpy(float), y.to_numpy(float), equal_nan=True))
        rows.append(dict(hybrid=k, old=str(o), compared=n, identical=same))
        print(rows[-1])
    OUT.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT / 'bitcheck.csv', index=False)


def analyse():
    m16 = a16()
    tci, wilson = m16.tci, m16.wilson
    df = table()
    t = df[df.year == 2020]
    t[['short', 'cell', 'prev', 'plhiv_agents', 'took_off']].to_csv(OUT / 'prev_2020_by_cell.csv', index=False)
    rows = []
    for k in HYBRIDS:
        x = t[t.short == k]
        if not len(x):
            continue
        kk, n = int(x.took_off.sum()), len(x)
        p, lo, hi = wilson(kk, n)
        rows.append(dict(hybrid=k, cells=n, took_off=kk, p=p, lo=lo, hi=hi, prev_2020_max_not=x[~x.took_off].prev.max(), prev_2020_min_taken_off=x[x.took_off].prev.min(),
                         prev_2020_mean_taken_off=x[x.took_off].prev.mean()))
    tk = pd.DataFrame(rows)
    tk.to_csv(OUT / 'takeoff.csv', index=False)
    print(tk.round(3).to_string())
    T = t.set_index(['short', 'cell'])
    rows, dp = [], []
    for opt, (steps, path) in PATHS.items():
        if not all(k in set(t.short) for k in path):
            continue
        allup = np.logical_and.reduce([T.loc[k].took_off.reindex(T.loc[path[0]].index).fillna(False).to_numpy(bool) for k in path])
        cells_all = T.loc[path[0]].index[allup]
        for i, step in enumerate(steps):
            a, b = path[i], path[i + 1]
            A, B = T.loc[a], T.loc[b].reindex(T.loc[a].index)
            n, m, lo, hi = tci(B.took_off.astype(float) - A.took_off.astype(float))
            dp.append(dict(option=opt, step=step, base=a, swap=b, p_base=A.took_off.mean(), p_swap=B.took_off.mean(), delta_p=m, lo=lo, hi=hi, both=int((A.took_off & B.took_off).sum()), path_all=len(cells_all)))
            both = A.took_off & B.took_off
            for o in HEADLINE:
                d = B[o] - A[o]
                for which, mask in (('all', np.ones(len(d), bool)), ('both', both.to_numpy(bool)), ('path', d.index.isin(cells_all))):
                    n, m, lo, hi = tci(d[mask])
                    rows.append(dict(option=opt, step=step, base=a, swap=b, output=o, which=which, n=n, delta=m, lo=lo, hi=hi))
    dl = pd.DataFrame(rows)
    dl.to_csv(OUT / 'headline_2020.csv', index=False)
    pd.DataFrame(dp).to_csv(OUT / 'takeoff_deltas.csv', index=False)
    print(pd.DataFrame(dp).round(3).to_string())
    # Means per hybrid (all, and taken off)
    rows = []
    for k in HYBRIDS:
        x = t[t.short == k]
        for o in HEADLINE + ['plhiv', 'prev']:
            n, m, lo, hi = tci(x[o])
            nc, mc, _, _ = tci(x.loc[x.took_off, o])
            rows.append(dict(hybrid=k, output=o, n=n, mean=m, lo=lo, hi=hi, n_took_off=nc, mean_took_off=mc))
    pd.DataFrame(rows).to_csv(OUT / 'means_2020.csv', index=False)
    # Wide table for the summary
    f = lambda r: f'{r.delta:+.3f} ({r.lo:+.3f}, {r.hi:+.3f}) n{r.n}' if r.output != 'new_infections' else f'{r.delta/1e3:+.1f}k ({r.lo/1e3:+.1f}, {r.hi/1e3:+.1f}) n{r.n}'
    dl['cell'] = [f(r) if np.isfinite(r.lo) else f'{r.delta:.3g} n{r.n}' for r in dl.itertuples()]
    w = dl.pivot_table(index=['option', 'which', 'output'], columns='step', values='cell', aggfunc='first', sort=False)
    w.to_csv(OUT / 'headline_2020_wide.csv')
    print(w.to_string())


if __name__ == '__main__':
    dict(rates=rates, runs=runs, bitcheck=bitcheck, analyse=analyse)[sys.argv[1]]()
