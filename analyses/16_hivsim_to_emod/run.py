"""
Analysis 16: the first swap path, HIVsim → EMOD′ (C, V, RN, D), uncompensated. See README.md.

    python run.py runs        # 5 hybrids x 200 cells, 1985-2025, into raw_results/swaps/
    python run.py identity    # endpoint identity at the full horizon, seed 1 (each end in its own process)
    python run.py analyse     # Δs per step (all runs, and conditional on take-off), take-off, waterfalls
"""
import sys
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))

OUT = here / 'outputs'
FIGS = here / 'figures'
SEEDS = 200
CONFIG = {'from': 'hivsim', 'to': 'emod', 'slots': ['C', 'V', 'RN', 'D'], 'pairs': {'hivsim': [0], 'emod': [0]}, 'seeds': SEEDS, 'stop': 2025}
YEARS = [2010, 2015, 2020, 2024]
COUNTS = ['plhiv', 'new_infections', 'hiv_deaths'] # Scaled to people
OUTPUTS = ['art_coverage', 'plhiv', 'new_infections', 'incidence_15_49_f', 'incidence_15_49_m', 'hiv_deaths', 'share_deaths_on_art',
           'tx_undiagnosed', 'tx_diagnosed_naive', 'tx_art_lt6m', 'tx_art_supp', 'tx_art_nonsupp', 'tx_interrupted']
TAKEOFF = 0.005 # Prevalence (all ages) at the end of 2020
STEPS = ['C', 'V', 'RN', 'D']


def swaps():
    from common import swaps as s
    return s


def short(choice):
    return ''.join('H' if choice[k] == 'hivsim' else 'E' for k in ('D', 'RN', 'V', 'C'))


def runs():
    s = swaps()
    rows = s.run(CONFIG, cores=int(sys.argv[2]) if len(sys.argv) > 2 else 100)
    print(f'{len(rows)} runs, {sum(bool(r["error"]) for r in rows)} errors')
    OUT.mkdir(exist_ok=True)
    man = pd.read_csv(s.runner.OUT / 'manifest.csv')
    ids = {s.runner.hybrid_id(h, CONFIG, dict(s.runner.code_version(), **s.runner.pins())): short(h) for h in s.hybrids(CONFIG)}
    man = man[man.hybrid.isin(ids)].drop_duplicates(['hybrid', 'cell'], keep='last').assign(short=lambda d: d.hybrid.map(ids))
    t = man.groupby('short', sort=False).agg(runs=('cell', 'size'), errors=('error', lambda e: int(e.fillna('').ne('').sum())),
                                             s_median=('seconds', 'median'), s_max=('seconds', 'max'), commit=('commit', 'last'), stisim=('stisim', 'last'), starsim=('starsim', 'last'))
    t.to_csv(OUT / 'run_times.csv')
    print(t.to_string())


# Endpoint identity
def identity_one(end):
    """ One endpoint at the full horizon, seed 1: the runner's build (EEEE resumed from its stored burn-in) against the model's own make_sim """
    if end == 'E':
        sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
        import harness # First, so `cascade` is EMOD′'s (as in common/swaps/tests)
    s = swaps()
    from common import outputs
    from common.outputs import Outputs
    cell = s.cells(CONFIG)[0]
    choice = dict(D='emod', RN='emod', V='emod', C='emod') if end == 'E' else dict(D='hivsim', RN='hivsim', V='hivsim', C='hivsim')
    folder = s.runner.OUT / s.runner.hybrid_id(choice, CONFIG, dict(s.runner.code_version(), **s.runner.pins()))
    new = s.runner.start(choice, cell, CONFIG, folder).run()
    if end == 'E':
        from cascade_host import EmodHost
        old = harness.make_sim(0, seed=cell['seed'], stop=CONFIG['stop'], cascade=True, condensed=True, record=False, analyzers=[Outputs(EmodHost, scale=100)]).run()
    else:
        sys.path[:0] = [str(repo / 'models' / 'hivsim' / 'malawi'), str(repo / 'models' / 'hivsim')]
        import swap
        import stisim_host
        old = swap.make_sim('own', seed=cell['seed'], stop=CONFIG['stop'], analyzers=[Outputs(stisim_host.StisimHost)]).run()
    d = s.differences(new, old)
    stored = pd.read_csv(folder / f'{cell["cell"]}.csv.gz')
    y = outputs.yearly(new.analyzers.outputs)
    same = stored.shape == y.shape and np.allclose(stored.to_numpy(float), y.to_numpy(float), equal_nan=True, rtol=1e-12, atol=0)
    print(json.dumps(dict(end=end, cell=cell['cell'], stop=CONFIG['stop'], differing=len(d), first=[str(x) for x in d[:5]], stored_reproduced=bool(same),
                          transmissions=len(new.analyzers.outputs.df))))


def identity():
    OUT.mkdir(exist_ok=True)
    rows = []
    for end in ('H', 'E'):
        r = subprocess.run([sys.executable, __file__, 'identity_one', end], capture_output=True, text=True)
        line = [l for l in r.stdout.splitlines() if l.startswith('{')]
        rows.append(json.loads(line[-1]) if line else dict(end=end, error=(r.stderr or r.stdout)[-2000:]))
        print(rows[-1])
    pd.DataFrame(rows).to_csv(OUT / 'identity.csv', index=False)


# Analysis
def scales():
    """ People per agent, by D's model: EMOD′'s SCALE; HIVsim's sim.pop_scale (from its total_pop at the start) """
    s = swaps()
    cell = s.cells(CONFIG)[0]
    sim = s.build(dict(D='hivsim', RN='hivsim', V='hivsim', C='hivsim'), cell, CONFIG)
    sim.init()
    return dict(hivsim=float(sim.pop_scale), emod=float(s.load_model('emod').SCALE))


def table():
    s = swaps()
    df = s.load(CONFIG)
    sc_ = scales()
    df['short'] = [''.join('H' if v == 'hivsim' else 'E' for v in r) for r in df[['D', 'RN', 'V', 'C']].itertuples(index=False)]
    df['scale'] = df.D.map(sc_)
    for c in COUNTS:
        df[c + '_agents'] = df[c]
        df[c] = df[c]*df.scale
    df['prev'] = df.plhiv_agents/df['pop']
    take = df[df.year == 2020].set_index(['short', 'cell']).prev.ge(TAKEOFF).rename('took_off')
    df = df.join(take, on=['short', 'cell'])
    return df, sc_


def tci(x):
    """ Mean and 95% t-interval """
    from scipy import stats
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 2:
        return n, (x.mean() if n else np.nan), np.nan, np.nan
    h = stats.t.ppf(0.975, n - 1)*x.std(ddof=1)/np.sqrt(n)
    return n, x.mean(), x.mean() - h, x.mean() + h


def wilson(k, n):
    z = 1.96
    p = k/n
    c = (p + z*z/(2*n))/(1 + z*z/n)
    h = z*np.sqrt(p*(1 - p)/n + z*z/(4*n*n))/(1 + z*z/n)
    return p, c - h, c + h


def analyse():
    OUT.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)
    df, sc_ = table()
    path = ['HHHH', 'HHHE', 'HHEE', 'HEEE', 'EEEE']
    print('scales', sc_, 'runs', df.groupby('short').cell.nunique().to_dict())

    # (c) take-off
    t = df[df.year == 2020]
    rows = []
    for h in path:
        x = t[t.short == h]
        k, n = int(x.took_off.sum()), len(x)
        p, lo, hi = wilson(k, n)
        rows.append(dict(hybrid=h, cells=n, took_off=k, p=p, lo=lo, hi=hi, any_plhiv=int((x.plhiv_agents > 0).sum()),
                         prev_2020_median_taken_off=x[x.took_off].prev.median(), prev_2020_max_not=x[~x.took_off].prev.max(),
                         prev_2020_min_taken_off=x[x.took_off].prev.min()))
    tk = pd.DataFrame(rows)
    tk.to_csv(OUT / 'takeoff.csv', index=False)
    print(tk.round(4).to_string())
    t[['short', 'cell', 'prev', 'plhiv_agents', 'took_off']].to_csv(OUT / 'prev_2020_by_cell.csv', index=False)

    # Means per hybrid
    d = df[df.year.isin(YEARS)]
    rows = []
    for h in path:
        for y in YEARS:
            x = d[(d.short == h) & (d.year == y)]
            for o in OUTPUTS:
                n, m, lo, hi = tci(x[o])
                nc, mc, loc, hic = tci(x.loc[x.took_off, o])
                rows.append(dict(hybrid=h, year=y, output=o, n=n, mean=m, lo=lo, hi=hi, n_cond=nc, mean_cond=mc, lo_cond=loc, hi_cond=hic))
    means = pd.DataFrame(rows)
    means.to_csv(OUT / 'means.csv', index=False)

    # (a), (b) Δs per step, paired by cell
    rows = []
    for i, step in enumerate(STEPS):
        a, b = path[i], path[i + 1]
        A = d[d.short == a].set_index(['cell', 'year'])
        B = d[d.short == b].set_index(['cell', 'year'])
        both = A.took_off & B.took_off.reindex(A.index)
        pa, pb = tk.set_index('hybrid').p[a], tk.set_index('hybrid').p[b]
        for y in YEARS:
            for o in OUTPUTS:
                dd = (B[o] - A[o].reindex(B.index)).xs(y, level='year')
                cond = both.xs(y, level='year').reindex(dd.index).fillna(False).astype(bool)
                n, m, lo, hi = tci(dd)
                nc, mc, loc, hic = tci(dd[cond])
                rows.append(dict(step=step, base=a, swap=b, year=y, output=o, n=n, delta=m, lo=lo, hi=hi,
                                 n_cond=nc, delta_cond=mc, lo_cond=loc, hi_cond=hic, cond_minus_all=mc - m,
                                 base_mean=A[o].xs(y, level='year').mean(), takeoff_base=pa, takeoff_swap=pb))
    deltas = pd.DataFrame(rows)
    deltas.to_csv(OUT / 'deltas.csv', index=False)
    # Take-off change per step, paired: share of cells taking off in swap minus base, with a t-interval of the paired indicator difference
    rows = []
    for i, step in enumerate(STEPS):
        a, b = path[i], path[i + 1]
        A = t[t.short == a].set_index('cell').took_off.astype(float)
        B = t[t.short == b].set_index('cell').took_off.astype(float).reindex(A.index)
        n, m, lo, hi = tci(B - A)
        rows.append(dict(step=step, base=a, swap=b, n=n, delta_p=m, lo=lo, hi=hi, both=int(((A > 0) & (B > 0)).sum())))
    dt = pd.DataFrame(rows)
    dt.to_csv(OUT / 'takeoff_deltas.csv', index=False)
    print(dt.round(3).to_string())

    # Headline table (2020, 2024)
    h = deltas[deltas.year.isin([2020, 2024])].copy()
    fmt = lambda m, lo, hi, n: f'{m:.4g} ({lo:.3g}, {hi:.3g}) n={n}' if np.isfinite(lo) else f'{m:.4g} n={n}'
    h['all'] = [fmt(*r) for r in h[['delta', 'lo', 'hi', 'n']].itertuples(index=False)]
    h['conditional'] = [fmt(*r) for r in h[['delta_cond', 'lo_cond', 'hi_cond', 'n_cond']].itertuples(index=False)]
    h[['output', 'year', 'step', 'base_mean', 'all', 'conditional']].to_csv(OUT / 'deltas_2020_2024.csv', index=False)
    waterfalls(deltas, means, path)


def waterfalls(deltas, means, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels = ['HHHH', '+C', '+V', '+RN', '+D', 'EEEE']
    for o in OUTPUTS:
        fig, axs = plt.subplots(1, len(YEARS), figsize=(16, 4), sharey=True)
        for ax, y in zip(axs, YEARS):
            m = means[(means.output == o) & (means.year == y)].set_index('hybrid')
            dl = deltas[(deltas.output == o) & (deltas.year == y)].set_index('step').loc[STEPS]
            for k, (col, lo, hi, colr, off, lab) in enumerate([('delta', 'lo', 'hi', '#4C72B0', -0.18, 'all runs'), ('delta_cond', 'lo_cond', 'hi_cond', '#DD8452', 0.18, 'both arms took off')]):
                start = m['mean' if k == 0 else 'mean_cond']['HHHH']
                ax.bar(0 + off, start, width=0.34, color=colr, alpha=0.5)
                cum = start
                for j, step in enumerate(STEPS):
                    v = dl[col][step]
                    if not np.isfinite(v):
                        continue
                    ax.bar(j + 1 + off, v, bottom=cum, width=0.34, color=colr, label=lab if j == 0 else None)
                    ax.errorbar(j + 1 + off, cum + v, yerr=[[v - dl[lo][step]], [dl[hi][step] - v]], color='k', lw=1, capsize=2)
                    cum += v
                end = m['mean' if k == 0 else 'mean_cond']['EEEE']
                ax.bar(5 + off, end, width=0.34, color=colr, alpha=0.5)
            ax.axhline(0, color='grey', lw=0.5)
            ax.set_xticks(range(6))
            ax.set_xticklabels(labels)
            ax.set_title(f'{o}, {y}')
        axs[0].legend(fontsize=8)
        fig.suptitle(f'{o}: HIVsim → EMOD′ (D, RN, V, C), uncompensated; endpoints are hybrid means, bars the paired Δ per step (95% CI)', fontsize=9)
        fig.tight_layout()
        fig.savefig(FIGS / f'waterfall_{o}.png', dpi=110)
        plt.close(fig)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'identity_one':
        identity_one(sys.argv[2])
    else:
        dict(runs=runs, identity=identity, analyse=analyse)[cmd]()
