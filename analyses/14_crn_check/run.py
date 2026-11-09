"""
Analysis 14: common random numbers on step 1 (C). See README.md.

    python run.py runs      # the three hybrids x 60 seeds, through the runner (raw_results/swaps/)
    python run.py diverge   # 20 seeds, agent-by-agent divergence of paired runs (outputs/divergence.csv)
    python run.py rounds    # seed 1: rounds of sti.HIV's engine per step, base vs swap (outputs/engine_rounds_seed1.csv)
    python run.py analyse   # tables (outputs/) and figures (figures/)
"""
import os
import sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))
from common import swaps # noqa: E402

OUT, FIG = here / 'outputs', here / 'figures'
H = dict(D='hivsim', RN='hivsim', V='hivsim')
HYBRIDS = dict(base=H | dict(C='hivsim'), swap=H | dict(C='hivsim_cx'), control=H | dict(C='hivsim_b'))
CONFIG = dict(hybrids=list(HYBRIDS.values()), pairs={'hivsim': [0], 'hivsim_cx': [0], 'hivsim_b': [0]}, seeds=60, stop=2025,
              models={'hivsim_cx': 'analyses/14_crn_check/hivsim_cx.py', 'hivsim_b': 'models/hivsim/components.py'})
NDIV = 20 # Seeds for the divergence runs
YEARS = [1995, 2000, 2005, 2010, 2015, 2020, 2024]
TX = ['tx_undiagnosed', 'tx_diagnosed_naive', 'tx_art_lt6m', 'tx_art_supp', 'tx_art_nonsupp', 'tx_interrupted']
MEASURES = ['art_coverage', 'plhiv', 'new_infections', 'hiv_deaths', 'share_deaths_on_art', 'incidence_15_49_f', 'incidence_15_49_m'] + TX
BITS = dict(alive=0, infected=1, diagnosed=2, on_art=3) # Bits of each starting agent's recorded state


def runs():
    rows = swaps.run(CONFIG, cores=160)
    print(f'{len(rows)} runs, {sum(bool(r["error"]) for r in rows)} errors')


class AgentStates:
    """ Each step, for the starting agents: alive, infected, diagnosed and on ART (as bits), and number of sexual partners. Draws nothing """
    def __init__(self):
        import starsim as ss
        class _A(ss.Analyzer):
            def init_post(a):
                super().init_post()
                a.n0 = len(a.sim.people.auids)
                a.state, a.partners, a.years = [], [], []
            def step(a):
                sim = a.sim
                ppl, hiv, host = sim.people, sim.diseases.hiv, sim.analyzers.outputs.host
                u = ss.uids(np.arange(a.n0))
                alive = ppl.alive.raw[u]
                b = alive.astype(np.uint8)
                st = host.cascade_stage(u)
                b |= (host.infected(u) & alive).astype(np.uint8) << BITS['infected']
                b |= (st['diagnosed'] & alive).astype(np.uint8) << BITS['diagnosed']
                b |= (st['on_art'] & alive).astype(np.uint8) << BITS['on_art']
                e = sim.networks.structuredsexual.edges
                n = np.bincount(np.concatenate([e.p1, e.p2]).astype(int), minlength=len(ppl.uid.raw))[:a.n0]
                a.state.append(b)
                a.partners.append(np.minimum(n, 255).astype(np.uint8))
                a.years.append(a.t.now('year'))
        self.cls = _A


def diverge_one(seed):
    """ The three hybrids for one seed, one after another; per step, the share of starting agents whose state differs from base """
    from common import outputs
    cell = next(c for c in swaps.cells(CONFIG) if c['seed'] == seed)
    rec, yearly = {}, []
    for k, choice in HYBRIDS.items():
        sim = swaps.build(choice, cell, CONFIG)
        sim.pars.analyzers = list(sim.pars.analyzers) + [AgentStates().cls(name='agentstates')]
        sim.run()
        a = sim.analyzers.agentstates
        rec[k] = (np.array(a.state), np.array(a.partners), np.array(a.years))
        yearly.append(outputs.yearly(sim.analyzers.outputs).assign(arm=k, seed=seed))
    s0, p0, t = rec['base']
    rows = []
    for k in ('swap', 'control'):
        s1, p1, _ = rec[k]
        d = dict(year=t, arm=k, seed=seed, any_state=(s0 != s1).mean(1), partners=(p0 != p1).mean(1))
        for name, bit in BITS.items():
            d[name] = (((s0 >> bit) & 1) != ((s1 >> bit) & 1)).mean(1)
        rows.append(pd.DataFrame(d))
    return pd.concat(rows), pd.concat(yearly)


def diverge():
    with Pool(NDIV) as pool:
        res = pool.map(diverge_one, range(1, NDIV + 1))
    pd.concat([r[0] for r in res]).to_csv(OUT / 'divergence.csv.gz', index=False)
    pd.concat([r[1] for r in res]).to_csv(OUT / 'divergence_yearly.csv.gz', index=False)


def rounds(seed=1, stop=1993):
    """ Where pairing first breaks (seed 1): per step, the number of rounds of moves in sti.HIV's agent engine (its draws
    are keyed by step and round) and the first infections whose drawn acute duration differs between base and swap """
    import starsim as ss
    cfg = dict(CONFIG, stop=2005) # The condensed cascade has rounds from 2001, and Starsim rc3 won't start one after the stop
    cell = next(c for c in swaps.cells(cfg) if c['seed'] == seed)
    class Rounds(ss.Analyzer):
        def init_post(a):
            super().init_post()
            a.rows = []
        def step(a):
            if a.t.now('year') <= stop:
                hiv = a.sim.diseases.hiv
                new = (hiv.ti_infected == a.ti).uids
                a.rows.append(dict(ti=a.ti, year=a.t.now('year'), rounds=hiv.engine.round[1] + 1, new=';'.join(map(str, new)),
                                   exit_acute=';'.join(f'{x:g}' for x in hiv.ti_exit_acute.raw[new])))
    out = []
    for k in ('base', 'swap'):
        sim = swaps.build(HYBRIDS[k], cell, cfg)
        sim.pars.analyzers = list(sim.pars.analyzers) + [Rounds(name='rounds')]
        sim.run()
        out.append(pd.DataFrame(sim.analyzers.rounds.rows).set_index(['ti', 'year']).add_suffix(f'_{k}'))
    df = pd.concat(out, axis=1).reset_index()
    df['same_draws'] = (df.new_base == df.new_swap) & (df.exit_acute_base == df.exit_acute_swap)
    df.to_csv(OUT / 'engine_rounds_seed1.csv', index=False)
    print(df[(df.rounds_base != df.rounds_swap) | ~df.same_draws].head(12).to_string())


def table():
    """ Per arm, seed and year: the measures, from the runner's store """
    df = swaps.load(CONFIG)
    arm = {swaps.name(c): k for k, c in HYBRIDS.items()}
    df['arm'] = df.name.map(arm)
    return df[['arm', 'seed', 'year'] + MEASURES]


def boot_vrf(a, b, n=2000, rng=np.random.default_rng(0)):
    i = rng.integers(0, len(a), (n, len(a)))
    va = a[i].var(1, ddof=1) + b[i].var(1, ddof=1)
    vp = (a[i] - b[i]).var(1, ddof=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.nanpercentile(va / vp, [2.5, 97.5])


def analyse():
    df = table()
    df.to_csv(OUT / 'yearly_by_seed.csv.gz', index=False)
    wide = df.set_index(['arm', 'seed', 'year'])
    base, swap, ctrl = (wide.loc[k] for k in ('base', 'swap', 'control'))
    seeds = sorted(base.index.get_level_values('seed').unique())
    assert len(seeds) == CONFIG['seeds'], seeds
    # The control: exactly 0
    d = (ctrl - base).loc[(slice(None), slice(1985, 2025)), :]
    nz = (d.fillna(0) != 0).sum().sum()
    pd.DataFrame(dict(nonzero=[int(nz)], cells=[d.size], max_abs=[float(np.nanmax(np.abs(d.to_numpy())))])).to_csv(OUT / 'control.csv', index=False)
    print(f'control: {nz} nonzero differences of {d.size}')
    rows = []
    for m in MEASURES:
        for y in YEARS:
            a = swap.xs(y, level='year')[m].reindex(seeds).to_numpy(float)
            b = base.xs(y, level='year')[m].reindex(seeds).to_numpy(float)
            ok = ~(np.isnan(a) | np.isnan(b))
            a, b = a[ok], b[ok]
            dp = a - b
            du = a - np.roll(b, 1) # Different seeds
            vp, vu = dp.var(ddof=1), a.var(ddof=1) + b.var(ddof=1)
            lo, hi = boot_vrf(a, b)
            md = dp.mean()
            r = dict(measure=m, year=y, n=len(a), base_mean=b.mean(), swap_mean=a.mean(), delta=md, sd_paired=np.sqrt(vp), sd_unpaired=np.sqrt(vu),
                     sd_unpaired_empirical=du.std(ddof=1), corr=np.corrcoef(a, b)[0, 1], vrf=vu/vp if vp > 0 else np.inf, vrf_lo=lo, vrf_hi=hi)
            for tag, h in (('rel10', 0.1*abs(md)), ('base5', 0.05*abs(b.mean()))):
                for kind, v in (('paired', vp), ('unpaired', vu)):
                    r[f'n_{tag}_{kind}'] = int(np.ceil(v*(1.96/h)**2)) if h > 0 else np.nan
            rows.append(r)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / 'vrf.csv', index=False, float_format='%.5g')
    print(res[['measure', 'year', 'delta', 'corr', 'vrf', 'vrf_lo', 'vrf_hi', 'n_rel10_paired', 'n_rel10_unpaired']].to_string())
    divergence()
    plots(df, res)


def divergence():
    f = OUT / 'divergence.csv.gz'
    if not f.exists():
        return
    d = pd.read_csv(f)
    g = d.groupby(['arm', 'year']).median(numeric_only=True).drop(columns='seed').reset_index()
    first = d[d.any_state > 0].groupby(['arm', 'seed']).year.min().groupby('arm').describe()
    first.to_csv(OUT / 'first_divergence.csv')
    print('first step any starting agent differs:\n', first)
    s = g[g.arm == 'swap'].set_index('year')
    tab = s.loc[[s.index[np.abs(s.index - y).argmin()] for y in [1986, 1988, 1990, 1991, 1992, 1995, 2000, 2005, 2010, 2015, 2020, 2024]]]
    tab.to_csv(OUT / 'divergence_by_year.csv', float_format='%.4g')
    print(tab.to_string())
    # Check: the divergence runs (with the extra analyzer) equal the runner's stored runs
    dy = pd.read_csv(OUT / 'divergence_yearly.csv.gz')
    st = table()
    m = dy.merge(st, on=['arm', 'seed', 'year'], suffixes=('', '_store'))
    bad = sum((~np.isclose(m[c], m[c + '_store'], equal_nan=True, rtol=0, atol=0)).sum() for c in MEASURES)
    print(f'divergence runs vs runner store: {bad} differing values of {len(m)*len(MEASURES)}')
    pd.DataFrame(dict(differing=[int(bad)], values=[len(m)*len(MEASURES)])).to_csv(OUT / 'divergence_vs_store.csv', index=False)


def plots(df, res):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    blue, orange, grey = '#1f5fa8', '#d9822b', '#888888'
    # 1. VRF per measure and year
    fig, ax = plt.subplots(figsize=(9, 6))
    ms = MEASURES[::-1]
    for j, y in enumerate(YEARS):
        r = res[res.year == y].set_index('measure').loc[ms]
        yy = np.arange(len(ms)) + (j - len(YEARS)/2)*0.1
        ax.errorbar(r.vrf, yy, xerr=[r.vrf - r.vrf_lo, r.vrf_hi - r.vrf], fmt='o', ms=3, lw=0.8, color=plt.cm.viridis(j/(len(YEARS)-1)), label=str(y))
    ax.axvline(1, color=grey, lw=0.8)
    ax.axvline(2, color=grey, lw=0.8, ls='--')
    ax.set_xscale('log')
    ax.set_yticks(range(len(ms)), ms)
    ax.set_xlabel('Variance reduction factor, var(unpaired Δ) / var(paired Δ) (95% bootstrap)')
    ax.set_title('Pairing seeds: C swap (HIVsim → condensed cascade), 60 seeds')
    ax.legend(title='Year', fontsize=7, loc='lower right')
    fig.tight_layout()
    fig.savefig(FIG / 'vrf.png', dpi=150)
    # 2. Paired trajectories, a few seeds
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.5))
    w = df.set_index(['arm', 'seed', 'year'])
    for ax, m in zip(axs.flat, ['plhiv', 'art_coverage', 'new_infections', 'hiv_deaths', 'tx_undiagnosed', 'tx_interrupted']):
        for s, c in zip([1, 2, 3], ['#1f5fa8', '#2a9d8f', '#9b5de5']):
            b = w.loc[('base', s)][m]
            ax.plot(b.index, b, color=c, lw=1, label=f'seed {s}, HIVsim')
            ax.plot(b.index, w.loc[('swap', s)][m], color=c, lw=1, ls='--', label=f'seed {s}, C swapped')
        ax.set_title(m, fontsize=9)
        ax.set_xlim(1985, 2025)
    axs[0, 0].legend(fontsize=6)
    fig.suptitle('Paired runs (same seed): HIVsim (solid) and HIVsim with the condensed cascade (dashed)')
    fig.tight_layout()
    fig.savefig(FIG / 'paired_trajectories.png', dpi=150)
    # 3. Divergence of starting agents
    f = OUT / 'divergence.csv.gz'
    if f.exists():
        d = pd.read_csv(f)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for col, c in [('any_state', 'k'), ('alive', grey), ('infected', blue), ('diagnosed', '#2a9d8f'), ('on_art', orange), ('partners', '#9b5de5')]:
            q = d[d.arm == 'swap'].groupby('year')[col].quantile([0.1, 0.5, 0.9]).unstack()
            ax.plot(q.index, q[0.5], color=c, lw=1.2, label=col)
            ax.fill_between(q.index, q[0.1], q[0.9], color=c, alpha=0.15, lw=0)
        ctrl = d[d.arm == 'control'].any_state.max()
        ax.axvline(1990, color=grey, lw=0.6, ls=':')
        ax.set_yscale('symlog', linthresh=1e-4)
        ax.set_ylim(0, 0.3)
        ax.set_ylabel('Share of starting agents that differ')
        ax.set_title(f'Paired runs, base vs C swapped (median, 10–90%, 20 seeds); control max = {ctrl:g}', fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG / 'divergence.png', dpi=150)


if __name__ == '__main__':
    {'runs': runs, 'diverge': diverge, 'rounds': rounds, 'analyse': analyse}[sys.argv[1]]()
