"""
Tables for analysis 24 (run.py tables): outputs/emod_cells.csv, optima_ode.csv, optima_abm.csv, deltas.csv, gap.csv,
open_item.csv and report.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))
import run # noqa: E402

OUT = HERE / 'outputs'
YEARS = [2020, 2024]
SHARES = run.FIG3A
Q = ['new_infections', 'sexual', 'mtct'] + SHARES + ['plhiv', 'art_coverage']


def ci(x):
    """ Mean and 95% t interval """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    m, se = x.mean(), x.std(ddof=1)/np.sqrt(len(x))
    h = stats.t.ppf(0.975, len(x) - 1)*se
    return m, m - h, m + h, len(x)


def ci2(a, b):
    """ Unpaired difference a − b: mean and 95% interval (Welch) """
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b)
    se = np.sqrt(va + vb)
    df = (va + vb)**2/(va**2/(len(a) - 1) + vb**2/(len(b) - 1))
    d = a.mean() - b.mean()
    h = stats.t.ppf(0.975, df)*se
    return d, d - h, d + h, len(a) + len(b)


def from_yearly(y, scale):
    """ The quantities per year from common.outputs.yearly(): counts scaled to people; shares in % """
    y = y.set_index('year')
    out = pd.DataFrame(index=y.index)
    out['new_infections'] = y.new_infections*scale
    out['sexual'] = y.tx_n*scale
    out['mtct'] = out.new_infections - out.sexual # New infections that aren't sexual transmissions (mother to child; and seeds)
    out['undiagnosed'] = 100*y.tx_undiagnosed
    out['diagnosed_naive'] = 100*y.tx_diagnosed_naive
    out['on_art'] = 100*(y.tx_art_lt6m + y.tx_art_supp + y.tx_art_nonsupp)
    out['interrupted'] = 100*y.tx_interrupted
    out['plhiv'] = y.plhiv*scale
    out['art_coverage'] = 100*y.art_coverage
    return out


# %% Direction 1: EMOD′-fixed with its own cascade, and with Optima's

def emod_cells():
    from common import swaps
    df = swaps.load(run.emod_config())
    rows = []
    for (name, cell), d in df.groupby(['name', 'cell']):
        q = from_yearly(d.drop(columns=['name', 'cell']), 100) # EMOD′: 100 people per agent (components.SCALE)
        for yr in YEARS:
            rows.append(dict(arm='optima_C' if 'C=optima' in name else 'own', cell=cell, pair=d.pair.iloc[0], seed=d.seed.iloc[0], year=yr, **q.loc[yr]))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / 'emod_cells.csv', index=False)
    return out


# %% Direction 2: Optima′ with its own cascade, and with EMOD′-fixed's

def ode_quantities(o):
    yrs = list(o.years)
    st = o.by_state
    sh = 100*st.div(st.sum(axis=1), axis=0)
    rows = {}
    for yr in YEARS:
        i = yrs.index(yr)
        plhiv = np.asarray(o.plhiv).sum(axis=0)[i]
        rows[yr] = dict(new_infections=o.inci_year[i], sexual=o.sexinci_year[i], mtct=o.inci_year[i] - o.sexinci_year[i],
                        undiagnosed=sh.undx[yr], diagnosed_naive=sh.dx[yr] + sh.care[yr], on_art=sh.usvl[yr] + sh.svl[yr], interrupted=sh.lost[yr],
                        plhiv=plhiv, art_coverage=100*np.asarray(o.treat).sum(axis=0)[i]/plhiv)
    return pd.DataFrame(rows).T


def optima_ode():
    rows = []
    for f in sorted(OUT.glob('ode_*.obj')):
        q = ode_quantities(sc.load(f))
        for yr in YEARS:
            rows.append(dict(variant=f.stem.removeprefix('ode_'), year=yr, **q.loc[yr]))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / 'optima_ode.csv', index=False)
    return out


def optima_abm():
    rows = []
    for f in sorted((OUT / 'abm').glob('*.obj')):
        r = sc.load(f)
        q = from_yearly(r.yearly, r.scale)
        for yr in YEARS:
            rows.append(dict(variant=r.variant, seed=r.seed, year=yr, **q.loc[yr]))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / 'optima_abm.csv', index=False)
    return out


# %% Δs, the gap, symmetry

def fmt(m, lo, hi, q):
    if q in ('new_infections', 'sexual', 'mtct', 'plhiv'):
        return f'{m/1e3:+.1f}k ({lo/1e3:+.1f}, {hi/1e3:+.1f})'
    return f'{m:+.1f} ({lo:+.1f}, {hi:+.1f})'


def level(m, q):
    return f'{m/1e3:.1f}k' if q in ('new_infections', 'sexual', 'mtct', 'plhiv') else f'{m:.1f}'


def main():
    em, od, ab = emod_cells(), optima_ode(), optima_abm()
    lines, drows = [], []
    # Direction 1: paired by cell
    for yr in YEARS:
        a = em[(em.year == yr) & (em.arm == 'own')].set_index('cell')
        b = em[(em.year == yr) & (em.arm == 'optima_C')].set_index('cell').reindex(a.index)
        for q in Q:
            m, lo, hi, n = ci(b[q] - a[q])
            drows.append(dict(direction='1: Optima C into EMOD′-fixed', basis=f'agents, paired cells (sets {run.PIDS} x {run.EMOD_SEEDS} seeds)', year=yr, quantity=q,
                              own=a[q].mean(), swapped=b[q].mean(), delta=m, lo=lo, hi=hi, n=n))
            for pid, (ai, bi) in enumerate(zip([a[a.pair == p] for p in range(3)], [b[a.pair == p] for p in range(3)])):
                m, lo, hi, n = ci(bi[q] - ai[q])
                drows.append(dict(direction=f'1, set {run.PIDS[pid]}', basis='agents, paired', year=yr, quantity=q, own=ai[q].mean(), swapped=bi[q].mean(), delta=m, lo=lo, hi=hi, n=n))
    # Direction 2: ODE (exact) per table; agents (unpaired: own seeds 1-64 against each table's 16)
    for yr in YEARS:
        own = od[(od.variant == 'own') & (od.year == yr)].iloc[0]
        for v in sorted(set(od.variant) - {'own'}):
            s = od[(od.variant == v) & (od.year == yr)].iloc[0]
            for q in Q:
                drows.append(dict(direction=f'2: EMOD C into Optima′ ({v})', basis='ODE', year=yr, quantity=q, own=own[q], swapped=s[q], delta=s[q] - own[q], lo=np.nan, hi=np.nan, n=1))
        a = ab[(ab.variant == 'own') & (ab.year == yr)]
        for v in ['emodfix0', 'emodfix25', 'emodfix50', 'emodfix']:
            b = ab[(ab.variant.str.startswith(v)) & (ab.year == yr)] if v == 'emodfix' else ab[(ab.variant == v) & (ab.year == yr)]
            if not len(b):
                continue
            for q in Q:
                m, lo, hi, n = ci2(b[q], a[q])
                drows.append(dict(direction=f'2: EMOD C into Optima′ ({v if v != "emodfix" else "emodfix, sets 0/25/50"})', basis=f'agents ({len(b)} vs own {len(a)}, unpaired)',
                                  year=yr, quantity=q, own=a[q].mean(), swapped=b[q].mean(), delta=m, lo=lo, hi=hi, n=n))
    D = pd.DataFrame(drows)
    D.to_csv(OUT / 'deltas.csv', index=False)

    # The gap between our ends, and the fraction of it each direction's swap carries
    grows = []
    for yr in YEARS:
        e_own = D[(D.direction == '1: Optima C into EMOD′-fixed') & (D.year == yr)].set_index('quantity')
        o_ode = D[(D.direction == '2: EMOD C into Optima′ (emodfix0)') & (D.basis == 'ODE') & (D.year == yr)].set_index('quantity')
        ode_all = {v: D[(D.direction == f'2: EMOD C into Optima′ (emodfix{p})') & (D.basis == 'ODE') & (D.year == yr)].set_index('quantity') for v, p in zip(run.PIDS, run.PIDS)}
        o_abm = D[(D.direction == '2: EMOD C into Optima′ (emodfix, sets 0/25/50)') & (D.year == yr)].set_index('quantity')
        for q in Q:
            E, O = e_own.loc[q, 'own'], o_ode.loc[q, 'own']
            gap = O - E # Optima′ − EMOD′-fixed
            d1 = e_own.loc[q]
            d2_ode = np.mean([ode_all[p].loc[q, 'delta'] for p in run.PIDS])
            row = dict(year=yr, quantity=q, emod_fixed=E, optima_ode=O, optima_abm=o_abm.loc[q, 'own'] if len(o_abm) else np.nan, gap=gap,
                       d1=d1.delta, d1_lo=d1.lo, d1_hi=d1.hi, frac1=d1.delta/gap, frac1_lo=d1.lo/gap, frac1_hi=d1.hi/gap,
                       d2_ode=d2_ode, frac2_ode=d2_ode/-gap)
            if len(o_abm):
                d2 = o_abm.loc[q]
                gap_abm = d2.own - E
                row.update(d2_abm=d2.delta, d2_lo=d2.lo, d2_hi=d2.hi, frac2_abm=d2.delta/-gap_abm, frac2_abm_lo=d2.lo/-gap_abm, frac2_abm_hi=d2.hi/-gap_abm)
            grows.append(row)
    G = pd.DataFrame(grows)
    G.to_csv(OUT / 'gap.csv', index=False)

    for yr in YEARS:
        lines.append(f'===== {yr} =====')
        lines.append('Direction 1, EMOD′-fixed: own cascade → Optima\'s table (paired, 48 cells): own; with Optima\'s C; Δ (95% CI)')
        for q in Q:
            r = D[(D.direction == '1: Optima C into EMOD′-fixed') & (D.year == yr) & (D.quantity == q)].iloc[0]
            lines.append(f'  {q:16s} {level(r.own, q):>8s} {level(r.swapped, q):>8s}  {fmt(r.delta, r.lo, r.hi, q)}')
        lines.append('Direction 2, Optima′: own cascade → EMOD′-fixed\'s table. ODE per set (exact): own; with EMOD\'s C; Δ.  Agents, sets pooled (unpaired): Δ (95% CI)')
        for q in Q:
            o = [D[(D.direction == f'2: EMOD C into Optima′ (emodfix{p})') & (D.basis == 'ODE') & (D.year == yr) & (D.quantity == q)].iloc[0] for p in run.PIDS]
            s = f'  {q:16s} {level(o[0].own, q):>8s} ' + ' / '.join(level(x.swapped, q) for x in o) + '  Δ ' + ' / '.join(level(x.delta, q) if q in SHARES + ['art_coverage'] else f'{x.delta/1e3:+.1f}k' for x in o)
            r = D[(D.direction == '2: EMOD C into Optima′ (emodfix, sets 0/25/50)') & (D.year == yr) & (D.quantity == q)]
            if len(r):
                r = r.iloc[0]
                s += f';  agents own {level(r.own, q)}, swapped {level(r.swapped, q)}, Δ {fmt(r.delta, r.lo, r.hi, q)}'
            lines.append(s)
        for v in ['emodbug0', 'emodfix0_noflows']:
            o = D[(D.direction == f'2: EMOD C into Optima′ ({v})') & (D.year == yr)].set_index('quantity')
            if len(o):
                lines.append(f'  ODE {v}: ' + ', '.join(f'{q} {level(o.loc[q, "swapped"], q)}' for q in Q))
        lines.append('Gap (Optima′ ODE − EMOD′-fixed) and the fraction each swap carries: direction 1 (95% CI); direction 2 ODE (mean of sets); direction 2 agents (95% CI)')
        for _, r in G[G.year == yr].iterrows():
            s = f'  {r.quantity:16s} gap {level(r.gap, r.quantity):>8s}: 1 → {r.frac1:+.2f} ({r.frac1_lo:+.2f}, {r.frac1_hi:+.2f}); 2 ODE → {r.frac2_ode:+.2f}'
            if 'frac2_abm' in r and np.isfinite(r.frac2_abm):
                s += f'; 2 agents → {r.frac2_abm:+.2f} ({r.frac2_abm_lo:+.2f}, {r.frac2_abm_hi:+.2f})'
            lines.append(s)
        lines.append('')
    lines += open_item(ab)
    (OUT / 'report.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def open_item(ab):
    """ Analysis 23's open item: Optima′ agents, Optima's table against its own cascade, same seeds and code (64 each) """
    a20 = run.optima_imports()
    yrs = np.arange(2000, 2021)
    files = lambda v: sorted((OUT / 'abm').glob(f'{v}_*.obj'), key=lambda f: int(f.stem.split('_')[-1]))
    own, tab = [sc.load(f) for f in files('own')], [sc.load(f) for f in files('optima')]
    if not own or not tab:
        return []
    seeds = sorted(set(x.seed for x in own) & set(x.seed for x in tab))
    own = [x for x in own if x.seed in seeds]
    tab = [x for x in tab if x.seed in seeds]
    ode = sc.load(OUT / 'ode_own.obj')
    Y = list(own[0].extract.years)
    idx = [Y.index(y) for y in yrs]
    tot = lambda o, k: np.asarray(o[k]).sum(axis=0)[idx] if np.ndim(o[k]) == 2 else np.asarray(o[k])[idx]
    lines = [f'===== Analysis 23\'s open item: Optima′ agents ({run.N_AGENTS:,}), Optima\'s table vs own cascade, {len(seeds)} seeds each (same seeds, same code) =====',
             'mean over 2000-2020 of (table/own − 1), with its 95% CI from the seeds (paired by seed, and unpaired); and each against the ODE']
    rows = []
    for k, lab in dict(plhiv='PLHIV', diag='Diagnosed', treat='On ART', supp='Suppressed', inci_year='New infections (year)', deaths_year='HIV deaths (year)').items():
        A, B = np.array([tot(x.extract, k) for x in tab]), np.array([tot(x.extract, k) for x in own])
        o = tot(ode, k)
        rel_seed = (A/B - 1).mean(axis=1) # Per seed: mean over years
        mp, lp, hp, _ = ci(100*rel_seed)
        ra, rb = 100*(A/o - 1).mean(axis=1), 100*(B/o - 1).mean(axis=1)
        mu, lu, hu, _ = ci2(ra, rb)
        z = (A.mean(0) - B.mean(0))/np.sqrt(A.var(0, ddof=1)/len(A) + B.var(0, ddof=1)/len(B))
        rows.append(dict(output=lab, paired=mp, paired_lo=lp, paired_hi=hp, unpaired=mu, unpaired_lo=lu, unpaired_hi=hu, table_vs_ode=ra.mean(), own_vs_ode=rb.mean(),
                         years_absz_gt2=int((np.abs(z) > 2).sum()), r_seed=np.corrcoef(A.mean(1), B.mean(1))[0, 1]))
        lines.append(f'  {lab:22s} paired {mp:+.2f}% ({lp:+.2f}, {hp:+.2f}); unpaired {mu:+.2f}% ({lu:+.2f}, {hu:+.2f}); vs ODE: table {ra.mean():+.2f}%, own {rb.mean():+.2f}%; |z| > 2 in {int((np.abs(z) > 2).sum())}/21 years')
    pd.DataFrame(rows).to_csv(OUT / 'open_item.csv', index=False)
    for yr in YEARS:
        a, b = ab[(ab.variant == 'optima') & (ab.year == yr)], ab[(ab.variant == 'own') & (ab.year == yr)]
        lines.append(f'  Figure 3A {yr}, table vs own: ' + '; '.join(f'{q} {a[q].mean():.1f} vs {b[q].mean():.1f} ({fmt(*ci2(a[q], b[q])[:3], q)})' for q in SHARES) +
                     f'; new infections {fmt(*ci2(a.new_infections, b.new_infections)[:3], "new_infections")}')
    return lines


if __name__ == '__main__':
    main()
