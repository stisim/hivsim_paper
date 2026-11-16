"""
Tables for analysis 26 (run.py tables): modules against declared, as distributions over seeds, by year; and both
against the ODE. outputs/compare_N.txt and outputs/by_year_N.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run # noqa: E402

OUT = run.OUT
Q = run.COUNTS + run.FIG3A[:1] + ['diagnosed_naive', 'on_art_tx', 'interrupted']
SHARES = ['undiagnosed', 'diagnosed_naive', 'on_art_tx', 'interrupted']
FLOWS = ['new_infections', 'hiv_deaths', 'sexual']


def load(arm, n):
    files = sorted(run.RUNS.glob(f'{arm}_{n}_*.obj'), key=lambda f: int(f.stem.split('_')[-1]))
    return [sc.load(f) for f in files]


def welch(a, b):
    """ a − b, its SE, z, two-sided p (Welch), and the ratio of SDs with the F test's p """
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se = np.sqrt(va/len(a) + vb/len(b))
    d = a.mean() - b.mean()
    t = d/se if se > 0 else 0.0
    df = (va/len(a) + vb/len(b))**2/((va/len(a))**2/(len(a) - 1) + (vb/len(b))**2/(len(b) - 1)) if se > 0 else 1
    p = 2*stats.t.sf(abs(t), df)
    f = va/vb if vb > 0 else np.nan
    pf = 2*min(stats.f.sf(f, len(a) - 1, len(b) - 1), stats.f.cdf(f, len(a) - 1, len(b) - 1)) if np.isfinite(f) else np.nan
    return d, se, t, p, np.sqrt(f), pf


def ode_table():
    f = OUT / 'ode.obj'
    return sc.load(f) if f.exists() else None


def main(n=run.N_AGENTS, start=2000, end=2020):
    dec, mod = load('declared', n), load('modules', n)
    k = min(len(dec), len(mod))
    years = [y for y in dec[0].yearly.index if start <= y <= end and y in mod[0].yearly.index]
    lines = [f'Optima′ as modules ({len(mod)} seeds) against Optima′ declared as agents ({len(dec)} seeds), {n:,} agents in 1990; years {years[0]}-{years[-1]}',
             f'Run time: declared {np.mean([r.time for r in dec]):.0f} s, modules {np.mean([r.time for r in mod]):.0f} s (mean per run)', '',
             'Per output: mean over years of (modules mean / declared mean − 1) in % (shares: difference in points), its largest |value|, mean z (Welch),',
             'number of years with p < 0.05 (of N; about 5% expected by chance) and the smallest p; SD ratio (modules/declared), mean over years, and years with F-test p < 0.05', '']
    lines.append(f'{"output":16s}{"mean Δ":>9s}{"max |Δ|":>9s}{"mean z":>8s}{"p<.05":>7s}{"min p":>8s}{"SD ratio":>9s}{"F p<.05":>8s}{"SE %":>7s}')
    rows = []
    for q in Q:
        A = np.array([r.yearly.loc[years, q] for r in mod])
        B = np.array([r.yearly.loc[years, q] for r in dec])
        st = np.array([welch(A[:, i], B[:, i]) for i in range(len(years))]) # Years × (d, se, z, p, sdr, pf)
        share = q in SHARES
        rel = st[:, 0] if share else 100*st[:, 0]/B.mean(0)
        se = st[:, 1] if share else 100*st[:, 1]/B.mean(0)
        for i, y in enumerate(years):
            rows.append(dict(output=q, year=y, modules=A[:, i].mean(), declared=B[:, i].mean(), diff=st[i, 0], rel=rel[i], z=st[i, 2], p=st[i, 3], sd_ratio=st[i, 4], p_sd=st[i, 5]))
        lines.append(f'{q:16s}{rel.mean():9.2f}{np.abs(rel).max():9.2f}{st[:, 2].mean():8.2f}{(st[:, 3] < 0.05).sum():4d}/{len(years):<2d}{st[:, 3].min():8.3f}{np.nanmean(st[:, 4]):9.2f}{(st[:, 5] < 0.05).sum():5d}/{len(years):<2d}{se.mean():7.2f}')
    by = pd.DataFrame(rows)
    by.to_csv(OUT / f'by_year_{n}.csv', index=False)

    # PLHIV by place in care × CD4 band, 2000-2020 (people), as above
    lines += ['', 'People with HIV by place in care × CD4 band: mean % difference over the years (modules/declared − 1), mean z, years with p < 0.05; and the 2020 share of PLHIV (declared)']
    yrs = dec[0].bp_years
    sel = (yrs >= start) & (yrs <= end)
    A = np.array([r.bandplace[sel] for r in mod]) # Seeds × years × place × band
    B = np.array([r.bandplace[sel] for r in dec])
    tot2020 = B[:, -1].sum(axis=(1, 2)).mean()
    lines.append(f'{"":6s}' + ''.join(f'{b:>22s}' for b in ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50']))
    bp_rows = []
    for i, pname in enumerate(run.PLACES):
        cells = []
        for j in range(6):
            a, b = A[:, :, i, j], B[:, :, i, j]
            st = np.array([welch(a[:, t], b[:, t]) for t in range(a.shape[1])])
            ok = b.mean(0) > 0
            rel = 100*st[ok, 0]/b.mean(0)[ok] if ok.any() else np.array([np.nan])
            share = 100*b[:, -1].mean()/tot2020
            bp_rows.append(dict(place=pname, band=j, mean_rel=np.nanmean(rel), mean_z=np.nanmean(st[:, 2]), n_p05=int((st[:, 3] < 0.05).sum()), share_2020=share))
            cells.append(f'{np.nanmean(rel):+6.1f}% z{np.nanmean(st[:, 2]):+5.1f} {int((st[:, 3] < 0.05).sum()):2d} {share:4.1f}%')
        lines.append(f'{pname:6s}' + ''.join(f'{c:>22s}' for c in cells))
    pd.DataFrame(bp_rows).to_csv(OUT / f'bandplace_{n}.csv', index=False)

    # Selected years
    lines += ['', 'Selected years: declared mean ± SE | modules mean ± SE']
    for q in Q:
        s = []
        for y in [2000, 2005, 2010, 2015, 2020]:
            if y not in years: continue
            a = np.array([r.yearly.loc[y, q] for r in mod]); b = np.array([r.yearly.loc[y, q] for r in dec])
            fmt = (lambda v: f'{v:.1f}') if q in SHARES else (lambda v: f'{v:,.0f}')
            s.append(f'{y}: {fmt(b.mean())}±{fmt(b.std(ddof=1)/np.sqrt(len(b)))} | {fmt(a.mean())}±{fmt(a.std(ddof=1)/np.sqrt(len(a)))}')
        lines.append(f'  {q:15s} ' + '  '.join(s))

    # Against the ODE
    ode = ode_table()
    if ode is not None:
        lines += ['', 'Against the ODE (optima_prime.make_sim with Optima\'s table, 2022 project; acceptance.py\'s 2% standard): mean of seeds / ODE − 1, %, by year (max |%| over the years, and the year), with the agent SE %']
        for q in [c for c in Q if c in ode.yearly.columns]:
            o = ode.yearly.loc[years, q].to_numpy(float)
            for label, runs in (('declared', dec), ('modules', mod)):
                m = np.array([r.yearly.loc[years, q] for r in runs])
                if q in SHARES:
                    d = m.mean(0) - o
                    se = m.std(0, ddof=1)/np.sqrt(len(m))
                    unit = ' pts'
                else:
                    d = 100*(m.mean(0)/o - 1)
                    se = 100*m.std(0, ddof=1)/np.sqrt(len(m))/o
                    unit = '%'
                j = np.argmax(np.abs(d))
                lines.append(f'  {q:15s} {label:9s} mean {d.mean():+6.2f}{unit}  max |{abs(d[j]):.2f}|{unit} ({years[j]})  2020 {d[-1]:+6.2f}{unit}  SE {se.mean():.2f}{unit}')
    text = '\n'.join(lines)
    print(text)
    (OUT / f'compare_{n}.txt').write_text(text + '\n')
    return by
