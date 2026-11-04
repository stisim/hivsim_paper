"""Run the Optima HIV Malawi project (Malawi_20220912.prj) and export reference outputs."""
import sys, time
import numpy as np, pandas as pd, sciris as sc
import optima as op

P = op.loadproj('Malawi_20220912.prj', verbose=0)
T = sc.timer()
res = P.runsim(name='20220912', end=2030, verbose=0)
print(f'runsim: {T.toc(output=True):.1f} s')

years = res.tvec
keys = ['prev', 'numplhiv', 'numinci', 'numdeath', 'numdiag', 'numtreat', 'numsuppressed', 'popsize']
rows = {}
for k in keys:
    r = res.main[k]
    rows[k] = r.tot[0]
    print(f'{k:14s} {r.name:40s} 2000={r.tot[0][np.argmin(abs(years-2000))]:.4g}  2020={r.tot[0][np.argmin(abs(years-2020))]:.4g}')
df = pd.DataFrame(rows, index=pd.Index(years, name='year'))
df.to_csv('optima_malawi_totals.csv')
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 4, figsize=(16, 7))
for ax, k in zip(axs.flat, keys):
    r = res.main[k]
    ax.plot(years, r.tot[0], label='Optima')
    dt = getattr(r, 'datatot', None)
    if dt is not None and res.datayears is not None:
        y = np.array(dt, dtype=float)
        while y.ndim > 1: y = y[0]  # best estimate (rows are best/low/high)
        ok = ~np.isnan(y)
        ax.scatter(np.array(res.datayears)[ok], y[ok], s=12, c='k', label='Data')
    ax.set_title(r.name); ax.set_xlim(1990, 2030); ax.set_ylim(bottom=0)
axs.flat[0].legend()
fig.suptitle('Optima HIV, Malawi_20220912.prj, parset 20220912')
fig.tight_layout(); fig.savefig('../../../analyses/01_reference_fits/figures/optima_totals.png', dpi=100)
