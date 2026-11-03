"""How much do Optima's Malawi outputs depend on its time step? Bounds how close a continuous-time ODE can get."""
import numpy as np, sciris as sc
import optima as op

P = op.loadproj('Malawi_20220912.prj', verbose=0)
keys = ['numplhiv', 'numinci', 'numdeath', 'numdiag', 'numtreat', 'popsize']
out = {}
for dt in [0.2, 0.1, 0.05, 0.02]:
    T = sc.timer()
    res = P.runsim(name='20220912', end=2030, dt=dt, verbose=0)
    out[dt] = {k: res.main[k].tot[0] for k in keys}
    print(f'dt={dt}: {T.toc(output=True):.1f} s')
years = res.tvec
sel = (years >= 2000) & (years <= 2030)
print(f'\nMax relative difference from dt=0.2, years 2000-2030 (%):')
print('dt    ' + ''.join(f'{k:>10s}' for k in keys))
for dt in [0.1, 0.05, 0.02]:
    print(f'{dt:<6}' + ''.join(f'{100*np.max(np.abs(out[dt][k][sel]/out[0.2][k][sel]-1)):10.1f}' for k in keys))
