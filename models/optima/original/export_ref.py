"""Export the Optima Malawi inputs (simpars at dt=0.2) and the calibrated run's compartment-level outputs, for building and checking Optima′."""
import numpy as np, sciris as sc
import optima as op

P = op.loadproj('Malawi_20220912.prj', verbose=0)
END = 2030  # Run past the comparison window (2000-2020) so the window is not at the end of the run
kw = dict(projectversion=P.version, start=P.settings.start, end=END, dt=0.2, settings=P.settings, name='20220912')
simpars = op.makesimpars(P.parsets['20220912'].pars, **kw)
res = P.runsim(name='20220912', simpars=simpars, verbose=0, keepraw=True)
raw = res.raw[0] if isinstance(res.raw, list) else list(res.raw.values())[0]

s = P.settings
settings = dict(statelabels=s.statelabels, hivstates=s.hivstates, healthstates=s.healthstates,
                **{k: np.array(getattr(s, k)) for k in ['susreg', 'progcirc', 'sus', 'undx', 'dx', 'care', 'usvl', 'svl', 'lost', 'acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50', 'allplhiv', 'alltx', 'notonart']})
def plain(x):
    """Convert Optima's odicts to plain dicts, so the file loads without Optima installed"""
    if isinstance(x, dict): return {k: plain(v) for k, v in x.items()}
    if isinstance(x, list): return [plain(v) for v in x]
    return x

simpars = plain(dict(simpars))
out = dict(
    version=P.version, optima_version=op.version, dt=0.2,
    simpars=simpars,
    settings=settings,
    raw={k: v for k, v in raw.items() if k != 'tvec'} | {'tvec': raw['tvec']},
    totals={k: res.main[k].tot[0] for k in res.main.keys()},
    bypop={k: res.main[k].pops[0] for k in res.main.keys()},
    years=res.tvec,
)
out = plain(out)
sc.save('optima_malawi_ref.obj', out)
print('people', raw['people'].shape, 'saved', sc.prettyfilesize('optima_malawi_ref.obj') if hasattr(sc, 'prettyfilesize') else '')
