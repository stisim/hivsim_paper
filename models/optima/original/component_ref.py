"""Reference runs for checking Optima′ one component at a time: Optima with the components not yet built switched off,
from the step-independent harness (model_dt, FIX={'linear','inputs'}) at a small step, and from the original at dt=0.2"""
import sys, numpy as np, sciris as sc
import optima as op
import model_dt, dt_fix

END = 2030

def switch_off(sp, what):
    """Switch off components in a simpars dict, from the second time point on (so the initial distribution is Optima's own)"""
    sp = sc.dcp(sp)
    if 'transmission' in what: sp['force'] = 0*sp['force']
    if 'mtct' in what:
        sp['mtctbreast'] = sp['mtctnobreast'] = 0.
        sp['numpmtct'] = 0*sp['numpmtct']
    if 'cascade' in what:
        for key in ['hivtest', 'leavecare', 'returntocare']: sp[key][..., 1:] = 0.
        for key in ['aidstest', 'aidsleavecare']: sp[key][1:] = 0.
        sp['linktocare'][:, 1:] = 1e9; sp['aidslinktocare'][1:] = 1e9
    if 'treatment' in what:
        sp['numtx'] = 0*sp['numtx']
    return sp

TESTS = sc.objdict(
    A = ['transmission', 'mtct', 'cascade', 'treatment'], # Demography and natural history only
    B = ['mtct', 'cascade', 'treatment'], # Plus transmission
    C = ['mtct', 'treatment'], # Plus the cascade (testing, linkage, loss to follow-up, return to care)
    D = ['mtct'], # Plus treatment (ART targets, suppression, held proportions)
    E = [], # The whole model
)

def run(simpars, fix):
    model_dt.FIX = set(fix)
    r = dt_fix.P.runsim(name='20220912', simpars=simpars, verbose=0, keepraw=True)
    raw = r.raw[0] if isinstance(r.raw, list) else list(r.raw.values())[0]
    t = raw['tvec']
    yearly = [np.argmin(abs(t - y)) for y in range(1990, END+1)]
    return dict(people=raw['people'][:, :, yearly], death=raw['death'][:, :, yearly], otherdeath=raw['otherdeath'][:, yearly],
                births=raw['births'][:, yearly], inci=raw['inci'][:, yearly], years=np.arange(1990, END+1))

if __name__ == '__main__':
    dt = 0.02
    dt_fix.END = END
    out = sc.objdict()
    for name in sys.argv[1:] or TESTS.keys():
        what = TESTS[name]
        kw = dict(projectversion=dt_fix.P.version, start=1990, end=END, settings=dt_fix.P.settings, name='20220912')
        orig = switch_off(op.makesimpars(dt_fix.P.parsets['20220912'].pars, dt=0.2, **kw), what)
        fine = switch_off(dt_fix.ref_inputs(dt, END), what)
        out[name] = dict(what=what, dt=dt, original=run(orig, []), harness=run(fine, ['linear']))
        o, h = out[name]['original'], out[name]['harness']
        plhiv = lambda x: x['people'][2:].sum(axis=(0, 1))
        print(f'{name}: PLHIV 2000/2020 original {plhiv(o)[10]:.0f}/{plhiv(o)[30]:.0f}, harness {plhiv(h)[10]:.0f}/{plhiv(h)[30]:.0f}')
    try:    old = sc.load('component_refs.obj') # Keep the tests not rerun
    except FileNotFoundError: old = {}
    sc.save('component_refs.obj', sc.mergedicts(old, out))
