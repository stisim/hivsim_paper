"""Find Optima's step-dependent rules: convert groups of them to a step-independent form and measure the gap to the calibrated run (dt=0.2).

Compared over the years with treatment data (2000-2020), since Optima's projections beyond the data amplify small differences.
"""
import sys, numpy as np
import optima as op
import model_dt

op.project.model = model_dt.model
P = op.loadproj('Malawi_20220912.prj', verbose=0)
keys = ['numplhiv', 'numinci', 'numdeath', 'numdiag', 'numtreat', 'numsuppressed']
START, END = 2000, 2020

TIME  = ['supp', 'acute', 'dwell']                     # Keep the time spent in a state (exponential stays); r*dt rules left as rate r
SURV  = ['supp', 'cascade', 'death', 'usvl', 'age']    # Keep the survival over the step (r*dt rules); exponential stays left as 1/tau
RULESETS = {
    'none':             [],
    'supp only':        ['supp'],
    'time in state':    TIME,
    'survival':         SURV,
    'mixed':            TIME + SURV[1:],
    'survival + trans': SURV + ['trans'],
    'linear':           ['linear'],
    'linear + inputs':  ['linear', 'inputs'],
}

def ref_inputs(dt, end=None):
    """Inputs as Optima makes them at REF_DT, carried to a finer grid: linear between REF_DT points, or the previous point where a series has NaN or inf markers"""
    kw = dict(projectversion=P.version, start=P.settings.start, end=end or END, settings=P.settings, name='20220912')
    pars = P.parsets['20220912'].pars
    ref, fine = op.makesimpars(pars, dt=model_dt.REF_DT, **kw), op.makesimpars(pars, dt=dt, **kw)
    tr, tf = ref['tvec'], fine['tvec']
    prev = np.clip(np.searchsorted(tr, tf + 1e-9) - 1, 0, len(tr)-1)
    for key, val in ref.items():
        if key == 'tvec' or not isinstance(val, np.ndarray) or val.dtype.kind != 'f' or val.ndim == 0 or val.shape[-1] != len(tr):
            continue
        if np.all(np.isfinite(val)):
            fine[key] = np.apply_along_axis(lambda y: np.interp(tf, tr, y), -1, val) if val.ndim > 1 else np.interp(tf, tr, val)
        else:
            fine[key] = val[..., prev]
    return fine

def run(dt, fix=()):
    model_dt.FIX = set(fix)
    if 'inputs' in fix and dt != model_dt.REF_DT:
        r = P.runsim(name='20220912', simpars=ref_inputs(dt), verbose=0)
    else:
        r = P.runsim(name='20220912', end=END, dt=dt, verbose=0)
    return r, {k: r.main[k].tot[0] for k in keys}

ref_res, ref = run(0.2)
years = ref_res.tvec
sel = (years >= START) & (years <= END)
def gap(out):
    return {k: 100*np.max(np.abs(out[k][sel]/ref[k][sel] - 1)) for k in keys}
def show(label, g):
    print(f'{label:32s}' + ''.join(f'{g[k]:9.2f}' for k in keys) + f'{max(g.values()):9.2f}')

if __name__ == '__main__':
    print(f'Max % difference from the calibrated run (dt=0.2), {START}-{END}')
    print(f'{"":32s}' + ''.join(f'{k[3:]:>9s}' for k in keys) + f'{"max":>9s}')
    for dt in [float(a) for a in sys.argv[1:]] or [0.02]:
        print(f'-- dt={dt}')
        for name, fix in RULESETS.items():
            show(name, gap(run(dt, fix)[1]))
            show(name + ' + new ageing', gap(run(dt, fix + ['agebug'])[1]))
