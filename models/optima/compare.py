"""Compare Optima′ with Optima, one component test at a time (references from original/component_ref.py)"""
import sys
from pathlib import Path
import numpy as np, sciris as sc
import optima_prime as opp

REFS = Path(__file__).resolve().parent / 'original' / 'component_refs.obj'
COMPONENTS = dict(A=('demography', 'natural_history'), B=('demography', 'natural_history', 'transmission'), C=('demography', 'natural_history', 'transmission', 'cascade'), D=('demography', 'natural_history', 'transmission', 'cascade', 'treatment'), E=('demography', 'natural_history', 'transmission', 'cascade', 'treatment', 'mtct'))
WINDOWS = dict(A=(1991, 2000), B=(1991, 2020), C=(1991, 2020), D=(1991, 2020), E=(1991, 2020)) # Without transmission (A), people with HIV die out: compare while there are some

def summarise(people, years, inci=None, deaths=None):
    """ Stocks from a (states × pops × years) array, and new infections per year if given """
    hiv = people[2:]
    out = dict(popsize=people.sum(axis=(0, 1)), plhiv=hiv.sum(axis=(0, 1)))
    for i, k in enumerate(opp.CD4):
        out[f'plhiv_{k}'] = hiv[i::6].sum(axis=(0, 1))
    out['plhiv_15_49'] = None
    out['diag'] = people[8:].sum(axis=(0, 1)) # Diagnosed: in dx, care, lost, or on ART
    out['treat'] = people[26:].sum(axis=(0, 1)) # On ART
    out['supp'] = people[32:].sum(axis=(0, 1)) # Virally suppressed
    out['inci'] = inci
    out['deaths'] = deaths
    return {k: v for k, v in out.items() if v is not None}

def prime_people(sim, years):
    """ Optima′'s compartments as a (states × pops × years) array, at whole years """
    r = sim.results.optima
    names = [s.replace('-', '_') for s in opp.load().settings['statelabels']]
    out = []
    for n in names:
        res = r[n].by('pop') # One series per population (on its own time points)
        t = np.array([d.years for d in res.timevec])
        out.append(np.asarray(res)[[np.argmin(abs(t - y)) for y in years]].T)
    return np.array(out)

def prime_inci(sim, years):
    """ Optima′'s new infections per year, in the step that starts at each whole year """
    r = sim.results.optima
    flows = [r[k] for k in r.keys() if k.startswith('infection_') or k in ('birth_undx', 'birth_dx')] # Optima counts mother-to-child transmission as new infections
    if not flows: return None
    t = np.array([d.years for d in flows[0].timevec])
    total = sum(np.asarray(f) for f in flows)
    return np.array([total[min(np.argmin(abs(t - y)) + 1, len(t) - 1)] for y in years])/sim.pars.dt.years if hasattr(sim.pars.dt, 'years') else None

def prime_deaths(sim, years):
    """ Optima′'s HIV deaths per year, in the step that starts at each whole year, as Optima reports them (deaths on ART scaled by rrcomorbiditydeathtx) """
    r = sim.results.optima
    rr = opp.load().simpars['rrcomorbiditydeathtx']
    tvec = opp.load().simpars['tvec']
    total = 0
    for k in [k for k in r.keys() if k.startswith('hivdeath_')]:
        res = r[k].by('pop')
        t = np.array([d.years for d in res.timevec])
        inds = [min(np.argmin(abs(t - y)) + 1, len(t) - 1) for y in years]
        flows = np.asarray(res)[inds] # Years × pops
        if k.split('_')[1] in opp.TREATED:
            flows = flows*np.array([[np.interp(y, tvec, row) for row in rr] for y in years])
        total = total + flows.sum(axis=1)
    return total/(t[1] - t[0])

def gaps(a, b, years, window):
    """ Max difference over the window, as a % of the largest value of the reference in the window """
    sel = (years >= window[0]) & (years <= window[1])
    return {k: 100*np.max(np.abs(b[k][sel] - a[k][sel]))/np.max(np.abs(a[k][sel])) for k in a if b.get(k) is not None and np.max(np.abs(a[k][sel])) > 0}

if __name__ == '__main__':
    refs = sc.load(REFS)
    for name in sys.argv[1:] or COMPONENTS:
        ref = refs[name]
        years = ref['original']['years']
        T = sc.timer()
        sim = opp.make_sim(components=COMPONENTS[name], dt=ref['dt'], end=int(years[-1])).run()
        print(f'Test {name} ({", ".join(COMPONENTS[name])}): Optima′ ran in {T.toc(output=True):.1f} s')
        inci = lambda x: x['inci'].sum(axis=0) if 'inci' in x and 'transmission' in COMPONENTS[name] else None
        deaths = lambda x: x['death'].sum(axis=(0, 1))
        o, h = summarise(ref['original']['people'], years, inci(ref['original']), deaths(ref['original'])), summarise(ref['harness']['people'], years, inci(ref['harness']), deaths(ref['harness']))
        p = summarise(prime_people(sim, years), years, prime_inci(sim, years) if 'transmission' in COMPONENTS[name] else None, prime_deaths(sim, years))
        window = WINDOWS[name]
        print(f'Max difference, {window[0]}-{window[1]}, as % of the reference:')
        print(f'{"":28s}' + ''.join(f'{k[6:] or k:>9s}' for k in o))
        for label, a, b in [('harness vs original', o, h), ("Optima′ vs harness", h, p), ("Optima′ vs original", o, p)]:
            g = gaps(a, b, years, window)
            print(f'{label:28s}' + ''.join(f'{g.get(k, np.nan):9.2f}' for k in o))
