"""
Optima's own cascade for Malawi (the project's inputs, `simpars` at dt = 0.2, models/optima/reference) as a
condensed-cascade table of rate and level pathways (common/cascade/params.py), with each row naming what it comes from.

    python extract_cascade.py    # writes common/cascade/tables/malawi_optima_params.csv

Rates are Optima′'s: Optima's probability of a move in a step of 0.2 years, divided by 0.2 (optima_prime.cascade_probs),
so Optima′ running the table moves people as with its own cascade. Optima combines some rates for those with AIDS
(CD4 < 200, its `aidsind`): the table gives the combined rate for those bands, as its own row. Series are Optima's
values at its time points (1990-2030 by 0.2), linear between them, with points on a straight line left out; rows by
population are grouped where the populations' series are the same.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import optima_prime as opp # noqa: E402

OUT = Path(__file__).resolve().parents[3] / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv'
DT = opp.REF_DT
NOT_AIDS, AIDS = 'cd4:acute|gt500|gt350|gt200', 'cd4:gt50|lt50'


def series(tvec, v):
    """ "year: value, ..." at Optima's time points, without points on a straight line between their neighbours """
    v = np.asarray(v, dtype=float)
    keep = np.ones(len(v), bool)
    for i in range(1, len(v) - 1):
        line = v[i-1] + (v[i+1] - v[i-1])*(tvec[i] - tvec[i-1])/(tvec[i+1] - tvec[i-1])
        keep[i] = abs(line - v[i]) > 1e-13*max(1, abs(v[i]))
    return ', '.join(f'{t!r}: {x!r}' for t, x in zip(np.asarray(tvec, float)[keep].tolist(), v[keep].tolist()))


def extract(ref=None):
    ref = ref or opp.load()
    sp = ref.simpars
    tvec, pops = np.asarray(sp['tvec'], float), list(sp['popkeys'])
    P, T = len(pops), len(tvec)
    full = lambda x: np.broadcast_to(np.asarray(x, dtype=float), (P, T))
    rows = []
    def row(pathway, parameter, form, value, source, group='', places='', per=''):
        rows.append(dict(pathway=pathway, parameter=parameter, sex='all', group=group, places=places, per=per, form=form, value=value, default=None, source=source))
    def by_pop(pathway, parameter, arr, source, cd4='', **kw):
        """ One row per set of populations with the same series (no pop key if all are the same) """
        arr = full(arr)
        uniq, inv = np.unique(arr, axis=0, return_inverse=True)
        for k in range(len(uniq)):
            who = [p for p, j in zip(pops, inv.ravel()) if j == k]
            grp = '; '.join(x for x in [cd4, '' if len(who) == P else 'pop:' + '|'.join(who)] if x)
            row(pathway, parameter, 'linear', series(tvec, uniq[k]), source, group=grp, **kw)

    # Rates (optima_prime.cascade_probs, after model.py "Calculate probabilities of shifting along cascade")
    link = lambda d: (1 - np.exp(-DT/np.maximum(1e-3, full(d))))/DT
    test, aidstest, ret = full(sp['hivtest']), full(sp['aidstest']), full(sp['returntocare'])
    by_pop('rate', 'hivtest', test, 'hivtest (per year; Optima: hivtest*dt a step)', NOT_AIDS, places='undx>dx')
    by_pop('rate', 'aidstest', np.maximum(aidstest, test), 'max(aidstest, hivtest) for those with AIDS (Optima: the larger probability a step)', AIDS, places='undx>dx')
    by_pop('rate', 'linktocare', link(sp['linktocare']), 'linktocare (years to link): (1 - exp(-dt/linktocare))/dt at dt = 0.2', NOT_AIDS, places='dx>care')
    by_pop('rate', 'aidslinktocare', np.maximum(link(sp['aidslinktocare']), link(sp['linktocare'])), 'aidslinktocare: the larger of its and linktocare\'s (1 - exp(-dt/years))/dt at dt = 0.2, for those with AIDS', AIDS, places='dx>care')
    by_pop('rate', 'leavecare', full(sp['leavecare']), 'leavecare: from care, and from ART (Optima\'s loss from ART is leavecare)', NOT_AIDS, places='care|svl|usvl>lost')
    by_pop('rate', 'aidsleavecare', np.minimum(full(sp['aidsleavecare']), full(sp['leavecare'])), 'min(aidsleavecare, leavecare) for those with AIDS: from care and from ART', AIDS, places='care|svl|usvl>lost')
    by_pop('rate', 'returntocare', ret, 'returntocare', NOT_AIDS, places='lost>care')
    by_pop('rate', 'aidsreturntocare', (1 - (1 - aidstest*DT)*(1 - ret*DT))/DT, 'returntocare and aidstest for those lost with AIDS: (1 - (1 - aidstest*dt)(1 - returntocare*dt))/dt at dt = 0.2', AIDS, places='lost>care')
    held = tvec <= sp['fixpropsupp'] # Failure and regain run only while the suppressed share is not held (Optima's userate(propsupp))
    row('rate', 'treatfail', 'linear', series(tvec, np.where(held, sp['treatfail'], 0.0)), f'treatfail, until fixpropsupp ({sp["fixpropsupp"]:g})', places='svl>usvl')
    row('rate', 'regainvs', 'linear', series(tvec, np.where(held, np.asarray(sp['regainvs'])*np.asarray(sp['numvlmon']), 0.0)),
        f'regainvs × numvlmon (viral load tests a year finding failure and switching), spread over everyone on ART, until fixpropsupp ({sp["fixpropsupp"]:g})', places='usvl>svl', per='svl|usvl')

    # Levels (optima_prime.CascadeTargets, after model.py "Proportions")
    row('target', 'on_art', 'linear', series(tvec, sp['numtx']), 'numtx: from those in care; back to care when over')
    row('target', 'lowest_cd4_before', 'const', f'{sp["allcd4eligibletx"]!r}', 'allcd4eligibletx: lowest CD4 band first before it')
    row('target', 'p_suppressive', 'const', f'{1 - np.exp(-DT/max(1e-3, sp["treatvs"]))!r}', 'treatvs (years to suppression): 1 - exp(-dt/treatvs) at dt = 0.2, at the start')
    for h, k in (('diagnosed', 'propdx'), ('in_care', 'propcare'), ('on_art', 'proptx'), ('suppressed', 'propsupp')):
        assert np.all(np.isnan(sp[k])), f'Given values of {k} are not handled'
        row('hold', h, 'const', f'{float(sp["fix" + k])!r}', f'fix{k}: {k} held from then at its value then')

    # Coverages Optima applies inside its births and force of infection
    row('pmtct', 'number', 'linear', series(tvec, sp['numpmtct']), 'numpmtct: diagnosed women giving birth with PMTCT, a year')
    row('pmtct', 'efficacy', 'const', f'{float(sp["effpmtct"])!r}', 'effpmtct')
    for k in ('prep', 'pep'):
        by_pop(k, 'coverage', sp[k], f'{k}: coverage by population')
        row(k, 'efficacy', 'const', f'{float(sp["eff" + k])!r}', f'eff{k}')
    return pd.DataFrame(rows)


if __name__ == '__main__':
    df = extract()
    df.to_csv(OUT, index=False)
    print(f'{len(df)} rows to {OUT}')
    print(df.groupby(['pathway', 'parameter']).size().to_string())
