"""
Measure the Optima HIV Malawi reference run against common/targets.csv and draw the shared target figures
(common/targetfit.py). Runs in the `hivsim_paper` (or `starsim`) env: it reads the committed `optima_malawi_ref.obj`
(written by export_ref.py in the `optima` env), not Optima itself.

    python optima_targets.py   # -> analyses/01_reference_fits/{outputs,figures}/optima_*

Optima is deterministic, so this is one run (parset 20220912). Everything is read from the raw outputs at dt = 0.2:
stocks are interpolated to the target year; flows (new infections, including mother-to-child, and HIV deaths) are
Optima's annualised rates averaged over the five steps of the calendar year. Incidence is new infections over
susceptible people, both averaged over the year to the survey.

Populations are F/M by 5-year age band (0-4 to 85+), plus three key populations aged 15-49 with no age band: FSW
(women), Clients and MSM (men). For targets by age band, each key population is spread over the 15-49 bands of its
sex in proportion to the general population of each band at each time step, so totals are unchanged. FSW
prevalence is the FSW population's.

Care states: diagnosed = diagnosed, in care, lost to follow-up, or on ART; on ART = on ART suppressed or unsuppressed;
suppressed = on ART suppressed.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import sciris as sc

here = Path(__file__).resolve().parent
root = here.parents[2]
sys.path.insert(0, str(root / 'common'))
import targetfit  # noqa: E402


def by_band(ref):
    """ Population, PLHIV, diagnosed, on ART, suppressed, susceptible, new infections and HIV deaths for each general
    population (sex, age band), with the key populations spread over the 15-49 bands of their sex """
    raw, s = ref['raw'], ref['settings']
    people = np.asarray(raw['people'])  # (state, pop, time)
    plhiv = s['allplhiv']
    dx = np.concatenate([s[k] for k in ['dx', 'care', 'lost', 'usvl', 'svl']])
    x = dict(population=people.sum(0), plhiv=people[plhiv].sum(0), diagnosed=people[dx].sum(0),
             on_art=people[s['alltx']].sum(0), suppressed=people[s['svl']].sum(0), susceptible=people[s['sus']].sum(0),
             new_infections=np.asarray(raw['inci']), aids_deaths=np.asarray(raw['death']).sum(0))
    keys = list(raw['popkeys'])
    gen = [k for k in keys if re.match(r'[FM]\d', k)]
    sex = np.array([k[0].lower() for k in gen])
    lo = np.array([int(re.findall(r'\d+', k)[0]) for k in gen])
    hi = np.array([int(k.split('-')[1]) + 1 if '-' in k else 100 for k in gen])
    g = [keys.index(k) for k in gen]
    out = {k: v[g].copy() for k, v in x.items()}
    for kp, ks in [('FSW', 'f'), ('Clients', 'm'), ('MSM', 'm')]:
        band = (sex == ks) & (lo >= 15) & (hi <= 50)
        w = x['population'][g][band] / x['population'][g][band].sum(0)
        for k, v in x.items():
            out[k][band] += w * v[keys.index(kp)]
    return out, sex, lo, hi, x, keys


def measure(ref, targets):
    out, sex, lo, hi, x, keys = by_band(ref)
    t = np.asarray(ref['raw']['tvec'])
    fsw = {k: v[keys.index('FSW')][None] for k, v in x.items()}
    vals = {}
    for tid, r in targets.iterrows():
        if r.group == 'fsw':
            d = fsw
        elif r.group == 'all':
            m = (lo >= r.age_lo) & (hi <= r.age_hi) & ((sex == r.sex) | (r.sex == 'all'))
            d = {k: v[m] for k, v in out.items()}
        else:
            raise ValueError(r.group)
        at = lambda k: np.interp(r.year, t, d[k].sum(0))
        if r.target in ('new_infections', 'aids_deaths'):
            vals[tid] = d[r.target].sum(0)[np.floor(t) == np.floor(r.year)].mean()
        elif r.target == 'incidence':
            yr = (t > r.year - 1) & (t <= r.year)
            vals[tid] = d['new_infections'].sum(0)[yr].mean() / d['susceptible'].sum(0)[yr].mean()
        else:
            vals[tid] = {
                'prevalence': lambda: at('plhiv') / at('population'),
                'plhiv': lambda: at('plhiv'),
                'on_art': lambda: at('on_art'),
                'population': lambda: at('population'),
                'art_coverage': lambda: at('on_art') / at('plhiv'),
                'diagnosed': lambda: at('diagnosed') / at('plhiv'),
                'on_ART_among_diagnosed': lambda: at('on_art') / at('diagnosed'),
                'suppressed_among_on_ART': lambda: at('suppressed') / at('on_art'),
                'vls_among_plhiv': lambda: at('suppressed') / at('plhiv'),
            }[r.target]()
    return {int(k): float(v) for k, v in vals.items()}


def main():
    targets = targetfit.load_targets()
    ref = sc.load(here / 'optima_malawi_ref.obj')
    expdir = root / 'analyses' / '01_reference_fits'
    (expdir / 'outputs').mkdir(parents=True, exist_ok=True)
    row = dict(i=0, kind='sample', pars=dict(parset='20220912'), values=measure(ref, targets))
    (expdir / 'outputs' / 'optima_results.jsonl').write_text(json.dumps(row) + '\n')
    targetfit.plot(expdir, 'optima')


if __name__ == '__main__':
    main()
