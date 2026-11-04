"""
The cascade's parameters as a table, one row per pathway parameter (and sex where it differs):

    pathway, parameter, sex, form, value, default, source

`sex` is all, f or m. `form` says how `value` reads:
- `const`: a number (a probability, or a count of days);
- `sigmoid`: "min, max, mid, rate", the value in year y being min + (max - min)/(1 + exp(-rate (y - mid)));
- `step`: "year: value, ..." constant from each year on, `default` before the first;
- `linear`: "year: value, ..." linear between years and constant after the last, `default` before the first;
- `exp_days`, `uniform_days`, `days`: a delay, "mean", "min, max", or "days".

Flows (VMMC, PrEP) are a second table, one row per group and period: intervention, node, sex, age_lo, age_hi,
start, end, n (people over the period, at full scale), restrict.
"""
import numpy as np
import pandas as pd

COLUMNS = ['pathway', 'parameter', 'sex', 'form', 'value', 'default', 'source']


def points(text):
    pairs = [p.split(':') for p in str(text).split(',') if p.strip()]
    return np.array([float(a) for a, b in pairs]), np.array([float(b) for a, b in pairs])


def by_year(form, value, default=0.0):
    """ A function of the year for one row """
    if form == 'const':
        v = float(value)
        return lambda year: v
    if form == 'sigmoid':
        lo, hi, mid, rate = (float(x) for x in str(value).split(','))
        return lambda year: lo + (hi - lo)/(1 + np.exp(-rate*(year - mid)))
    t, v = points(value)
    if form == 'step':
        return lambda year: default if year < t[0] else float(v[np.searchsorted(t, year, side='right') - 1])
    if form == 'linear':
        return lambda year: default if year < t[0] else float(np.interp(year, t, v))
    raise ValueError(f'form {form} is not a function of the year')


class Params:
    """ The table, read: `p(pathway, parameter, year, sex='all')` for values by year, `delay(...)` for delays """
    def __init__(self, table, flows=None):
        self.table = table if isinstance(table, pd.DataFrame) else pd.read_csv(table)
        self.flows = flows if flows is None or isinstance(flows, pd.DataFrame) else pd.read_csv(flows)
        self.fns = {}
        for r in self.table.itertuples():
            default = 0.0 if pd.isna(r.default) else float(r.default)
            self.fns[(r.pathway, r.parameter, r.sex)] = (r.form, r.value, by_year(r.form, r.value, default) if r.form in ('const', 'sigmoid', 'step', 'linear') else None)

    def row(self, pathway, parameter, sex='all'):
        key = (pathway, parameter, sex)
        if key not in self.fns:
            key = (pathway, parameter, 'all')
        return self.fns[key]

    def has(self, pathway, parameter):
        return any(k[:2] == (pathway, parameter) for k in self.fns)

    def p(self, pathway, parameter, year, sex='all'):
        return float(np.clip(self.row(pathway, parameter, sex)[2](year), 0, 1))

    def value(self, pathway, parameter, year=None, sex='all'):
        return self.row(pathway, parameter, sex)[2](year)

    def delay(self, pathway, parameter):
        """ (kind, numbers) of a delay row """
        form, value, _ = self.row(pathway, parameter)
        return form, [float(x) for x in str(value).split(',')]
