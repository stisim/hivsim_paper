"""
The cascade's parameters as a table, one row per pathway parameter (and sex where it differs):

    pathway, parameter, sex, form, value, default, source

`sex` is all, f or m. `form` says how `value` reads:
- `const`: a number (a probability, or a count of days);
- `sigmoid`: "min, max, mid, rate", the value in year y being min + (max - min)/(1 + exp(-rate (y - mid)));
- `step`: "year: value, ..." constant from each year on, `default` before the first;
- `linear`: "year: value, ..." linear between years and constant after the last, `default` before the first;
- `exp_days`, `uniform_days`, `days`: a delay, "mean", "min, max", or "days".

Optional columns, for pathways by group and place in care (a missing column, or an empty cell, is `all` or none):
- `group`: who the row applies to, "key:value|value; ..." with keys `sex` (f, m), `age` ("lo-hi", years, hi
  exclusive), `cd4` (bands: acute, gt500, gt350, gt200, gt50, lt50) and `pop` (Optima populations, as named in
  models/optima/v4/poplabel.py); a person is in the group if they match every key;
- `places`: for rate pathways, "src|src>dst" over the places in care undx, dx, care, lost, svl and usvl (rates.py);
- `per`: for rate pathways given as a number a year, the places whose people the number is spread over.

Pathway types for a cascade given by rates and levels (Optima's), read by pathways.py (agents) and rates.py (compartments):
- `rate`: one row per rate (parameter: its name), a rate per person-year for those in `places`' sources and `group`,
  moving them to its destination; or, with `per`, people a year (at full scale), spread over everyone in `per`'s places;
- `target`: `on_art`, the number on ART by year, filled from those in care (`lowest_cd4_before`: lowest CD4 band first
  before that year, else in proportion; `p_suppressive`: the share starting suppressed) and removed back to care when over;
- `hold`: `diagnosed`, `in_care`, `on_art`, `suppressed`: the share held from the row's year (const) at its value then
  (diagnosed among PLHIV, in care or on ART among the diagnosed, on ART among those in care or on ART, suppressed among
  those on ART); a held share of on ART replaces the target;
- `pmtct` (`number`, `efficacy`), `prep` and `pep` (`coverage` by group, `efficacy`): coverages a model applies itself.

Flows (VMMC, PrEP) are a second table, one row per group and period: intervention, node, sex, age_lo, age_hi,
start, end, n (people over the period, at full scale), restrict.
"""
import numpy as np
import pandas as pd

COLUMNS = ['pathway', 'parameter', 'sex', 'form', 'value', 'default', 'source']
OPTIONAL = ['group', 'places', 'per'] # Pathways by group and place in care


def text(x):
    """ An optional cell: '' if missing """
    return '' if x is None or (isinstance(x, float) and np.isnan(x)) else str(x).strip()


def parse_group(g):
    """ {key: [values]} from "key:value|value; ..." ('' or 'all': everyone) """
    g = text(g)
    if g in ('', 'all'):
        return {}
    return {k.strip(): v.strip().split('|') for k, v in (part.split(':', 1) for part in g.split(';') if part.strip())}


def parse_places(p):
    """ (sources, destination) from "src|src>dst" """
    src, dst = text(p).split('>')
    return src.split('|'), dst


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
        for c in OPTIONAL:
            if c not in self.table:
                self.table[c] = ''
        self.fns = {}
        self.records = [] # Rows by group (rate, coverage): dicts with pathway, parameter, group, places, per, fn, source
        for r in self.table.itertuples():
            default = 0.0 if pd.isna(r.default) else float(r.default)
            fn = by_year(r.form, r.value, default) if r.form in ('const', 'sigmoid', 'step', 'linear') else None
            if r.pathway == 'rate' or parse_group(r.group):
                self.records.append(dict(pathway=r.pathway, parameter=r.parameter, sex=r.sex, group=parse_group(r.group), places=text(r.places), per=text(r.per), fn=fn, source=r.source))
            else:
                self.fns[(r.pathway, r.parameter, r.sex)] = (r.form, r.value, fn)

    def row(self, pathway, parameter, sex='all'):
        key = (pathway, parameter, sex)
        if key not in self.fns:
            key = (pathway, parameter, 'all')
        return self.fns[key]

    def has(self, pathway, parameter=None):
        return any(k[0] == pathway and (parameter is None or k[1] == parameter) for k in self.fns) or any(r['pathway'] == pathway and (parameter is None or r['parameter'] == parameter) for r in self.records)

    def rows(self, pathway):
        """ The rows by group of a pathway (see `records`) """
        return [r for r in self.records if r['pathway'] == pathway]

    def p(self, pathway, parameter, year, sex='all'):
        return float(np.clip(self.row(pathway, parameter, sex)[2](year), 0, 1))

    def value(self, pathway, parameter, year=None, sex='all'):
        return self.row(pathway, parameter, sex)[2](year)

    def delay(self, pathway, parameter):
        """ (kind, numbers) of a delay row """
        form, value, _ = self.row(pathway, parameter)
        return form, [float(x) for x in str(value).split(',')]
