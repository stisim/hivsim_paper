"""
Exact comparison of two runs: every result array of `sim.results`, and the tables the analyzers keep (`df`), and
EMOD′'s infection and death logs. For the identity tests (tests/test_components.py) and the burn-in check (analysis 11).
"""
import numpy as np
import pandas as pd


def arrays(sim):
    """ name -> array, for every result and every table the run keeps """
    out = {k: np.asarray(v) for k, v in sim.results.flatten(only_results=True).items()}
    for mod in list(sim.analyzers.values()) + list(sim.diseases.values()) + list(sim.networks.values()):
        for attr in ('df', 'infections', 'deaths', 'starts', 'ends'):
            v = getattr(mod, attr, None)
            if isinstance(v, pd.DataFrame):
                for c in v.columns:
                    out[f'{mod.name}.{attr}.{c}'] = v[c].to_numpy()
    return out


def differences(a, b):
    """ The names whose arrays differ between two runs (or that are in only one), with the largest difference; empty if identical """
    a, b = (x if isinstance(x, dict) else arrays(x) for x in (a, b))
    rows = []
    for k in sorted(set(a) | set(b)):
        if k not in a or k not in b:
            rows.append(dict(name=k, why='only in ' + ('second' if k not in a else 'first')))
        elif a[k].shape != b[k].shape:
            rows.append(dict(name=k, why=f'shape {a[k].shape} vs {b[k].shape}'))
        elif not np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == 'f'):
            d = np.abs(a[k].astype(float) - b[k].astype(float)) if a[k].dtype.kind in 'biuf' else None
            rows.append(dict(name=k, why='values', max_abs=np.nanmax(d) if d is not None else np.nan,
                             first=int(np.argmax(a[k] != b[k])) if a[k].ndim == 1 else -1))
    return pd.DataFrame(rows, columns=['name', 'why', 'max_abs', 'first'])
