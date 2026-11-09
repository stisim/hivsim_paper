"""Optima HIV Malawi (Malawi_20220912.prj): new infections by the source's cascade state, with the ever-treated tracer (env optima)"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import optima as op

HERE = Path(__file__).resolve().parent
REF = HERE.parents[1] / 'models' / 'optima' / 'reference'
sys.path.insert(0, str(HERE))
import optima_attr  # noqa: E402

END = 2030
STATES = dict(undx=range(2, 8), dx=range(8, 14), care=range(14, 20), lost=range(20, 26), usvl=range(26, 32), svl=range(32, 38))
CATS = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted_art']


def by_source(T, key):
    """ Infections per year by source, (rows: Optima's states, care and lost split into never/ever treated) × time steps """
    inc, ppl, ever = T[key], T['people'], T['ever']
    frac = np.where(ppl > 0, ever/np.where(ppl > 0, ppl, 1), 0) # Ever-treated share of each compartment
    rows = {}
    for s, idx in STATES.items():
        idx = list(idx)
        if s in ('care', 'lost'):
            rows[f'{s}_never'] = (inc[idx]*(1 - frac[idx])).sum(axis=(0, 1))
            rows[f'{s}_ever'] = (inc[idx]*frac[idx]).sum(axis=(0, 1))
        else:
            rows[s] = inc[idx].sum(axis=(0, 1))
    return pd.DataFrame(rows, index=T['tvec'])


MAPPINGS = dict(
    states=dict(undx='undiagnosed', dx='diagnosed_never_treated', care_never='diagnosed_never_treated', care_ever='diagnosed_never_treated',
                lost_never='interrupted_art', lost_ever='interrupted_art', usvl='on_art', svl='on_art'),
    history=dict(undx='undiagnosed', dx='diagnosed_never_treated', care_never='diagnosed_never_treated', care_ever='interrupted_art',
                 lost_never='diagnosed_never_treated', lost_ever='interrupted_art', usvl='on_art', svl='on_art'),
)


def annual(df, how):
    """ Optima's convention ('step': the rate in the step starting at each whole year) or the mean over the year's steps """
    yr = np.floor(df.index.values + 1e-9).astype(int)
    if how == 'step':
        keep = np.isclose(df.index.values, yr)
        return df[keep].set_index(yr[keep])
    return df.groupby(yr).mean()


if __name__ == '__main__':
    P = op.loadproj(str(REF / 'Malawi_20220912.prj'), verbose=0)
    base = P.runsim(name='20220912', end=END, verbose=0)
    op.project.model = optima_attr.model
    res = P.runsim(name='20220912', end=END, verbose=0)
    T = optima_attr.TRACE
    for k in ['numinci', 'numplhiv', 'numtreat', 'numdeath']:
        d = np.max(np.abs(res.main[k].tot[0] - base.main[k].tot[0]))
        print(f'patched vs original {k}: max abs diff {d:.3g}')
        assert d < 1e-6*np.max(base.main[k].tot[0])
    # Checks on the tracer: never above the compartment, and ever-treated stock reconciles
    ppl, ever = T['people'], T['ever']
    care, lost = list(STATES['care']), list(STATES['lost'])
    over = (ever[care + lost] - ppl[care + lost]).max()
    print(f'tracer: max (ever - people) in care/lost {over:.3g}; min ever {ever.min():.3g}')
    assert over < 1e-6 and ever.min() > -1e-6
    ti = np.argmin(abs(T['tvec'] - 2020))
    print(f"2020: in care {ppl[care, :, ti].sum():.0f} (ever treated {ever[care, :, ti].sum():.0f}); lost {ppl[lost, :, ti].sum():.0f} (ever treated {ever[lost, :, ti].sum():.0f})")
    # The total matches Optima's numinci
    tot = T['incibypop'].sum(axis=(0, 1))
    yrs = np.round(res.tvec).astype(int)
    print('numinci check (whole-year steps):', np.allclose(annual(pd.DataFrame(dict(t=tot), index=T['tvec']), 'step').t.loc[yrs].values, res.main['numinci'].tot[0]))
    out = []
    for counted, key in [('sexual', 'incibypop_sex'), ('all', 'incibypop')]:
        raw = by_source(T, key)
        for how in ['step', 'mean']:
            a = annual(raw, how)
            for c in a.columns:
                out.append(pd.DataFrame(dict(year=a.index, counted=counted, annual=how, mapping='raw', category=c, number=a[c].values)))
            for m, mp in MAPPINGS.items():
                g = a.T.groupby(pd.Series(mp)).sum().T
                for c in CATS:
                    out.append(pd.DataFrame(dict(year=g.index, counted=counted, annual=how, mapping=m, category=c, number=g[c].values)))
    df = pd.concat(out, ignore_index=True)
    df['share'] = 100*df.number/df.groupby(['year', 'counted', 'annual', 'mapping']).number.transform('sum')
    df.loc[df.mapping == 'raw', 'share'] = np.nan
    df['projection'] = df.year > 2020
    df.to_csv(HERE / 'outputs' / 'optima_by_source.csv', index=False)
    sel = df[(df.annual == 'step') & (df.mapping != 'raw') & df.year.isin([2000, 2005, 2010, 2015, 2020, 2024])]
    print(sel.pivot_table(index=['counted', 'mapping', 'year'], columns='category', values='share')[CATS].round(1).to_string())
    print(df[(df.annual == 'step') & (df.mapping == 'states')].groupby(['counted', 'year']).number.sum().unstack(0).loc[[2000, 2010, 2020, 2024, 2030]].round(0))
