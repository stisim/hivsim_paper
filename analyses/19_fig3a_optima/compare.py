"""Our Optima (and Optima′) against the published Optima HIV panel of Figure 3A (common/fig3a_malawi.csv): table and figure"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CATS = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted_art']
LABELS = ['Undiagnosed', 'Diagnosed, treatment naive', 'Receiving treatment', 'Not receiving treatment after starting']
COLORS = ['#8fd3e8', '#f6b321', '#9cc96b', '#cd5c5c'] # As in the paper's figure
STATE_CAT = dict(undx='undiagnosed', dx='diagnosed_never_treated', care='diagnosed_never_treated', lost='interrupted_art', usvl='on_art', svl='on_art')


def shares(df):
    return 100*df[CATS].div(df[CATS].sum(axis=1), axis=0)


if __name__ == '__main__':
    pub = pd.read_csv(HERE.parents[1] / 'common' / 'fig3a_malawi.csv', comment='#')
    pub = pub[pub.model == 'Optima HIV'].pivot(index='year', columns='category', values='share')[CATS]
    o = pd.read_csv(HERE / 'outputs' / 'optima_by_source.csv')
    o = o[(o.annual == 'step') & (o.mapping != 'raw')]
    ours = {(c, m): o[(o.counted == c) & (o.mapping == m)].pivot(index='year', columns='category', values='number')[CATS] for c in ['sexual', 'all'] for m in ['history', 'states']}
    panels = {'Published (Figure 3A)': pub, 'Our Optima, sexual': shares(ours['sexual', 'history']), 'Our Optima, incl. MTCT': shares(ours['all', 'history'])}
    pf = HERE / 'outputs' / 'prime_by_source.csv'
    if pf.exists():
        p = pd.read_csv(pf, index_col='year')
        p = p.T.groupby(pd.Series(STATE_CAT)).sum().T.reindex(columns=CATS)
        panels["Optima′, sexual (states)"] = shares(p)
    rows = []
    for name, d in panels.items():
        for y in [2000, 2005, 2010, 2015, 2020, 2024]:
            if y in d.index:
                rows.append(dict(source=name, year=y, **d.loc[y].round(1).to_dict()))
    for (c, m), d in ours.items():
        if m == 'states':
            for y in [2000, 2020, 2024]:
                rows.append(dict(source=f'Our Optima, {"sexual" if c == "sexual" else "incl. MTCT"} (states: lost = interrupted)', year=y, **shares(d).loc[y].round(1).to_dict()))
    tab = pd.DataFrame(rows)
    tab.to_csv(HERE / 'outputs' / 'shares_vs_published.csv', index=False)
    print(tab.to_string(index=False))
    tot = pd.DataFrame({'sexual': ours['sexual', 'history'].sum(axis=1), 'all': ours['all', 'history'].sum(axis=1)})
    print('\nOur Optima, new infections (step at each whole year):\n', tot.loc[[2000, 2005, 2010, 2015, 2018, 2020, 2022, 2024]].round(0))

    fig, axs = plt.subplots(len(panels), 1, figsize=(9, 2.6*len(panels)), sharex=True)
    for ax, (name, d) in zip(axs, panels.items()):
        d = d.loc[2000:2040]
        bottom = np.zeros(len(d))
        for c, lab, col in zip(CATS, LABELS, COLORS):
            ax.bar(d.index, d[c], bottom=bottom, color=col, width=0.7, label=lab)
            bottom += d[c].values
        ax.set_title(name, fontsize=10, loc='left')
        ax.set_ylim(0, 100); ax.set_ylabel('% of new infections')
        ax.axvline(2020.5, color='k', lw=0.8, ls=':')
    axs[0].legend(ncol=2, fontsize=8, loc='lower left')
    axs[-1].set_xlabel('Year (dotted: end of treatment data; our runs after 2020 are projections with inputs held)')
    fig.tight_layout()
    fig.savefig(HERE / 'figures' / 'fig3a_optima.png', dpi=110)
