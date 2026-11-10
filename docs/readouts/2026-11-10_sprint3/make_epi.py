"""
Epi comparison figures for the sprint 3 readout (10 Nov 2026): targets matrix, core indicators, Figure 3A rebuilt.

    python make_epi.py extract   # data/ extracts from analyses 17, 18, 19 and the Optima reference (run_hivsim.py makes data/hivsim_*)
    python make_epi.py           # figures/ from data/ and common/

Sources: EMOD binary with Daniel's two fixes (PR #4), 10 sets x 5 seeds (analysis 18, variant `both`); EMOD' with the same
fixes (`emod_fixed`, analysis 17; sets 0/25/50, 10 seeds each; agents x 100); our Optima (the 2022 project, not MIHPSA's;
models/optima/reference and analysis 19, sexual transmissions, inputs frozen at 2020 after 2020); HIVsim reference set
(models/hivsim/reference/pars.json, provisional), 20 seeds x 50k agents (run_hivsim.py). All EMOD versions are uncalibrated
after the fixes.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

here = Path(__file__).resolve().parent
repo = here.parents[2]
A = repo / 'analyses'
DATA, FIGS = here / 'data', here / 'figures'

# Colours (dataviz reference palette; validated: models all-pairs, stages adjacent and all-pairs, light surface)
EMOD, OPTIMA, HIVSIM = '#2a78d6', '#eb6834', '#1baf7a'
DATAC = '#52514e' # Data / published: neutral dark grey
INK, MUTED, GRID, SURFACE = '#0b0b0b', '#898781', '#e1e0d9', '#fcfcfb'
STAGES = ['undiagnosed', 'diagnosed_never_treated', 'on_art', 'interrupted_art']
STAGE_LABEL = dict(undiagnosed='Undiagnosed', diagnosed_never_treated='Diagnosed, never treated', on_art='On ART', interrupted_art='Interrupted ART')
STAGE_COL = dict(undiagnosed='#6250d6', diagnosed_never_treated='#eda100', on_art='#008300', interrupted_art='#e87ba4')
CENSUS_2018 = 17_563_749 # Daniel's scaling (analysis 18)
PUBLISHED_2024 = dict(EMOD=30_582, Optima=8_459) # Figure 3 caption

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.axisbelow': True, 'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
                     'axes.titlesize': 10, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'lines.linewidth': 2,
                     'legend.frameon': False, 'svg.fonttype': 'none'})


def save(fig, name):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(FIGS / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


# ---------------------------------------------------------------- extracts
def extract():
    DATA.mkdir(exist_ok=True)
    # EMOD' with Daniel's fixes (analysis 17): sets 0/25/50, mean and 10-90% over 10 seeds
    y = pd.read_csv(A / '17_emod_fixed' / 'outputs' / 'yearly.csv')
    y = y[(y.variant == 'emod_fixed') & (y.year >= 1985)].copy()
    for q in ('mean', 'q10', 'q90'):
        y[f'tx_on_art_{q}'] = y[f'tx_art_lt6m_{q}'] + y[f'tx_art_supp_{q}'] + y[f'tx_art_nonsupp_{q}']
    cols = ['pid', 'year'] + [f'{c}_{q}' for c in ['pop', 'plhiv', 'new_infections', 'hiv_deaths', 'art_coverage', 'tx_undiagnosed',
                                                   'tx_diagnosed_naive', 'tx_on_art', 'tx_interrupted'] for q in ('mean', 'q10', 'q90')]
    y[cols].to_csv(DATA / 'emod_prime_fixed_yearly.csv', index=False, float_format='%.5g')

    # EMOD binary with Daniel's fixes (analysis 18, `both`): per year median over 50 runs of the shares; totals at 15-64 census-scaled
    s = pd.read_csv(A / '18_fig3a_emod' / 'outputs' / 'shares.csv.gz')
    s = s[(s.variant == 'both') & (s.year >= 1990)]
    sh = s.groupby('year')[STAGES].median()
    sh = sh.div(sh.sum(axis=1), axis=0).add_prefix('share_')
    t = pd.read_csv(A / '18_fig3a_emod' / 'outputs' / 'totals.csv.gz')
    t = t[(t.variant == 'both') & (t.year >= 1990)].copy()
    t['ni'] = t.ni_15_64*CENSUS_2018/t.pop2018*2/t.halves # Annualised where the run misses steps of 2024
    tot = t.groupby('year').ni.agg(ni_median='median', ni_q05=lambda v: v.quantile(0.05), ni_q95=lambda v: v.quantile(0.95), runs='size')
    sh.join(tot).reset_index().to_csv(DATA / 'emod_binary_fixed_fig3a.csv', index=False, float_format='%.5g')

    # Our Optima (2022 project): totals, and sexual transmissions by source stage (analysis 19's mapping D8)
    o = pd.read_csv(repo / 'models' / 'optima' / 'reference' / 'optima_malawi_totals.csv')
    o['year'] = o.year.astype(int)
    b = pd.read_csv(A / '19_fig3a_optima' / 'outputs' / 'optima_by_source.csv')
    b = b[(b.counted == 'sexual') & (b.annual == 'step') & (b.mapping == 'history') & b.category.isin(STAGES)]
    b = b.pivot_table(index='year', columns='category', values='number')[STAGES]
    b = b.div(b.sum(axis=1), axis=0).add_prefix('share_')
    o = o.merge(b.reset_index(), on='year', how='left')
    o['projection'] = o.year > 2020
    o.to_csv(DATA / 'optima_2022_project.csv', index=False, float_format='%.6g')
    print('extracts written')


# ---------------------------------------------------------------- loaders
def load():
    d = {}
    d['emodp'] = pd.read_csv(DATA / 'emod_prime_fixed_yearly.csv')
    d['emodb'] = pd.read_csv(DATA / 'emod_binary_fixed_fig3a.csv').set_index('year')
    d['opt'] = pd.read_csv(DATA / 'optima_2022_project.csv').set_index('year')
    h = pd.read_csv(DATA / 'hivsim_yearly.csv')
    h = h[h.year <= 2024] # 2025 is a partial year
    h['tx_on_art'] = h.tx_art_lt6m + h.tx_art_supp + h.tx_art_nonsupp
    d['hiv'] = h
    d['targets'] = pd.read_csv(repo / 'common' / 'targets.csv')
    d['pub'] = pd.read_csv(repo / 'common' / 'fig3a_malawi.csv', comment='#')
    return d


def hivsim_shares(h):
    """ Pooled over seeds: transmissions by stage / all sexual transmissions """
    n = h.tx_n
    num = {st: (h[c]*n).groupby(h.year).sum() for st, c in zip(STAGES, ['tx_undiagnosed', 'tx_diagnosed_naive', 'tx_on_art', 'tx_interrupted'])}
    return pd.DataFrame(num).div(n.groupby(h.year).sum(), axis=0).add_prefix('share_')


# ---------------------------------------------------------------- 1. targets matrix
def targets_matrix():
    rows = ['Prevalence by age and sex', 'PLHIV', 'New infections', 'AIDS deaths', 'Incidence 15–49', 'FSW prevalence',
            'Number on ART / ART coverage', 'Viral suppression', 'Population']
    # From docs/targets.md (HIVsim's set; "The other two models" for EMOD and Optima)
    M = {
        'EMOD': ['fitted', 'not used', 'not used', 'not used', 'fitted', 'not used', 'fitted', 'not used', 'fitted'],
        'Optima': ['fitted', 'checked', 'checked', 'checked', 'not used', 'fitted', 'input', 'checked', 'checked'],
        'HIVsim': ['fitted', 'fitted', 'fitted', 'fitted', 'fitted', 'fitted', 'input · fitted', 'fitted', 'input'],
    }
    notes = {('EMOD', 6): 'number on ART by sex 2004–15;\nART coverage MPHIA 2015–16',
             ('EMOD', 8): 'WPP (in EMOD’s files)',
             ('EMOD', 0): 'DHS 2004/10/15; MPHIA 2015',
             ('EMOD', 4): 'MPHIA 2015–16, 2020–21',
             ('Optima', 0): 'by population group (hivprev)',
             ('Optima', 5): 'FSW population in hivprev',
             ('Optima', 6): 'number on ART is a programme input',
             ('Optima', 7): '2016 cascade, for comparison',
             ('Optima', 1): 'Spectrum, for comparison',
             ('HIVsim', 6): 'number on ART input; coverage fitted',
             ('HIVsim', 7): 'MPHIA 2020–21 (VLS among PLHIV)',
             ('HIVsim', 2): 'UNAIDS, from 2000',
             ('HIVsim', 8): 'demography is an input; WPP checked'}
    fill = {'fitted': '#52514e', 'checked': '#c3c2b7', 'input': '#f0efec', 'not used': SURFACE, 'input · fitted': '#52514e'}
    txt = {'fitted': SURFACE, 'checked': INK, 'input': INK, 'not used': MUTED, 'input · fitted': SURFACE}
    models = list(M)
    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.set_xlim(0, len(models)); ax.set_ylim(len(rows), 0)
    ax.grid(False)
    for j, m in enumerate(models):
        for i, r in enumerate(rows):
            v = M[m][i]
            ax.add_patch(plt.Rectangle((j + 0.02, i + 0.04), 0.96, 0.92, color=fill[v], ec=GRID if v in ('not used', 'input') else fill[v], lw=0.8,
                                       hatch='////' if v == 'input' else None))
            note = notes.get((m, i))
            ax.text(j + 0.5, i + (0.38 if note else 0.5), v, ha='center', va='center', color=txt[v], fontsize=10, fontweight='bold')
            if note:
                ax.text(j + 0.5, i + 0.7, note, ha='center', va='center', color=txt[v], fontsize=7, linespacing=0.95)
    ax.set_xticks(np.arange(len(models)) + 0.5, ['EMOD', 'Optima HIV', 'HIVsim'], fontsize=10, color=INK)
    ax.xaxis.tick_top(); ax.tick_params(length=0)
    ax.set_yticks(np.arange(len(rows)) + 0.5, rows, fontsize=9.5, color=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    for j, m in enumerate(models): # Model identity: a coloured bar under each header
        ax.add_patch(plt.Rectangle((j + 0.02, -0.12), 0.96, 0.08, color=[EMOD, OPTIMA, HIVSIM][j], clip_on=False))
    fig.suptitle('Only prevalence by age and sex is fitted by all three models; HIVsim is fitted to\n'
                 'UNAIDS PLHIV, new infections and AIDS deaths, which EMOD never used and Optima only compared against', x=0.02, ha='left', fontsize=11.5, fontweight='bold', y=1.04)
    fig.text(0.02, -0.02, 'fitted = in the calibration likelihood; checked = compared, not fitted; input = given to the model; not used.  '
             'Source: docs/targets.md. EMOD: its calibration files (weights not yet read).\nOptima: our 2022 project (MIHPSA’s not yet received). '
             'HIVsim: target set of docs/targets.md (reference parameter set is provisional, hivsim_malawi exp. 03).', fontsize=7.5, color=MUTED, ha='left')
    save(fig, 'targets_matrix')


# ---------------------------------------------------------------- 2. core indicators
def core_indicators(d):
    tg = d['targets']
    un = lambda t: tg[(tg.target == t) & (tg.source == 'UNAIDS 2025')].assign(yr=lambda x: np.floor(x.year)).set_index('yr').value
    wpp = tg[tg.target == 'population'].assign(yr=lambda x: np.floor(x.year)).groupby('yr').value.sum()
    ep = d['emodp']; e0 = ep[ep.pid == 0].set_index('year')
    eb = ep.groupby('year') # Band: over sets 0/25/50 (10-90% over seeds within each)
    h = d['hiv']; hq = h.groupby('year')
    o = d['opt']
    h['prev'] = h.plhiv/h['pop']
    panels = [('PLHIV (millions)', 'plhiv', 1e-6, un('plhiv'), 'numplhiv'),
              ('New infections per year (thousands)', 'new_infections', 1e-3, un('new_infections'), 'numinci'),
              ('AIDS deaths per year (thousands)', 'hiv_deaths', 1e-3, un('aids_deaths'), 'numdeath'),
              ('ART coverage (% of PLHIV)', 'art_coverage', 100, un('art_coverage'), None),
              ('HIV prevalence, all ages (%)', 'prev', 100, (un('plhiv')/wpp).dropna(), None),
              ('Population (millions)', 'pop', 1e-6, wpp[wpp.index <= 2021], 'popsize')]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6), sharex=True)
    for ax, (title, col, k, data, ocol) in zip(axes.flat, panels):
        ax.set_title(title, fontsize=9.5)
        # EMOD'
        if col == 'prev':
            em, lo, hi = e0.plhiv_mean/e0.pop_mean, eb.apply(lambda g: (g.plhiv_q10/g.pop_q10).min()), eb.apply(lambda g: (g.plhiv_q90/g.pop_q90).max())
        else:
            em, lo, hi = e0[f'{col}_mean'], eb[f'{col}_q10'].min(), eb[f'{col}_q90'].max()
        em, lo, hi = em[em.index <= 2024], lo[lo.index <= 2024], hi[hi.index <= 2024] # 2025 is a partial year
        ax.fill_between(lo.index, lo*k, hi*k, color=EMOD, alpha=0.15, lw=0)
        ax.plot(em.index, em*k, color=EMOD)
        # HIVsim
        ax.fill_between(hq[col].quantile(0.1).index, hq[col].quantile(0.1)*k, hq[col].quantile(0.9)*k, color=HIVSIM, alpha=0.18, lw=0)
        ax.plot(hq[col].median().index, hq[col].median()*k, color=HIVSIM)
        # Optima
        ov = {'art_coverage': o.numtreat/o.numplhiv, 'prev': o.numplhiv/o.popsize}.get(col, o[ocol] if ocol else None)
        ov = ov[ov.index <= 2024]
        ax.plot(ov[ov.index <= 2020].index, ov[ov.index <= 2020]*k, color=OPTIMA)
        ax.plot(ov[ov.index >= 2020].index, ov[ov.index >= 2020]*k, color=OPTIMA, ls=(0, (2, 1.5)))
        # Data
        ax.plot(data.index + 0.5, data*k, 'o', ms=4, color=DATAC, mfc=DATAC, mec=SURFACE, mew=0.6, zorder=5)
        if col == 'art_coverage':
            m = tg[(tg.target == 'art_coverage') & (tg.source == 'MPHIA 2020-21')]
            for _, r in m.iterrows():
                ax.errorbar(r.year + (0.25 if r.sex == 'f' else -0.25), r.value*k, yerr=[[(r.value - r.lower)*k], [(r.upper - r.value)*k]],
                            fmt='s', ms=4, color=DATAC, mfc=SURFACE, capsize=2, lw=1, zorder=6)
        ax.set_xlim(1990, 2025)
        ax.set_ylim(bottom=0)
    handles = [Line2D([], [], color=EMOD, label='EMOD′ + Daniel’s 2 fixes (uncalibrated; agents × 100, not census-rescaled; set 0, band: sets 0/25/50, 10–90% of seeds)'),
               Line2D([], [], color=OPTIMA, label='Optima, our 2022 project (dashed: inputs frozen at 2020)'),
               Line2D([], [], color=HIVSIM, label='HIVsim reference set, provisional (median, 10–90% of 20 seeds × 50k agents)'),
               Line2D([], [], color=DATAC, marker='o', ls='', label='UNAIDS 2025 / WPP (no bounds in targets.csv)'),
               Line2D([], [], color=DATAC, marker='s', mfc=SURFACE, ls='', label='MPHIA 2020–21 men / women, 95% CI')]
    fig.legend(handles=handles, loc='lower left', ncol=2, bbox_to_anchor=(0.01, -0.09), fontsize=8.2)
    o20, e20, h20, u20 = o.numinci[2020], e0.new_infections_mean[2020], hq.new_infections.median()[2020], un('new_infections')[2020]
    d20, hd20 = un('aids_deaths')[2020], hq.hiv_deaths.median()[2020]
    fig.suptitle(f'In 2020, new infections: EMOD′ {e20/1e3:.0f}k, Optima {o20/1e3:.0f}k, HIVsim {h20/1e3:.0f}k vs UNAIDS {u20/1e3:.0f}k; '
                 f'HIVsim has almost no AIDS deaths after 2015 ({hd20/1e3:.1f}k vs {d20/1e3:.0f}k)',
                 x=0.01, ha='left', fontsize=11.5, fontweight='bold', y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, 'core_indicators')
    return dict(e20=e20, o20=o20, h20=h20, u20=u20, hd20=hd20, d20=d20)


# ---------------------------------------------------------------- 3. Figure 3A rebuilt
def stack(ax, sh, title, y0=2000, y1=2024, proj=None, label_year=2020, note=None):
    sh = sh.loc[(sh.index >= y0) & (sh.index <= y1)]
    x = sh.index.values
    base = np.zeros(len(x))
    for st in STAGES:
        v = sh[f'share_{st}'].values*100
        ax.fill_between(x, base, base + v, color=STAGE_COL[st], lw=0, step=None)
        ax.plot(x, base + v, color=SURFACE, lw=1) # 1px surface gap between bands
        if label_year in sh.index:
            i = list(x).index(label_year)
            if v[i] >= 6:
                ax.text(label_year, base[i] + v[i]/2, f'{v[i]:.0f}', ha='center', va='center', fontsize=8, fontweight='bold',
                        color=INK if st in ('diagnosed_never_treated', 'interrupted_art') else SURFACE)
        base = base + v
    if proj is not None:
        ax.axvspan(proj, y1, color=SURFACE, alpha=0.45, lw=0)
        ax.text(proj + 0.2, 101, 'projection', fontsize=7, color=MUTED, va='bottom')
    ax.axvline(label_year, color=INK, lw=0.6, ls=':')
    ax.set_xlim(2000, 2024); ax.set_ylim(0, 100)
    ax.set_title(title, fontsize=9.5, pad=12)
    ax.grid(False)
    if note:
        ax.text(0, -0.2, note, transform=ax.transAxes, fontsize=7.2, color=MUTED, va='top')


def fig3a(d):
    pub = d['pub'].pivot_table(index=['model', 'year'], columns='category', values='share')/100
    pe = pub.loc['EMOD'][STAGES].add_prefix('share_')
    po = pub.loc['Optima HIV'][STAGES].add_prefix('share_')
    pe, po = pe.div(pe.sum(axis=1), axis=0), po.div(po.sum(axis=1), axis=0)
    eb = d['emodb'][[f'share_{s}' for s in STAGES]]
    ep = d['emodp']; e0 = ep[ep.pid == 0].set_index('year')
    epf = pd.DataFrame({f'share_{s}': e0[f'{c}_mean'] for s, c in zip(STAGES, ['tx_undiagnosed', 'tx_diagnosed_naive', 'tx_on_art', 'tx_interrupted'])})
    o = d['opt'][[f'share_{s}' for s in STAGES]].dropna()
    hs = hivsim_shares(d['hiv'])
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.6), sharey=True)
    stack(axes[0, 0], pe, 'Published EMOD', note='Bansi-Matharu 2025 Fig 3A (digitised)')
    stack(axes[0, 1], eb, 'EMOD binary + Daniel’s 2 fixes', note='analysis 18; median of 10 sets × 5 seeds;\nuncalibrated after the fixes')
    stack(axes[0, 2], epf, 'EMOD′ + Daniel’s 2 fixes', note='analysis 17; set 0, mean of 10 seeds;\nuncalibrated after the fixes')
    stack(axes[0, 3], hs, 'HIVsim (not in the paper)', note='reference set (provisional), 20 seeds × 50k;\npooled over seeds')
    stack(axes[1, 0], po, 'Published Optima HIV', note='Bansi-Matharu 2025 Fig 3A (digitised)')
    stack(axes[1, 1], o, 'Our Optima: 2022 project', proj=2020, note='not MIHPSA’s project (awaiting Rowan); analysis 19;\nsexual transmissions; inputs frozen after 2020')
    for ax in axes[:, 0]:
        ax.set_ylabel('% of new infections by source')
    # Panel slots 6-7: legend and reading notes
    for ax in axes[1, 2:]:
        ax.axis('off')
    axes[1, 2].legend(handles=[Patch(color=STAGE_COL[s], label=STAGE_LABEL[s]) for s in STAGES][::-1], loc='upper left',
                      title='Source partner’s stage (top to bottom)', title_fontsize=9, fontsize=9, alignment='left')
    r = lambda df, s, y=2020: 100*df.loc[y, f'share_{s}']
    txt = ('2020 shares, undiagnosed / dx never treated / on ART / interrupted:\n'
           + '\n'.join(f'{n:<16} ' + ' / '.join(f'{r(df, s):.0f}' for s in STAGES)
                       for n, df in [('Published EMOD', pe), ('EMOD binary', eb), ('EMOD′', epf), ('Published Optima', po), ('Our Optima', o), ('HIVsim', hs)])
           + '\nNumbers in the bands are 2020 (dotted line).')
    axes[1, 3].text(-0.05, 1.0, txt, transform=axes[1, 3].transAxes, fontsize=7.6, va='top', family='monospace', color=INK)
    fig.suptitle('With Daniel’s fixes EMOD reproduces the published source split (within ~5 points); our Optima (2022 project) does not:\n'
                 f'it has {r(o, "diagnosed_never_treated"):.0f}% from diagnosed-never-treated in 2020 vs {r(po, "diagnosed_never_treated"):.0f}% published. '
                 f'HIVsim is {r(hs, "undiagnosed"):.0f}% undiagnosed, never-treated ~0 after 2009 and interrupted ~0 until 2021 ({r(hs, "interrupted_art", 2024):.0f}% in 2024)',
                 x=0.01, ha='left', fontsize=11.5, fontweight='bold')
    fig.subplots_adjust(left=0.06, right=0.98, top=0.84, bottom=0.1, wspace=0.18, hspace=0.62)
    save(fig, 'fig3a_reconstruction')
    return dict(pe=pe, eb=eb, epf=epf, po=po, o=o, hs=hs)


def fig3a_totals(d):
    eb = d['emodb']; eb = eb[(eb.index >= 2000) & (eb.index <= 2024)]
    ep = d['emodp']; e0 = ep[ep.pid == 0].set_index('year'); e0 = e0[(e0.index >= 2000) & (e0.index <= 2024)]
    o = d['opt']; o = o[(o.index >= 2000) & (o.index <= 2024)]
    h = d['hiv']; hq = h[h.year >= 2000].groupby('year').new_infections
    tg = d['targets']; un = tg[(tg.target == 'new_infections') & (tg.source == 'UNAIDS 2025') & (tg.year >= 2000)]
    k = 1e-3
    fig, ax = plt.subplots(figsize=(11, 5.6))
    ax.fill_between(eb.index, eb.ni_q05*k, eb.ni_q95*k, color=EMOD, alpha=0.13, lw=0)
    ax.plot(eb.index, eb.ni_median*k, color=EMOD, label='EMOD binary + Daniel’s 2 fixes (ages 15–64, census-scaled; median, 5–95% of 50 runs)')
    ax.plot(e0.index, e0.new_infections_mean*k, color=EMOD, ls=(0, (1, 1.2)), lw=2, label='EMOD′ + Daniel’s 2 fixes (all ages, agents × 100; set 0)')
    ax.plot(o.index[o.index <= 2020], o.numinci[o.index <= 2020]*k, color=OPTIMA, label='Optima, our 2022 project (all ages; dashed: inputs frozen at 2020)')
    ax.plot(o.index[o.index >= 2020], o.numinci[o.index >= 2020]*k, color=OPTIMA, ls=(0, (2, 1.5)))
    ax.fill_between(hq.median().index, hq.quantile(0.1)*k, hq.quantile(0.9)*k, color=HIVSIM, alpha=0.18, lw=0)
    ax.plot(hq.median().index, hq.median()*k, color=HIVSIM, label='HIVsim reference set (all ages; median, 10–90% of 20 seeds)')
    ax.plot(un.year, un.value*k, 'o', ms=4.5, color=DATAC, mec=SURFACE, mew=0.6, label='UNAIDS 2025', zorder=5)
    for m, c, xy in (('EMOD', EMOD, (2024.3, 68)), ('Optima', OPTIMA, None)):
        v = PUBLISHED_2024[m]
        dy = xy is not None
        ax.plot(2024, v*k, 'D', ms=9, color=c, mec=INK, mew=1.2, zorder=6)
        ax.annotate(f'Published {m} 2024: {v:,}', (2024, v*k), xytext=xy or (-10, 0), textcoords='data' if dy else 'offset points',
                    ha='right', va='center',
                    fontsize=8.5, color=INK, fontweight='bold', bbox=dict(fc=SURFACE, ec='none', pad=1),
                    arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8) if dy else None, zorder=7)
    ax.plot([], [], 'D', ms=7, color=MUTED, mec=INK, mew=1.2, label='Published 2024 totals (Fig 3 caption)')
    ax.set_xlim(2000, 2024.8); ax.set_ylim(0, None)
    ax.set_ylabel('New infections per year (thousands)')
    ax.legend(loc='upper left', bbox_to_anchor=(0, -0.08), ncol=2, fontsize=7.8)
    e24, o24, h24 = eb.ni_median[2024], o.numinci[2024], hq.median()[2024]
    ax.set_title(f'2024 new infections: EMOD binary {e24/1e3:.1f}k vs published 30.6k; our Optima {o24/1e3:.1f}k vs published 8.5k;\n'
                 f'HIVsim {h24/1e3:.1f}k. Ours are {e24/PUBLISHED_2024["EMOD"]:.2f}× (EMOD) and {o24/PUBLISHED_2024["Optima"]:.1f}× (Optima) the published totals', fontsize=11, loc='left')
    fig.tight_layout()
    save(fig, 'fig3a_totals')
    return dict(e24=e24, o24=o24, h24=h24, ep24=e0.new_infections_mean[2024])


if __name__ == '__main__':
    if 'extract' in sys.argv:
        extract()
    else:
        d = load()
        targets_matrix()
        print(core_indicators(d))
        fig3a(d)
        print(fig3a_totals(d))
