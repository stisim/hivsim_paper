"""
Analysis 17: EMOD′ with its copied bugs fixed. See README.md.

    python run.py identity   # build(all EMOD′) to 2025 (set 0, seed 1), every result and table, saved to compare with the parent commit's
    python run.py literal    # reenrol in the condensed cascade against Daniel's event in the literal cascade (set 0, 5 seeds)
    python run.py runs       # every variant x sets 0, 25, 50 x 10 seeds through the runner (resumed from the 1985 burn-in)
    python run.py tables     # outputs/*.csv from the stored runs
    python run.py plot       # outputs/variants.png
"""
import os
import sys
import json
import pickle
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))
from common import swaps, targetfit # noqa: E402
from common.swaps import runner # noqa: E402

OUT = here / 'outputs'
PIDS, SEEDS, STOP = [0, 25, 50], 10, 2025
FIXES = ['mortality', 'reenrol', 'art_male', 'risk_nodes', 'anc', 'pmtct']
MODELS = dict(emod=None, emod_fixed='models/emod/v4/components_fixed.py', **{f'fix_{f}': f'analyses/17_emod_fixed/variants/{f}.py' for f in FIXES},
              fix_all='analyses/17_emod_fixed/variants/all.py') # Variant: factory file (None: EMOD′'s own)
LABELS = dict(emod='EMOD′ (bugs copied)', emod_fixed='Daniel\'s fixes (mortality + reenrol)', fix_all='all six fixes', **{f'fix_{f}': f'fix {f}' for f in FIXES})
WPP_YEARS = [1990, 2000, 2010, 2020]


def config(m):
    c = dict(hybrids=[{s: m for s in ['D', 'RN', 'V', 'C']}], pairs={m: PIDS}, seeds=SEEDS, stop=STOP)
    if MODELS[m]:
        c['models'] = {m: MODELS[m]}
    return c


def report_extras(sim):
    """ By year (EMOD's report at .92): the shares of PLHIV OffART (ARTStatus) and on ART, from EmodReport """
    r = sim.analyzers.emodreport.df
    r = r[(r.HasHIV == 1) & (r.Year % 1 > 0.5)].assign(year=lambda d: np.floor(d.Year).astype(int))
    g = r.groupby('year')
    tot = g.Population.sum()
    off = r[r.ARTStatus == 'OffART'].groupby('year').Population.sum().reindex(tot.index, fill_value=0)
    return pd.DataFrame(dict(offart_share=off/tot, report_art_coverage=g.On_ART.sum()/tot)).reset_index()


def run_one(job):
    m, cell = job
    cfg = config(m)
    choice = cfg['hybrids'][0]
    folder = runner.OUT / runner.hybrid_id(choice, cfg, dict(runner.code_version(), **runner.pins()))
    folder.mkdir(parents=True, exist_ok=True)
    extra = folder / f'{cell["cell"]}.extras.csv'
    if (folder / f'{cell["cell"]}.csv.gz').exists() and extra.exists():
        return m, cell['cell'], 'stored'
    try:
        (folder / 'hybrid.json').write_text(json.dumps(dict(name=swaps.name(choice), choice=choice, stop=STOP, models=cfg.get('models', {})), indent=1))
        sim = runner.start(choice, cell, cfg, folder)
        sim.run()
        runner.save(sim, folder, cell['cell'])
        report_extras(sim).to_csv(extra, index=False)
        return m, cell['cell'], ''
    except Exception as e:
        import traceback
        traceback.print_exc()
        return m, cell['cell'], repr(e)


def runs(cores=40):
    jobs = [(m, c) for m in MODELS for c in swaps.cells(config(m))]
    with Pool(cores) as pool:
        for m, c, err in pool.imap_unordered(run_one, jobs):
            print(m, c, err or 'ok', flush=True)


def load():
    """ Every stored run: the yearly outputs with the report's extras, and {tid: value} per run """
    ys, ts = [], []
    for m in MODELS:
        cfg = config(m)
        folder = runner.OUT / runner.hybrid_id(cfg['hybrids'][0], cfg, dict(runner.code_version(), **runner.pins()))
        for c in swaps.cells(cfg):
            f = folder / f'{c["cell"]}.csv.gz'
            if not f.exists():
                continue
            y = pd.read_csv(f).merge(pd.read_csv(folder / f'{c["cell"]}.extras.csv'), on='year', how='left')
            ys.append(y.assign(variant=m, pid=PIDS[c['pair']], seed=c['seed']))
            t = json.loads((folder / f'{c["cell"]}.targets.json').read_text())
            ts.append(dict(variant=m, pid=PIDS[c['pair']], seed=c['seed'], values={int(k): v for k, v in t.items()}))
    return pd.concat(ys, ignore_index=True), ts


def tables():
    OUT.mkdir(exist_ok=True)
    y, ts = load()
    targets = targetfit.load_targets()
    scale = 100 # EMOD′'s people per agent
    for k in ['pop', 'plhiv', 'on_art', 'hiv_deaths', 'new_infections']:
        y[k] = y[k]*scale
    keep = ['pop', 'plhiv', 'art_coverage', 'offart_share', 'tx_interrupted', 'tx_undiagnosed', 'tx_diagnosed_naive', 'tx_art_lt6m', 'tx_art_supp', 'tx_art_nonsupp',
            'new_infections', 'hiv_deaths', 'share_deaths_on_art', 'incidence_15_49_f', 'incidence_15_49_m']
    q = y.groupby(['variant', 'pid', 'year'])[keep]
    s = pd.concat({'mean': q.mean(), 'q10': q.quantile(0.1), 'q90': q.quantile(0.9)}, axis=1)
    s.columns = [f'{b}_{a}' for a, b in s.columns]
    s.reset_index().to_csv(OUT / 'yearly.csv', index=False, float_format='%.5g')
    # Population against WPP (all ages, both sexes)
    pop = targets[targets.target == 'population'].assign(year=lambda d: np.floor(d.year).astype(int)).groupby('year').value.sum()
    # Scores against the targets, and the effect table at 2010 and 2020
    rows = []
    for t in ts:
        rows.append(dict(variant=t['variant'], pid=t['pid'], seed=t['seed'], score=targetfit.score(t['values'], targets),
                         score_primary=targetfit.score(t['values'], targets, roles=('primary',))))
    sc_ = pd.DataFrame(rows)
    sc_.to_csv(OUT / 'scores.csv', index=False, float_format='%.4g')
    out = []
    for (m, pid), d in y.groupby(['variant', 'pid']):
        r = dict(variant=m, pid=pid, runs=d.seed.nunique(), took_off=int((d[d.year == 2000].plhiv > 0).sum()))
        for yr in (2010, 2020):
            e = d[d.year == yr]
            for k in ['pop', 'plhiv', 'art_coverage', 'offart_share', 'tx_interrupted', 'hiv_deaths', 'share_deaths_on_art']:
                r[f'{k}_{yr}'] = e[k].mean()
            r[f'pop_vs_wpp_{yr}'] = e['pop'].mean()/pop.get(yr, np.nan)
        r['tx_interrupted_2020_q10'], r['tx_interrupted_2020_q90'] = d[d.year == 2020].tx_interrupted.quantile([0.1, 0.9])
        g = sc_[(sc_.variant == m) & (sc_.pid == pid)]
        r['score'], r['score_primary'] = g.score.mean(), g.score_primary.mean()
        out.append(r)
    eff = pd.DataFrame(out)
    eff['order'] = eff.variant.map({m: i for i, m in enumerate(MODELS)})
    eff = eff.sort_values(['pid', 'order']).drop(columns='order')
    eff.to_csv(OUT / 'effects.csv', index=False, float_format='%.4g')
    # Bansi-Matharu 2025 Figure 3A: sources of new infections in four groups, and new infections (national, at 1% scale x 100)
    y['src_undiagnosed'], y['src_diagnosed_never_treated'] = y.tx_undiagnosed, y.tx_diagnosed_naive
    y['src_on_art'], y['src_interrupted'] = y.tx_art_lt6m + y.tx_art_supp + y.tx_art_nonsupp, y.tx_interrupted
    src = ['src_undiagnosed', 'src_diagnosed_never_treated', 'src_on_art', 'src_interrupted', 'new_infections']
    f = y[y.year.isin([2020, 2024])].groupby(['pid', 'variant', 'year'])[src]
    fig = f.mean().join(f.quantile(0.1).add_suffix('_q10')).join(f.quantile(0.9).add_suffix('_q90')).reset_index()
    fig['order'] = fig.variant.map({m: i for i, m in enumerate(MODELS)})
    fig = fig.sort_values(['pid', 'year', 'order']).drop(columns='order')
    fig.to_csv(OUT / 'fig3a.csv', index=False, float_format='%.4g')
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(fig[fig.pid == 0][['variant', 'year'] + src].round(3).to_string(index=False))
    # Against the UNAIDS series (2020)
    un = targets[(targets.source == 'UNAIDS 2025') & (np.floor(targets.year) == 2020)][['target', 'value']]
    print(un.to_string())
    print(pop.to_string())
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(eff[eff.pid == 0].round(3).to_string(index=False))


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s = pd.read_csv(OUT / 'yearly.csv')
    targets = targetfit.load_targets()
    targets = targets.assign(yr=np.floor(targets.year).astype(int))
    panels = [('pop', 'Population', None), ('plhiv', 'PLHIV', 'plhiv'), ('art_coverage', 'ART coverage (of PLHIV)', 'art_coverage'),
              ('offart_share', 'OffART share of PLHIV', None), ('tx_interrupted', 'Share of transmissions from interrupted ART', None),
              ('new_infections', 'New infections', 'new_infections'), ('hiv_deaths', 'HIV deaths', 'aids_deaths'), ('share_deaths_on_art', 'Share of HIV deaths on ART', None)]
    show = ['emod', 'fix_mortality', 'fix_reenrol', 'emod_fixed', 'fix_all', 'fix_art_male']
    for pid in sorted(s.pid.unique()):
        fig, axes = plt.subplots(2, 4, figsize=(18, 8))
        d = s[(s.pid == pid) & (s.year >= 1985)]
        for ax, (k, title, tg) in zip(axes.flat, panels):
            for m in show:
                e = d[d.variant == m]
                ax.plot(e.year, e[f'mean_{k}'], label=LABELS[m])
                ax.fill_between(e.year, e[f'q10_{k}'], e[f'q90_{k}'], alpha=0.15)
            if tg:
                t = targets[(targets.target == tg) & (targets.source == 'UNAIDS 2025')]
                ax.plot(t.yr, t.value, 'k.', label='UNAIDS')
            if k == 'pop':
                t = targets[targets.target == 'population'].groupby('yr').value.sum()
                ax.plot(t.index, t.values, 'ks', label='WPP')
            if k == 'tx_interrupted':
                ax.axhspan(0.084, 0.201, color='grey', alpha=0.2, label='MIHPSA range')
            ax.set_title(title)
            ax.spines[['top', 'right']].set_visible(False)
        axes.flat[0].legend(fontsize=7)
        fig.suptitle(f'EMOD′ variants, parameter set {pid}: mean and 10-90% over {SEEDS} seeds')
        fig.tight_layout()
        fig.savefig(OUT / f'variants_set{pid:03d}.png', dpi=110)


def identity(path):
    """ build(all EMOD′), set 0 seed 1, to 2025: results and analyzer tables, pickled to `path` """
    cfg = config('emod')
    cfg['pairs'], cfg['seeds'] = {'emod': [0]}, 1
    sim = swaps.build(cfg['hybrids'][0], swaps.cells(cfg)[0], cfg).run()
    res = {k: np.asarray(v) for k, v in sim.results.flatten().items()}
    an = {f'{k}.{a}': getattr(x, a) for k, x in sim.analyzers.items() for a in ('df', 'strata') if isinstance(getattr(x, a, None), pd.DataFrame)}
    pickle.dump(dict(res=res, an=an), open(path, 'wb'))


def literal_one(job):
    sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
    import harness
    import inputs
    cascade, fixes, seed = job
    sim = harness.make_sim(0, seed=inputs.load(0).run_number*10 + seed, cascade=True, condensed=cascade == 'condensed', record=False, fixes=fixes)
    sim.run()
    return report_extras(sim).assign(cascade=cascade, fixes='+'.join(fixes) or 'none', seed=seed)


def literal(cores=20):
    OUT.mkdir(exist_ok=True)
    jobs = [(c, f, s) for c in ('interpreted', 'condensed') for f in ((), ('reenrol',)) for s in range(1, 6)]
    with Pool(cores) as pool:
        df = pd.concat(pool.map(literal_one, jobs))
    df.to_csv(OUT / 'literal_vs_condensed.csv', index=False, float_format='%.5g')
    g = df[df.year.isin([2005, 2010, 2015, 2020, 2024])].groupby(['fixes', 'cascade', 'year'])[['offart_share', 'report_art_coverage']]
    t = pd.concat({'mean': g.mean(), 'sd': g.std()}, axis=1).unstack('cascade')
    print(t.round(3).to_string())


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'runs'
    if what == 'identity':
        identity(sys.argv[2])
    else:
        globals()[what]()
