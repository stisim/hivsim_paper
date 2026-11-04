"""
Analysis 07: HIVsim with its own care against HIVsim with EMOD's condensed cascade. See README.md.
Needs Starsim v4 head and STIsim branch rs-cascade-host on PYTHONPATH.

    python run.py          # 50 seeds × (own, condensed), 20 × condensed_nopmtct, to outputs/
    python run.py --plot   # figures/
"""
import sys
import json
import time
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd
import starsim as ss

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path[:0] = [str(repo / 'models' / 'hivsim')]
import swap # noqa: E402
sys.path.insert(0, str(repo / 'common'))
import targetfit # noqa: E402  (common's, for scoring and plotting)
import importlib.util # noqa: E402
spec = importlib.util.spec_from_file_location('hivfit', repo / 'models' / 'hivsim' / 'malawi' / 'targetfit.py')
hivfit = importlib.util.module_from_spec(spec) # hivsim_malawi's, which measures the common targets in a HIVsim run
spec.loader.exec_module(hivfit)

SEEDS = dict(own=range(1, 51), condensed=range(1, 51), condensed_nopmtct=range(1, 21))
CARE = ['none', 'uptake', 'testing', 'staging', 'preart', 'linking', 'on_art', 'out', 'lost']


class Check(ss.Analyzer):
    """ Yearly, the host's events and states against HIVsim's own (people with the condensed cascade) """
    def __init__(self, **kwargs):
        super().__init__(name='check', **kwargs)
        self.rows, self.acc = [], {}

    def step(self):
        sim, ti = self.sim, self.ti
        hiv, art, ppl, preg = sim.diseases.hiv, sim.interventions.art, sim.people, sim.demographics.pregnancy
        cx = sim.interventions.get('cascade')
        a = self.acc
        add = lambda k, n: a.__setitem__(k, a.get(k, 0) + int(n))
        add('new_infections', len(hiv.new_infected().uids))
        add('art_starts', np.count_nonzero((art.ti_art.raw == ti) & art.on_treatment.raw))
        add('art_stops', np.count_nonzero(art.ti_stop_art.raw == ti) - np.count_nonzero((art.ti_stop_art.raw == ti) & art.on_treatment.raw))
        add('conceptions', np.count_nonzero(preg.ti_P.raw[preg.P.uids] == ti - 1)) # Conceptions come after the analyzers
        add('births', np.count_nonzero((ppl.age.raw[ppl.auids] >= 0) & (ppl.age.raw[ppl.auids] < float(sim.t.dt_year)))) # Born in the step
        if cx is not None:
            h = cx.host
            add('ev_debut', len(h.debuts()))
            add('ev_symptoms', len(h.newly_symptomatic()))
            add('ev_twelve_weeks', len(h.twelve_weeks_pregnant()))
            add('ev_six_weeks', len(h.six_weeks_old()))
        if (ti + 1) % 12: # Yearly, at the end of the year
            return
        u = (hiv.infected & ppl.alive).uids
        on = art.on_treatment.raw[u]
        row = dict(year=int(round(self.t.now('year'))), plhiv=len(u), on_art=int(on.sum()), diagnosed=int(hiv.diagnosed.raw[u].sum()),
                   effective=int(art.on_effective_art.raw[u].sum()), off_after_art=int(art.art_discontinued.raw[u].sum()),
                   pmtct_mothers=int(np.count_nonzero(hiv.ti_pmtct_stop.raw[ppl.auids] > ti)), on_prep=int(hiv.on_prep.raw[ppl.auids].sum()),
                   circ_15_49=float(np.mean(hiv.circumcised.raw[((~ppl.female) & (ppl.age >= 15) & (ppl.age < 50)).uids])), **a)
        if cx is not None:
            care = cx.care.raw[u]
            row |= {f'care_{c}': int(np.sum(care == k)) for k, c in enumerate(CARE)}
            row |= dict(bad_on_art_not_in_care=int(np.sum(on & (care != CARE.index('on_art')))), bad_in_care_not_on_art=int(np.sum(~on & (care == CARE.index('on_art')))),
                        bad_on_art_undiagnosed=int(np.sum(on & ~hiv.diagnosed.raw[u])), bad_staged_undiagnosed=int(np.sum((care >= CARE.index('staging')) & ~hiv.diagnosed.raw[u])))
        self.rows.append(row)
        self.acc = {}


def table_without_pmtct():
    P = swap.emod_table()
    t = P.table.copy()
    t.loc[(t.pathway == 'pmtct') & t.parameter.str.startswith('eff_'), 'value'] = '0'
    return swap.Params(t, P.flows)


def run_one(job):
    variant, seed = job
    params = table_without_pmtct() if variant == 'condensed_nopmtct' else None
    sim = swap.make_sim('own' if variant == 'own' else 'condensed', seed=seed, params=params, analyzers=[hivfit.TargetRecorder(targets=targetfit.load_targets()), Check()])
    t = time.time()
    sim.run()
    secs = time.time() - t
    r, years = sim.results.hiv, np.floor(np.asarray(sim.t.yearvec) + 1e-6).astype(int)
    flows = pd.DataFrame(dict(year=years, deaths=np.asarray(r.new_deaths), mtct=np.asarray(r.new_infections_mtct), diagnoses=np.asarray(r.new_diagnoses))).groupby('year').sum()/sim.pars.pop_scale # In agents
    check = pd.DataFrame(sim.analyzers.check.rows).set_index('year').join(flows).reset_index()
    tag = dict(model=variant, seed=seed)
    return dict(tag=tag, secs=secs, values=hivfit.measure(sim, targetfit.load_targets()), check=check.assign(scale=float(sim.pars.pop_scale), **tag))


def reduce(cores):
    jobs = [(v, s) for v, seeds in SEEDS.items() for s in seeds]
    with Pool(cores) as pool:
        res = pool.map(run_one, jobs, chunksize=1)
    out = here / 'outputs'
    out.mkdir(exist_ok=True)
    for v in SEEDS:
        rows = [r for r in res if r['tag']['model'] == v]
        (out / f'{v}_results.jsonl').write_text(''.join(json.dumps(dict(i=i, kind='sample', pars=dict(seed=r['tag']['seed']), values=r['values'])) + '\n' for i, r in enumerate(rows)))
    pd.concat([r['check'] for r in res], ignore_index=True).to_csv(out / 'series.csv.gz', index=False)
    pd.DataFrame([dict(**r['tag'], secs=r['secs']) for r in res]).to_csv(out / 'runtime.csv', index=False)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_muted
    figdir, out = here / 'figures', here / 'outputs'
    figdir.mkdir(exist_ok=True)
    for m in ('own', 'condensed'):
        targetfit.plot(here, m)
    note = 'Lines: median over 50 seeds; bands: 5-95%. HIVsim at 10,000 agents, numbers scaled to the population.'

    def save(fig, name, text=note):
        fig.text(0.01, -0.03, text, fontsize=7, color=c_muted, va='top')
        fig.savefig(figdir / name, dpi=150, bbox_inches='tight')
        plt.close(fig)

    s = pd.read_csv(out / 'series.csv.gz')
    s = s[(s.year >= 1990) & (s.year < 2022)] # 2022 has one step
    for c in ('plhiv', 'on_art', 'deaths', 'new_infections', 'mtct'):
        s[c + '_n'] = s[c]*s.scale
    s['dx_share'] = s.diagnosed/s.plhiv
    s['art_cov'] = s.on_art/s.plhiv
    s['off_share'] = s.off_after_art/s.plhiv
    s['supp_share'] = s.effective/s.on_art.where(s.on_art > 0)
    looks = dict(own=dict(color=c_best, label='HIVsim, own care'), condensed=dict(color=c_run, label='HIVsim, condensed cascade'))
    ref = pd.read_csv(repo / 'analyses' / '06_condensed_cascade' / 'outputs' / 'series.csv.gz')
    ref = ref[ref.model == 'condensed'].assign(dx_share=lambda d: d.diagnosed/d.plhiv, art_cov=lambda d: d.on_art/d.plhiv, off_share=lambda d: d.art_OffART/d.plhiv,
                                              supp_share=lambda d: d.art_OnEffectiveART/d.on_art.where(d.on_art > 0))
    panels = [('plhiv_n', 'People with HIV'), ('on_art_n', 'On ART'), ('new_infections_n', 'New infections a year'), ('deaths_n', 'HIV deaths a year'),
              ('dx_share', 'Diagnosed, share of PLHIV'), ('art_cov', 'On ART, share of PLHIV'), ('off_share', 'Off ART after starting it, share of PLHIV'),
              ('supp_share', 'Suppressed, share of those on ART'), ('mtct_n', 'Infections of children a year'), ('circ_15_49', 'Circumcised, men 15-49')]
    fig, axs = plt.subplots(2, 5, figsize=(17, 6))
    for ax, (col, title) in zip(axs.flat, panels):
        for m, lk in looks.items():
            g = s[s.model == m].groupby('year')[col]
            ax.fill_between(g.median().index, g.quantile(0.05), g.quantile(0.95), color=lk['color'], alpha=0.15, lw=0)
            ax.plot(g.median().index, g.median(), lw=1.5, **lk)
        if col in ref:
            g = ref.groupby('year')[col].median()
            ax.plot(g.index, g, color=c_muted, lw=1, ls='--', label="EMOD′, condensed cascade (analysis 06)")
        style(ax, title)
    axs[0, 0].legend(fontsize=7, frameon=False)
    axs[1, 0].legend(fontsize=7, frameon=False)
    save(fig, 'cascade_over_time.png')
    cols = [c for c, _ in panels]
    s[s.year.isin([2005, 2010, 2015, 2020])].groupby(['model', 'year'])[cols].median().round(3).to_csv(out / 'series_medians.csv')
    ref[ref.year.isin([2005, 2010, 2015, 2020])].groupby('year')[['dx_share', 'art_cov', 'off_share', 'supp_share']].median().round(3).to_csv(out / 'emod_condensed_shares.csv')

    # Sanity check of the host
    c = s[s.model == 'condensed']
    ev = c.groupby('year')[['conceptions', 'ev_twelve_weeks', 'births', 'ev_six_weeks', 'ev_debut', 'ev_symptoms', 'art_starts', 'art_stops', 'diagnoses']].median()
    ev.loc[[2000, 2005, 2010, 2015, 2020]].round(1).to_csv(out / 'check_events.csv')
    bad = [x for x in c.columns if x.startswith('bad_')]
    c[bad].agg(['sum', 'max']).to_csv(out / 'check_states.csv')
    care = c.groupby('year')[[f'care_{k}' for k in CARE]].median().div(c.groupby('year').plhiv.median(), axis=0)
    care.loc[[2005, 2010, 2015, 2020]].round(3).to_csv(out / 'check_care.csv')
    m = s[s.model.isin(['condensed', 'condensed_nopmtct'])]
    pm = m[m.year.between(2001, 2021)].groupby(['model', 'seed'])[['mtct', 'pmtct_mothers']].sum().groupby('model').mean()
    pm.round(1).to_csv(out / 'check_pmtct.csv')
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.6))
    for ax, (a, b, title) in zip(axs, [('conceptions', 'ev_twelve_weeks', 'Conceptions, and 12-week pregnancies a year'), ('births', 'ev_six_weeks', 'Births, and 6-week-olds a year'),
                                         ('art_starts', 'art_stops', 'ART starts and stops a year')]):
        for col, ls in ((a, '-'), (b, '--')):
            g = c.groupby('year')[col].median()
            ax.plot(g.index, g, ls=ls, color=c_run, label=col.replace('_', ' '))
        style(ax, title, 'Agents')
        ax.legend(fontsize=7, frameon=False)
    save(fig, 'check_events.png', 'Condensed cascade, median over 50 seeds, agents (unscaled).')
    fig, ax = plt.subplots(figsize=(5, 3.2))
    for mdl, lk in (('condensed', dict(color=c_run, label='With PMTCT')), ('condensed_nopmtct', dict(color=c_muted, label='PMTCT efficacy 0'))):
        g = s[s.model == mdl].groupby('year').mtct_n.mean()
        ax.plot(g.index, g, lw=1.5, **lk)
    style(ax, 'Infections of children a year (mean of runs)')
    ax.legend(fontsize=7, frameon=False)
    save(fig, 'check_pmtct.png', 'Condensed cascade: 50 seeds with PMTCT, 20 without.')

    t = targetfit.load_targets()
    vals = {m: targetfit.load_runs(here, m)[1] for m in looks}
    rows = []
    for tid in t.index:
        if tid in vals['own'] and tid in vals['condensed']:
            a, b = vals['own'][tid].dropna(), vals['condensed'][tid].dropna()
            if len(a) and len(b):
                rows.append(dict(tid=tid, target=t.loc[tid, 'target'], sex=t.loc[tid, 'sex'], year=t.loc[tid, 'year'], source=t.loc[tid, 'source'], data=t.loc[tid, 'value'],
                                 own=a.median(), condensed=b.median()))
    cmp = pd.DataFrame(rows).assign(ratio=lambda d: d.condensed/d.own)
    cmp.to_csv(out / 'targets_compared.csv', index=False)
    sc_ = {m: np.array([targetfit.score(dict(r.dropna()), t) for _, r in vals[m].iterrows()]) for m in looks}
    pd.DataFrame({m: dict(median_score=np.median(v), best_score=v.min()) for m, v in sc_.items()}).round(2).to_csv(out / 'scores.csv')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=30)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
