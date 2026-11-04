"""
Analysis 08: the condensed cascade as rates in Optima′. See README.md.

    PYTHONPATH=~/starsim python run.py emod    # (i) EMOD′ with the condensed cascade, 3 seeds: care-state moves seen and expected
    python run.py optima                      # (ii) Optima′ with its own cascade and with EMOD's (Starsim 4.0.0rc1)
    python run.py --plot                      # figures/
"""
import sys
import json
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
out = here / 'outputs'
TABLES = [repo / 'common' / 'cascade' / 'tables' / f'malawi_emod_set000_{k}.csv' for k in ('params', 'flows')]
STATES = ['undx', 'dx', 'care', 'lost', 'svl', 'usvl']
BANDS = ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50']
GROUPS = [dict(sex=s, age_lo=5*b, age_hi=5*b + 5 if b < 17 else 100) for b in range(18) for s in 'fm'] # EMOD's people by sex and 5-year band (Optima's general populations)
START, STOP, DT = 1990, 2025, 1/12
SEEDS = [1, 2, 3]
VARIANTS = dict(own=dict(), emodcx=dict(cascade=True), emodcx_noflows=dict(cascade=True, flows=False))


def params():
    sys.path.insert(0, str(repo))
    from common.cascade.params import Params
    return Params(*TABLES)


# %% (i) EMOD′

def emod_one(seed):
    for p in ('models/emod/v4', ''):
        sys.path.insert(0, str(repo / p))
    import starsim as ss
    import harness, inputs
    from common.cascade.cascade import STAGING, PREART, LINKING, OUT, LOST, TESTING

    class CareFlows(ss.Analyzer):
        """ Each month from 1990: person-time and moves between care states of people with HIV, by CD4 band and group; entries into CD4 < 200; pregnancies; hidden shares """
        def init_pre(self, sim):
            super().init_pre(sim)
            T, S, B, G = int(round((STOP - START)/DT)), len(STATES), len(BANDS), len(GROUPS)
            self.pt, self.moves = np.zeros((T, S, B, G)), np.zeros((T, S, S, B, G))
            self.entries, self.preg, self.women = np.zeros((T, S, G)), np.zeros((T, G)), np.zeros((T, G))
            self.shares = np.full((T, 4), np.nan)
            self.prev = {}

        def step(self):
            sim = self.sim
            t = int(round((self.t.now('year') - START)/DT))
            if not 0 <= t < len(self.pt):
                return
            ppl, hiv, cx = sim.people, sim.diseases.hiv, sim.interventions.cascade
            au = ppl.auids
            age = ppl.emod_age.raw[au]/365
            g = 2*np.minimum(age//5, 17).astype(int) + ~ppl.female.raw[au]
            women = ppl.female.raw[au]
            np.add.at(self.women[t], g[women], DT)
            np.add.at(self.preg[t], g[np.isin(au, sim.demographics.emoddemog.twelve_weeks)], 1)
            care = cx.care.raw[au]
            inf = hiv.infected.raw[au]
            on = hiv.art.raw[au] != 0
            diag = hiv.diagnosed.raw[au].astype(int)
            code = np.select([on & (cx.arm.raw[au] == 1), on, care == STAGING, np.isin(care, [PREART, LINKING]), np.isin(care, [OUT, LOST])], [4, 5, -1, 2, 3], diag)
            cd4 = np.where(inf, hiv.cd4(au), np.inf)
            band = np.where(hiv.stage.raw[au] == 0, 0, 5 - np.searchsorted([50, 200, 350, 500], cd4, side='right'))
            u, c, b, gg = au[inf], code[inf], band[inf], g[inf]
            pc, pb, pg = self.lookup(u)
            c = np.where(c == -1, np.where(pc >= 0, pc, diag[inf]), c) # Staging lasts a day (to the next step): people in it are where they came from
            code[inf] = c
            np.add.at(self.pt[t], (c, b, gg), DT)
            seen = pc >= 0
            mv = seen & (pc != c)
            np.add.at(self.moves[t], (pc[mv], c[mv], pb[mv], pg[mv]), 1)
            ent = seen & (pc <= 1) & (pb <= 3) & (b >= 4)
            np.add.at(self.entries[t], (pc[ent], pg[ent]), 1)
            self.remember(u, c, b, gg)
            debuted = age >= sim.networks.emodnet.debut.raw[au]
            loop = care == TESTING
            sel = [debuted & ~inf, debuted & inf & (code == 0), debuted & inf & (code == 1), inf & (code == 3)]
            self.shares[t] = [loop[s].mean() if s.any() else np.nan for s in sel[:3]] + [(care[sel[3]] == OUT).mean() if sel[3].any() else np.nan]

        def lookup(self, u):
            n = len(self.prev.get('c', []))
            ok = u < n
            res = [np.full(len(u), -1) for _ in range(3)]
            for r, k in zip(res, 'cbg'):
                r[ok] = self.prev[k][u[ok]] if n else -1
            return res

        def remember(self, u, c, b, g):
            n = int(u.max()) + 1 if len(u) else 0
            for k, v in zip('cbg', (c, b, g)):
                arr = np.full(n, -1)
                arr[u] = v
                self.prev[k] = arr

    inp = inputs.load(0)
    sim = harness.make_sim(0, seed=inp.run_number*10 + seed, stop=STOP, cascade=True, condensed=True, record=False, analyzers=CareFlows(name='careflows'))
    sim.run()
    a = sim.analyzers.careflows
    return {k: getattr(a, k) for k in ('pt', 'moves', 'entries', 'preg', 'women', 'shares')}


def emod(cores):
    from common.cascade.rates import CompartmentalCascade
    with Pool(cores) as pool:
        res = pool.map(emod_one, SEEDS)
    tot = {k: sum(r[k] for r in res) for k in ('pt', 'moves', 'entries', 'preg', 'women')}
    years = START + DT*np.arange(len(tot['pt']))
    yr = np.floor(years + 1e-9).astype(int)
    preg = pd.DataFrame(tot['preg']).groupby(yr).sum()/pd.DataFrame(tot['women']).groupby(yr).sum().replace(0, np.nan)
    cx = CompartmentalCascade(params(), GROUPS, years, pregnancy=preg.fillna(0).loc[yr].values.T)
    rows = []
    for (src, dst), R in cx.rates.items():
        i, j = STATES.index(src), STATES.index(dst)
        exp = np.einsum('tbg,bgt->t', tot['pt'][:, i], R)
        if (src, dst) in cx.on_progression:
            exp = exp + np.einsum('tg,gt->t', tot['entries'][:, i], cx.on_progression[(src, dst)])
        rows.append(pd.DataFrame(dict(year=yr, src=src, dst=dst, person_years=tot['pt'][:, i].sum(axis=(1, 2)), seen=tot['moves'][:, i, j].sum(axis=(1, 2)), expected=exp)))
    for i, src in enumerate(STATES): # Moves seen that rates.py has no rate for
        for j, dst in enumerate(STATES):
            if i != j and (src, dst) not in cx.rates:
                rows.append(pd.DataFrame(dict(year=yr, src=src, dst=dst, person_years=tot['pt'][:, i].sum(axis=(1, 2)), seen=tot['moves'][:, i, j].sum(axis=(1, 2)), expected=0.0)))
    df = pd.concat(rows).groupby(['year', 'src', 'dst']).agg(person_years=('person_years', 'sum'), seen=('seen', 'sum'), expected=('expected', 'sum')).reset_index()
    df = df.assign(rate_seen=df.seen/df.person_years, rate_expected=df.expected/df.person_years)
    out.mkdir(exist_ok=True)
    df.round(6).to_csv(out / 'emod_rates.csv', index=False)
    sh = np.nanmean(np.array([r['shares'] for r in res]), axis=0)
    hidden = pd.DataFrame(dict(year=yr, seen_everyone=sh[:, 0], seen_undx=sh[:, 1], seen_dx=sh[:, 2], seen_out=sh[:, 3],
                               rates_everyone=cx.loop_share['everyone'], rates_undx=cx.loop_share['undx'], rates_dx=cx.loop_share['dx'], rates_out=cx.out_share)).groupby('year').mean()
    hidden.round(4).to_csv(out / 'emod_hidden_shares.csv')
    preg.round(4).to_csv(out / 'emod_pregnancy_rates.csv')
    return


# %% (ii) Optima′

def optima_one(variant):
    sys.path.insert(0, str(repo / 'models' / 'optima' / 'v4'))
    sys.path.insert(0, str(repo / 'models' / 'optima' / 'reference'))
    sys.path.insert(0, str(repo / 'common'))
    import optima_prime as opp, compare as cmp
    import optima_targets, targetfit
    ref = opp.load()
    kw = dict(VARIANTS[variant])
    if kw.pop('cascade', False):
        kw['cascade'] = params()
    sim = opp.make_sim(ref, components=cmp.COMPONENTS['E'], dt=0.02, end=STOP + 1, **kw).run()
    r = sim.results.optima
    tvec = ref.simpars['tvec']
    tv = tvec[tvec <= STOP + 0.2]
    names = [s.replace('-', '_') for s in ref.settings['statelabels']]
    t = np.array([d.years for d in r[names[0]].by('pop').timevec])
    at = [np.argmin(abs(t - y)) for y in tv]
    people = np.array([np.asarray(r[n].by('pop'))[at].T for n in names]) # States × pops × times
    window = lambda arr: np.array([arr[(t >= y - 1e-9) & (t < y + 0.2 - 1e-9)].mean(axis=0) for y in tv]).T/np.diff(t).mean() # Per year, over each 0.2 years (results may be kept less often than each step)
    inf = [k for k in r.keys() if k.startswith('infection_') or (k.startswith('birth_') and k != 'birth')]
    inci = window(sum(np.asarray(r[k].by('pop')) for k in inf))
    rr = np.asarray(ref.simpars['rrcomorbiditydeathtx'])
    deaths = 0
    for k in [k for k in r.keys() if k.startswith('hivdeath_')]:
        d = window(np.asarray(r[k].by('pop')))
        deaths = deaths + (d*rr[:, :len(tv)] if k.split('_')[1] in opp.TREATED else d)
    raw = dict(people=people, inci=inci, death=deaths[None], tvec=tv, popkeys=ref.raw['popkeys'])
    values = optima_targets.measure(dict(raw=raw, settings=ref.settings), targetfit.load_targets())
    s = ref.settings
    tot = lambda idx: people[np.atleast_1d(idx)].sum(axis=(0, 1))
    plhiv = tot(s['allplhiv'])
    series = pd.DataFrame(dict(year=tv, plhiv=plhiv, new_infections=inci.sum(0), hiv_deaths=deaths.sum(0), on_art=tot(s['alltx']),
                               **{f'share_{c}': tot(s[c])/plhiv for c in ['undx', 'dx', 'care', 'lost', 'usvl', 'svl']},
                               circumcised=people[1].sum(0)), index=None).assign(model=variant)
    return dict(variant=variant, values=values, series=series)


def optima(cores):
    with Pool(cores) as pool:
        res = pool.map(optima_one, list(VARIANTS))
    out.mkdir(exist_ok=True)
    for r in res:
        (out / f'optima_{r["variant"]}_results.jsonl').write_text(json.dumps(dict(i=0, kind='sample', pars=dict(variant=r['variant']), values=r['values'])) + '\n')
    pd.concat([r['series'] for r in res]).round(5).to_csv(out / 'optima_series.csv', index=False)
    return


# %% Figures and tables

def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(repo / 'common'))
    import targetfit
    from targetfit import style, c_run, c_best, c_muted, c_data
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)

    # (i) Rates seen in EMOD′ and given by rates.py
    df = pd.read_csv(out / 'emod_rates.csv')
    df = df[df.year < STOP]
    ex = df.groupby(['year', 'src']).agg(person_years=('person_years', 'first'), seen=('seen', 'sum'), expected=('expected', 'sum')).reset_index()
    ex = ex.assign(rate_seen=ex.seen/ex.person_years, rate_expected=ex.expected/ex.person_years)
    ex.round(5).to_csv(out / 'emod_exit_rates.csv', index=False)
    pairs = [('undx', None), ('undx', 'dx'), ('undx', 'svl'), ('undx', 'care'), ('dx', None), ('dx', 'svl'), ('care', 'dx'), ('care', 'svl'),
             ('svl', 'lost'), ('svl', 'dx'), ('lost', 'dx'), ('usvl', 'lost')]
    fig, axs = plt.subplots(3, 4, figsize=(14, 8.5), constrained_layout=True)
    for ax, (src, dst) in zip(axs.flat, pairs):
        d = ex[ex.src == src] if dst is None else df[(df.src == src) & (df.dst == dst)]
        d = d[d.person_years > 20]
        ax.plot(d.year, d.rate_seen, 'o-', ms=2.5, lw=1, color=c_best, label='Seen in EMOD′')
        ax.plot(d.year, d.rate_expected, lw=1.5, color=c_run, label='rates.py, same people')
        style(ax, f'{src} → {dst or "any"}', 'Per person-year')
    axs.flat[0].legend(fontsize=7, frameon=False)
    fig.text(0.01, -0.02, 'EMOD′ with the condensed cascade, set 0, 3 seeds: moves between care states of people with HIV each month, over person-time. rates.py: its rates for each month, CD4 band and group, on the same person-time.', fontsize=7, color=c_muted)
    fig.savefig(figdir / 'emod_rates.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    h = pd.read_csv(out / 'emod_hidden_shares.csv', index_col=0)
    h = h[h.index < STOP]
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.2), constrained_layout=True)
    for ax, (k, title) in zip(axs, [('everyone', 'In the testing loop: HIV-negative, past debut'), ('undx', '... undiagnosed'), ('dx', '... diagnosed, not in care'), ('out', 'Out (not lost for good), of lost')]):
        ax.plot(h.index, h[f'seen_{k}'], 'o-', ms=2.5, lw=1, color=c_best, label='Seen in EMOD′')
        ax.plot(h.index, h[f'rates_{k}'], lw=1.5, color=c_run, label='rates.py')
        style(ax, title)
    axs[0].legend(fontsize=7, frameon=False)
    fig.savefig(figdir / 'emod_hidden_shares.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # (ii) Optima′ with its own cascade and with EMOD's
    s = pd.read_csv(out / 'optima_series.csv')
    looks = dict(own=dict(color=c_best, label="Optima′, own cascade"), emodcx=dict(color=c_run, label="Optima′, EMOD's cascade"), emodcx_noflows=dict(color=c_run, ls='--', label="... without VMMC and PrEP"))
    t = targetfit.load_targets()
    _, e06, _ = targetfit.load_runs(repo / 'analyses' / '06_condensed_cascade', 'condensed')
    panels = [('plhiv', 'plhiv', 'PLHIV'), ('new_infections', 'new_infections', 'New infections a year'), ('hiv_deaths', 'aids_deaths', 'HIV deaths a year'), ('on_art', 'on_art', 'On ART')]
    shares = [('share_undx', 'Undiagnosed'), ('share_dx', 'Diagnosed, not in care'), ('share_care', 'In care, not on ART'), ('share_lost', 'Lost'), ('share_svl', 'On ART, suppressed'), ('share_usvl', 'On ART, not suppressed')]
    fig, axs = plt.subplots(2, 5, figsize=(17, 6.5), constrained_layout=True)
    for ax, (col, target, title) in zip(axs.flat[:4], panels):
        for m, lk in looks.items():
            d = s[s.model == m]
            ax.plot(d.year, d[col], lw=1.5, **lk)
        d = t[(t.target == target) & (t.sex == 'all') & (t.group == 'all') & (t.age_lo == 0)] if target != 'on_art' else t.iloc[0:0]
        ax.plot(d.year, d.value, 'o', ms=3, color=c_data, label='Data (UNAIDS)')
        style(ax, title)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:,.0f}k'))
    for ax, (col, title) in zip(list(axs.flat[4:]), shares):
        for m, lk in looks.items():
            d = s[s.model == m]
            ax.plot(d.year, d[col], lw=1.5, **lk)
        style(ax, f'{title}, share of PLHIV')
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    axs.flat[0].legend(fontsize=7, frameon=False)
    fig.savefig(figdir / 'optima_cascade.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Targets: each variant and EMOD′ (analysis 06, condensed: median of 30 runs), by kind
    vals = {m: targetfit.load_runs(here, f'optima_{m}')[1].iloc[0] for m in VARIANTS}
    rows = []
    for tid, r in t.iterrows():
        if tid not in vals['own']:
            continue
        rows.append(dict(tid=tid, target=r.target, year=r.year, sex=r.sex, age_lo=r.age_lo, age_hi=r.age_hi, group=r.group, data=r.value,
                         **{m: vals[m].get(tid, np.nan) for m in VARIANTS}, emod06=e06[tid].median() if tid in e06 else np.nan))
    tab = pd.DataFrame(rows)
    tab.round(5).to_csv(out / 'targets_compared.csv', index=False)
    kinds = tab.assign(**{f'{m}_ratio': tab[m]/tab.own for m in ['emodcx', 'emodcx_noflows', 'emod06']}).groupby('target')[[f'{m}_ratio' for m in ['emodcx', 'emodcx_noflows', 'emod06']]].median()
    kinds.loc[kinds.index.isin(['plhiv', 'on_art', 'population', 'new_infections', 'aids_deaths']), 'emod06_ratio'] = np.nan # EMOD′'s population is not Optima's: compare its shares only
    scores = {m: targetfit.score(dict(vals[m].dropna()), t) for m in VARIANTS}
    kinds.round(3).to_csv(out / 'targets_by_kind.csv')
    pd.Series(scores, name='score').round(3).to_csv(out / 'scores.csv')
    print(kinds.round(3), '\n', scores)
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('part', nargs='?', choices=['emod', 'optima'])
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=3)
    a = ap.parse_args()
    sys.path.insert(0, str(repo))
    if a.plot:
        plot()
    elif a.part == 'emod':
        emod(a.cores)
    else:
        optima(a.cores)
