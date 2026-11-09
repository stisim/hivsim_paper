"""
13: HIVsim's mortality on ART (see README.md)

    PYTHONPATH=/tmp/sti-g:/tmp/ss-nov09 conda run -n hp_rc3 python run.py hivsim   # plain hivsim.Sim, high ART coverage
    ... python run.py malawi   # Malawi exp 03 best run, and variants
    ... python run.py optima   # Optima's death rates by cascade state
    ... python run.py emod     # EMOD's ARTMortalityTable
    ... python run.py plot
"""
import sys
import json
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

here = Path(__file__).resolve().parent
repo = here.parents[1]
out, figs = here / 'outputs', here / 'figures'
malawi_repo = Path.home() / 'hivsim_malawi'

CD4_BINS = [0, 50, 200, 350, 500, np.inf]
CD4_LABELS = ['<50', '50-200', '200-350', '350-500', '500+']
DUR_BINS = [0, 0.5, 1, 2, 5, np.inf] # Years since first ART start
DUR_LABELS = ['0-6m', '6-12m', '1-2y', '2-5y', '5y+']
STATUS = {1: 'suppressed', 2: 'nonsuppressed', 3: 'off_after_art'}
CAUSES = ['art_death', 'post_art_death', 'hiv_death', 'aids_death', 'background']


class ARTCohort(ss.Analyzer):
    """
    Follows everyone from their first ART start: person-time and deaths each step by status at the start of the step
    (suppressed, non-suppressed, off ART after starting), CD4 at first start, and time since first start; deaths by cause.
    Also HIV deaths by year among the never treated, and mean CD4 on ART by time since first start.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.rows = {} # (year, status, cd4bin, durbin) -> [person-years, deaths by cause...]
        self.cd4_rows = {} # (cd4bin, durbin) -> [sum CD4, n]
        self.never = {} # year -> HIV deaths among never treated
        return

    def init_post(self):
        super().init_post()
        self.prev = None
        return

    def grow(self, n):
        for k in ['ti0', 'cd4_0']:
            a = getattr(self, k, np.full(0, np.nan))
            if len(a) < n:
                setattr(self, k, np.concatenate([a, np.full(n - len(a), np.nan)]))
        return

    def snapshot(self):
        sim = self.sim
        art, hiv = sim.interventions.art, sim.diseases.hiv
        n = len(sim.people.alive.raw)
        self.grow(n)
        status = np.zeros(n, int)
        status[art.on_effective_art.raw[:n]] = 1
        status[art.on_nonsuppressive_art.raw[:n]] = 2
        status[art.art_discontinued.raw[:n]] = 3
        alive = sim.people.alive.raw[:n].copy()
        status[~alive] = 0
        new = (status > 0) & np.isnan(self.ti0[:n])
        self.ti0[new] = self.ti
        self.cd4_0[new] = art.cd4_preart.raw[new] # Set at the start, from the CD4 count then
        return dict(n=n, alive=alive, status=status, infected=hiv.infected.raw[:n].copy())

    def step(self):
        sim = self.sim
        art, hiv = sim.interventions.art, sim.diseases.hiv
        year = int(np.floor(sim.t.now('year')))
        dt = float(sim.t.dt.years) if hasattr(sim.t.dt, 'years') else 1/12
        if self.prev is not None:
            p = self.prev
            n = p['n']
            died = p['alive'] & ~sim.people.alive.raw[:n]
            cause = np.full(n, 'background', dtype=object)
            for name, tr in [('art_death', art.art_death), ('post_art_death', art.post_art_death), ('hiv_death', hiv.hiv_death), ('aids_death', hiv.aids_death)]:
                u = np.asarray(tr.uids, dtype=int)
                u = u[u < n]
                cause[u] = name
            # Never treated: HIV deaths
            nt = died & (p['status'] == 0) & p['infected'] & (cause != 'background')
            self.never[year] = self.never.get(year, 0) + int(nt.sum())
            # The cohort
            idx = np.flatnonzero(p['status'] > 0)
            if len(idx):
                dur = (self.ti - 1 - self.ti0[idx]) * dt
                cb = np.digitize(self.cd4_0[idx], CD4_BINS[1:-1])
                db = np.digitize(dur, DUR_BINS[1:-1])
                st = p['status'][idx]
                d = died[idx]
                c = cause[idx]
                keys = st*100 + cb*10 + db
                for k in np.unique(keys):
                    m = keys == k
                    row = self.rows.setdefault((year, int(k//100), int(k//10 % 10), int(k % 10)), np.zeros(1 + len(CAUSES)))
                    row[0] += m.sum()*dt
                    for j, cname in enumerate(CAUSES):
                        row[1 + j] += np.sum(d[m] & (c[m] == cname))
                # Mean CD4 on ART (status 1-2, alive now)
                on = (st <= 2) & ~d
                cd4 = hiv.cd4.raw[idx[on]]
                for k in np.unique(cb[on]*10 + db[on]):
                    m = (cb[on]*10 + db[on]) == k
                    r = self.cd4_rows.setdefault((int(k//10), int(k % 10)), np.zeros(2))
                    r += [cd4[m].sum(), m.sum()]
        self.prev = self.snapshot()
        return

    def to_df(self, scale=1.0):
        rows = [dict(year=y, status=STATUS[s], cd4_0=CD4_LABELS[c], dur=DUR_LABELS[d], py=v[0]*scale, **{k: v[1+j]*scale for j, k in enumerate(CAUSES)})
                for (y, s, c, d), v in self.rows.items()]
        cohort = pd.DataFrame(rows)
        cd4 = pd.DataFrame([dict(cd4_0=CD4_LABELS[c], dur=DUR_LABELS[d], cd4_sum=v[0], n=v[1]) for (c, d), v in self.cd4_rows.items()])
        never = pd.Series(self.never, name='never_treated_hiv_deaths')*scale
        return cohort, cd4, never


# %% Plain hivsim.Sim

def coverage():
    return {'year': [2000, 2004, 2010, 2015, 2020, 2030], 'value': [0, 0, 0.3, 0.6, 0.8, 0.85]}


def make_plain(seed, n_agents=30e3, art_pars=None):
    import hivsim
    tests = [sti.HIVTest(name='general', test_prob_data=0.2, start=2004),
             sti.HIVTest(name='low_cd4', test_prob_data=0.8, start=2004, eligibility=lambda sim: ~sim.diseases.hiv.diagnosed & (sim.diseases.hiv.cd4 < 200))]
    intvs = tests + [sti.ART(coverage=coverage(), **(art_pars or {}))]
    return hivsim.Sim(n_agents=n_agents, start=2000, stop=2030, rand_seed=seed, verbose=0, interventions=intvs, analyzers=ARTCohort())


# %% Malawi

def best_pars():
    for line in open(malawi_repo / 'experiments/03_testing_art/outputs/results.jsonl'):
        d = json.loads(line)
        if d['kind'] == 'sample' and d['i'] == json.load(open(malawi_repo / 'experiments/03_testing_art/outputs/coverage_summary.json'))['best_run']:
            return d['seed'], d['pars']


BETA = 0.004


def make_malawi(seed, n_agents=50e3, art_pars=None, art_class=None):
    sys.path.insert(0, str(malawi_repo / 'hivsim'))
    ss.MaternalNet = ss.PrenatalNet # malawi.py was written for Starsim 4.0.0rc1, whose MaternalNet is now PrenatalNet
    import malawi
    _, pars = best_pars()
    pars = pars | dict(hiv_beta_m2f=BETA) # Under PR #612 the exp 03 best run's beta (0.008) gives 2.3M PLHIV in 2021; 0.004 gives ~1.1M, as UNAIDS and as that run had
    hiv_pars = {k[4:]: v for k, v in pars.items() if k.startswith('hiv_')}
    nw_pars = {k[3:]: v for k, v in pars.items() if k.startswith('nw_')}
    other = {k: v for k, v in pars.items() if not k.startswith(('hiv_', 'nw_'))}
    sim = malawi.make_sim(seed=seed, n_agents=n_agents, stop=2022, verbose=0, hiv_pars=hiv_pars, nw_pars=nw_pars, analyzers=ARTCohort(), **other)
    if art_pars or art_class:
        intvs = list(sim.pars.interventions)
        i = next(i for i, m in enumerate(intvs) if isinstance(m, sti.ART))
        old = intvs[i]
        cls = art_class or sti.ART
        intvs[i] = cls(coverage=old._raw_coverage, vls_coverage=old._raw_vls_coverage, **(art_pars or {}))
        sim.pars.interventions = intvs
    return sim


def run_one(job):
    kind, variant, seed, kw = job
    sim = (make_plain if kind == 'hivsim' else make_malawi)(seed, **kw)
    sim.run()
    scale = float(sim.pars.pop_scale) if sim.pars.get('pop_scale') else 1.0
    cohort, cd4, never = sim.analyzers[0].to_df(scale=scale)
    res = sim.results
    years = np.floor(np.asarray(sim.t.yearvec)).astype(int)
    annual = pd.DataFrame(dict(year=years, hiv_deaths=np.asarray(res.hiv.new_deaths), on_art=np.asarray(res.hiv.on_art), plhiv=np.asarray(res.hiv.infected)))
    annual = annual.groupby('year').agg(hiv_deaths=('hiv_deaths', 'sum'), on_art=('on_art', 'mean'), plhiv=('plhiv', 'mean'))
    annual = annual.join(never).fillna(0)
    for df in (cohort, cd4, annual):
        df['variant'], df['seed'] = variant, seed
    return cohort, cd4, annual.reset_index()


def run(kind, variants, seeds):
    jobs = [(kind, v, s, kw) for v, kw in variants.items() for s in seeds]
    with Pool(min(len(jobs), 60)) as pool:
        results = pool.map(run_one, jobs)
    for i, name in enumerate(['cohort', 'cd4', 'annual']):
        pd.concat([r[i] for r in results]).to_csv(out / f'{kind}_{name}.csv', index=False)
    return


# %% Variants of ART for the Malawi runs (analysis only: STIsim is unchanged)

def emod_table_arrays():
    """ EMOD's four ARTMortalityTables as arrays: {(sex, arm): (dur_bins_years, age_bins, cd4_bins, mort[dur, age, cd4])}; the male
    non-suppressed 6-12 month row (0.56 at CD4 75, a known error, docs/emod_inventory.md) is replaced by the female row """
    camp = json.loads(sc.runcommand(f'git -C {repo} show origin/starsim-rc3:models/emod/reference/source/InputFiles/Templates/campaign_Malawi_20251215.json'))
    tabs = {}
    def find(o, ctx):
        if isinstance(o, dict):
            ctx = ctx | {k: o[k] for k in ['Target_Gender'] if k in o}
            if o.get('class') == 'ARTMortalityTable':
                tabs[(ctx['Target_Gender'], o['Intervention_Name'])] = (np.array(o['ART_Duration_Days_Bins'], float)/365, np.array(o['Age_Years_Bins'], float),
                                                                       np.array(o['CD4_Count_Bins'], float), np.array(o['MortalityTable'], float))
                return
            for v in o.values():
                find(v, ctx)
        elif isinstance(o, list):
            for v in o:
                find(v, ctx)
    find(camp, {})
    tabs[('Male', 'NotSuppressedOnART')][3][1] = tabs[('Female', 'NotSuppressedOnART')][3][1]
    return tabs


class ARTByStartCD4(sti.ART):
    """
    Test variant: death on ART at EMOD's (IeDEA-derived) rates by CD4 at the latest ART start, time since it, age, sex and
    suppression, as a piecewise-constant hazard (not EMOD's redraw from ART start), in place of the off-ART rate at the
    current CD4 count × rel_art_mortality
    """
    def init_pre(self, sim):
        super().init_pre(sim)
        self.tabs = emod_table_arrays()
        return

    def on_art_death_rate(self, uids):
        hiv = self.hiv
        ppl = self.sim.people
        rate = np.zeros(len(uids))
        dur = (hiv.ti - self.ti_art[uids]) * float(hiv.dt.years)
        age = ppl.age[uids]
        cd4 = self.cd4_preart[uids]
        male = ppl.male[uids]
        eff = self.on_effective_art[uids]
        for (sex, arm), (db, ab, cb, mort) in self.tabs.items():
            m = (male == (sex == 'Male')) & (eff == (arm == 'SuppressedOnART'))
            if not m.any():
                continue
            di = np.minimum(np.searchsorted(db, dur[m], side='right'), len(db) - 1)
            ai = np.minimum(np.searchsorted(ab, age[m], side='right'), len(ab) - 1)
            rows = mort[di, ai] # n × CD4 bins
            rate[m] = [np.interp(c, cb, r) for c, r in zip(cd4[m], rows)]
        return rate


class ARTNoRecoveryUnsupp(sti.ART):
    """ Test variant: no CD4 recovery while on non-suppressive ART (the CD4 count stays where it was at the start, or falls no further) """
    def update_cd4(self):
        hiv = self.hiv
        before = hiv.cd4[self.on_nonsuppressive_art.uids].copy()
        super().update_cd4()
        u = self.on_nonsuppressive_art.uids
        hiv.cd4[u] = np.minimum(before, hiv.cd4[u])
        return


ELIGIBILITY = [(2011, 200), (2014, 350), (2016, 500)] # Malawi's CD4 thresholds for starting ART (first start only): <200 to 2010, <350 2011-13, <500 2014-15, all from 2016


class EligibilityMixin:
    """ Test variant: a first ART start only below Malawi's CD4 threshold of the year (or when pregnant, from Option B+ in 2011); restarts are not limited """
    def ready(self, uids=None):
        ok = super().ready(uids)
        if uids is None:
            return ok
        year = self.sim.t.now('year')
        limit = next((cd4 for y, cd4 in ELIGIBILITY if year < y), np.inf)
        preg = self.sim.demographics.pregnancy.P[uids] if year >= 2011 and hasattr(self.sim.demographics, 'pregnancy') else np.zeros(len(uids), bool)
        return ok & ((self.hiv.cd4[uids] < limit) | ~self.art_naive[uids] | preg)


class ARTEligibility(EligibilityMixin, sti.ART):
    pass


class ARTByStartCD4Eligibility(EligibilityMixin, ARTByStartCD4):
    pass


MALAWI_VARIANTS = {
    'base': {},
    'rel_eff_0.75': dict(art_pars=dict(rel_art_mortality_effective=0.75)), # The largest allowed by ART.init_pre() with the age multiplier 1.32
    'emod_table': dict(art_class=ARTByStartCD4),
    'unsupp_no_recovery': dict(art_class=ARTNoRecoveryUnsupp),
    'eligibility': dict(art_class=ARTEligibility),
    'emod_table_eligibility': dict(art_class=ARTByStartCD4Eligibility),
    'dropout_0.1': dict(art_pars=dict(rate_art_dropout=ss.peryear(0.1))), # Disengagement (Optima's leavecare is ~0.09/yr); a coverage target refills the places
    'emod_table_dropout_0.1': dict(art_class=ARTByStartCD4, art_pars=dict(rate_art_dropout=ss.peryear(0.1))),
}


# %% Comparators

def optima():
    """ Optima's calibrated Malawi run: HIV deaths per person-year by cascade state, and by CD4 band on suppressive ART """
    ref = sc.load(repo / 'models/optima/reference/optima_malawi_ref.obj')
    st, raw = ref['settings'], ref['raw']
    t = np.asarray(raw['tvec'])
    labels = list(st['statelabels'])
    people, death = raw['people'].sum(axis=1), raw['death'].sum(axis=1) # States × time; deaths are per year
    groups = dict(untreated=np.concatenate([st['undx'], st['dx'], st['care']]), lost=st['lost'], nonsuppressed=st['usvl'], suppressed=st['svl'])
    for band in ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50']:
        groups[f'suppressed_{band}'] = np.array([labels.index(f'svl-{band}')])
        groups[f'nonsuppressed_{band}'] = np.array([labels.index(f'usvl-{band}')])
    rows = []
    for year in range(2005, 2031):
        k = np.flatnonzero(np.floor(t) == year) # Steps of 0.2 years in the year
        for g, idx in groups.items():
            n = people[idx][:, k].sum(axis=0).mean()
            d = death[idx][:, k].sum(axis=0).mean()
            rows.append(dict(year=year, group=g, people=n, hiv_deaths=d, rate=d/n if n > 0 else np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'optima_rates.csv', index=False)
    return df


def emod_tables():
    """ EMOD's ARTMortalityTable (Malawi campaign) at CD4 counts at ART start that represent this analysis's bins, interpolated as EMOD does """
    camp = json.loads(sc.runcommand(f'git -C {repo} show origin/starsim-rc3:models/emod/reference/source/InputFiles/Templates/campaign_Malawi_20251215.json'))
    found = []
    def find(o, ctx):
        if isinstance(o, dict):
            ctx = ctx | {k: o[k] for k in ['Target_Gender'] if k in o}
            if o.get('class') == 'ARTMortalityTable':
                found.append((ctx['Target_Gender'], o['Intervention_Name'], o))
                return
            for v in o.values():
                find(v, ctx)
        elif isinstance(o, list):
            for v in o:
                find(v, ctx)
    find(camp, {})
    rows = []
    cd4_rep = dict(zip(CD4_LABELS, [25, 125, 275, 425, 600]))
    durs = ['0-6m', '6-12m', '1-2y', '2-3y', '3y+']
    for sex, arm, tab in found:
        cb = np.array(tab['CD4_Count_Bins'], float)
        mort = np.array(tab['MortalityTable'], float)
        for ai, age in enumerate(['<25', '25-35', '35-45', '45+']):
            for di, dur in enumerate(durs):
                for lab, cd4 in cd4_rep.items():
                    rows.append(dict(sex=sex, arm=arm, age=age, dur=dur, cd4_0=lab, rate=np.interp(cd4, cb, mort[di, ai])))
    df = pd.DataFrame(rows)
    df.to_csv(out / 'emod_rates.csv', index=False)
    return df


# %% Tables and figures

BLUES = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#0d366b'] # Sequential, by CD4 at ART start
CATS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
HIVD = ['art_death', 'post_art_death', 'hiv_death', 'aids_death']


def cohort_rates(kind, variant='base'):
    """ HIV deaths per 100 person-years on ART (suppressed or not), by CD4 at first start and time since it """
    c = pd.read_csv(out / f'{kind}_cohort.csv')
    c = c[(c.variant == variant) & (c.status != 'off_after_art')]
    g = c.assign(hivd=c[HIVD].sum(axis=1)).groupby(['cd4_0', 'dur'])[['py', 'hivd']].sum()
    rate = (100*g.hivd/g.py).unstack()[DUR_LABELS].reindex(CD4_LABELS)
    share = (100*g.py/g.py.sum()).unstack()[DUR_LABELS].reindex(CD4_LABELS)
    return rate, share, 100*g.hivd.sum()/g.py.sum()


def malawi_summary():
    a, c = pd.read_csv(out / 'malawi_annual.csv'), pd.read_csv(out / 'malawi_cohort.csv')
    unaids = pd.read_csv(malawi_repo / 'hivsim/data/malawi_hiv_calib.csv', index_col='year')['hiv.new_deaths']
    c['hivd'] = c[HIVD].sum(axis=1)
    rows = []
    for v, g in a.groupby('variant', sort=False):
        m = g.groupby('year').mean(numeric_only=True)
        cc = c[c.variant == v]
        ns = cc.seed.nunique()
        on = cc[cc.status != 'off_after_art']
        rate = 100*on.groupby('year').hivd.sum()/on.groupby('year').py.sum()
        treated = cc.groupby('year').hivd.sum()/ns # Agents, as never_treated_hiv_deaths
        share = treated/(treated + g.groupby('year').never_treated_hiv_deaths.mean())
        row = dict(variant=v, plhiv_2021=m.loc[2021, 'plhiv'], art_coverage_2021=m.loc[2021, 'on_art']/m.loc[2021, 'plhiv'])
        for y in [2010, 2015, 2018, 2021]:
            row[f'deaths_{y}'] = m.loc[y, 'hiv_deaths']
            row[f'vs_unaids_{y}'] = m.loc[y, 'hiv_deaths']/unaids[y]
        row['vs_unaids_2015_2021'] = m.loc[2015:2021, 'hiv_deaths'].sum()/unaids.loc[2015:2021].sum()
        row['on_art_rate_2015_2021'] = 100*on[on.year.between(2015, 2021)].hivd.sum()/on[on.year.between(2015, 2021)].py.sum()
        row['share_on_or_after_art_2021'] = share[2021]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(out / 'malawi_summary.csv', index=False)
    return df


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e6e3', 'grid.linewidth': 0.6,
                         'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e', 'ytick.color': '#52514e', 'font.size': 9})

    # 1. Death rate on ART by CD4 at start and time on ART: HIVsim (plain, Malawi), EMOD's table
    emod = pd.read_csv(out / 'emod_rates.csv')
    emod = emod[(emod.sex == 'Female') & (emod.arm == 'SuppressedOnART') & (emod.age == '25-35')]
    emod_dur = {'0-6m': '0-6m', '6-12m': '6-12m', '1-2y': '1-2y', '2-5y': '2-3y', '5y+': '3y+'}
    panels = [('HIVsim, plain sim (80% coverage by 2020)', cohort_rates('hivsim')[0]),
              ('HIVsim, Malawi (exp 03 best run)', cohort_rates('malawi')[0]),
              ("EMOD's table (IeDEA-SA): suppressed,\nwomen 25-35", emod.pivot_table(index='cd4_0', columns='dur', values='rate').mul(100).rename(columns={v: k for k, v in emod_dur.items()})[DUR_LABELS].reindex(CD4_LABELS))]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True, constrained_layout=True)
    x = np.arange(len(DUR_LABELS))
    for ax, (title, tab) in zip(axes, panels):
        for i, lab in enumerate(CD4_LABELS):
            y = tab.loc[lab].values.astype(float)
            ax.plot(x, np.where(y > 0, y, np.nan), '-o', color=BLUES[i], lw=2, ms=5, label=f'CD4 {lab}')
        ax.set_yscale('log')
        ax.set_xticks(x, ['0-6 m', '6-12 m', '1-2 y', '2-5 y\n(EMOD 2-3)', '5 y+\n(EMOD 3+)'])
        ax.set_title(title, fontsize=9, loc='left')
        ax.axhspan(1, 2, color='#f0efec', zorder=0)
    axes[0].set_ylabel('HIV deaths per 100 person-years on ART')
    axes[0].text(0.05, 1.45, 'cohorts after year 1: ~1-2', fontsize=7, color='#52514e')
    axes[2].legend(title='CD4 at ART start', fontsize=7, title_fontsize=7, frameon=False, loc='upper right')
    fig.savefig(figs / 'death_rate_on_art.png', dpi=150)
    plt.close(fig)

    # 2. CD4 on ART by time since start (HIVsim, Malawi)
    cd = pd.read_csv(out / 'malawi_cd4.csv')
    cd = cd[cd.variant == 'base'].groupby(['cd4_0', 'dur'])[['cd4_sum', 'n']].sum()
    cd = (cd.cd4_sum/cd.n).unstack()[DUR_LABELS].reindex(CD4_LABELS)
    fig, ax = plt.subplots(figsize=(5, 3.4), constrained_layout=True)
    for i, lab in enumerate(CD4_LABELS):
        ax.plot(x, cd.loc[lab].values, '-o', color=BLUES[i], lw=2, ms=5, label=f'CD4 {lab}')
    for b in [50, 200, 350, 500]:
        ax.axhline(b, color='#c3c2b7', lw=0.8, ls='--', zorder=0)
    ax.set_xticks(x, ['0-6 m', '6-12 m', '1-2 y', '2-5 y', '5 y+'])
    ax.set_ylabel('Mean CD4 on ART')
    ax.set_title('HIVsim Malawi: CD4 on ART (dashed: bands of the death rate)', fontsize=9, loc='left')
    ax.legend(title='CD4 at ART start', fontsize=7, title_fontsize=7, frameon=False, loc='lower right')
    fig.savefig(figs / 'cd4_on_art.png', dpi=150)
    plt.close(fig)

    # 3. Malawi HIV deaths by variant against UNAIDS and Optima
    a = pd.read_csv(out / 'malawi_annual.csv')
    unaids = pd.read_csv(malawi_repo / 'hivsim/data/malawi_hiv_calib.csv', index_col='year')['hiv.new_deaths']
    opt = pd.read_csv(repo / 'models/optima/reference/optima_malawi_totals.csv')
    unaids_plhiv = pd.read_csv(malawi_repo / 'hivsim/data/malawi_hiv_calib.csv', index_col='year')['hiv.infected']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for j, ax in enumerate(axes):
        per = (lambda d, n: d/n*100) if j else (lambda d, n: d/1e3)
        for i, (v, g) in enumerate(a.groupby('variant', sort=False)):
            if v == 'unsupp_no_recovery':
                continue # The same as base to within 1%
            g = g.assign(y=per(g.hiv_deaths, g.plhiv))
            m = g.groupby('year').y.agg(['mean', 'min', 'max']).loc[2000:2021]
            ax.plot(m.index, m['mean'], color=CATS[i], lw=2, label=v if v != 'base' else 'base (= unsupp_no_recovery)')
            ax.fill_between(m.index, m['min'], m['max'], color=CATS[i], alpha=0.12, lw=0)
        o = opt[(opt.year >= 2000) & (opt.year <= 2021)]
        ax.plot(o.year, per(o.numdeath, o.numplhiv), color='#52514e', ls='--', lw=1.5, label='Optima (calibrated)')
        u = unaids.loc[2000:2021]
        ax.plot(u.index, per(u, unaids_plhiv.loc[2000:2021]).values, 'o', color='#0b0b0b', ms=5, label='UNAIDS')
        ax.set_ylim(0, None)
    axes[0].set_ylabel('HIV deaths per year (thousands)')
    axes[1].set_ylabel('HIV deaths per 100 PLHIV per year')
    axes[0].set_title('Malawi: HIV deaths (HIVsim: 6 seeds, range shaded)', fontsize=9, loc='left')
    axes[1].set_title('Per 100 PLHIV (HIVsim PLHIV is below UNAIDS before 2010)', fontsize=9, loc='left')
    axes[1].legend(fontsize=7, frameon=False)
    fig.savefig(figs / 'malawi_deaths_variants.png', dpi=150)
    plt.close(fig)
    return


if __name__ == '__main__':
    out.mkdir(exist_ok=True)
    figs.mkdir(exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else 'hivsim'
    if what == 'test':
        print(run_one(('hivsim', 'base', 1, dict(n_agents=3000)))[0].groupby('status')[['py'] + CAUSES].sum())
    elif what == 'hivsim':
        run('hivsim', {'base': {}}, range(1, 7))
    elif what == 'malawi':
        run('malawi', MALAWI_VARIANTS, range(19, 25))
    elif what == 'plot':
        for kind in ['hivsim', 'malawi']:
            rate, share, total = cohort_rates(kind)
            rate.round(3).to_csv(out / f'{kind}_rate_by_cd4_dur.csv')
            share.round(2).to_csv(out / f'{kind}_py_share_by_cd4_dur.csv')
            print(kind, f'HIV deaths per 100 PY on ART: {total:.3f}'); print(rate.round(2)); print(share.round(1))
        print(malawi_summary().round(3).T.to_string())
        plot()
    elif what == 'optima':
        print(optima().query('year in [2010, 2015, 2021] and group in ["untreated", "lost", "nonsuppressed", "suppressed"]').to_string())
    elif what == 'emod':
        print(emod_tables().query('age == "25-35"').pivot_table(index=['sex', 'arm', 'cd4_0'], columns='dur', values='rate').round(4).to_string())
