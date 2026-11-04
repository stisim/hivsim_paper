"""
Analysis 02: the EMOD′ network (models/emod/v4) against EMOD's, for 10 of Daniel's parameter sets.

Each run (EMOD with HIV, EMOD without HIV, EMOD′) is reduced to the same summary tables in outputs/, one row per
run and stratum, with columns `model` and `pid`. See README.md.

    python run.py          # run EMOD′, reduce its runs and EMOD's (raw_results/emod/net10, net3_nohiv) to outputs/
    python run.py --demog  # EMOD′'s demography alone, 20 seeds, to outputs/population.csv.gz
    python run.py --plot   # figures/ from outputs/
"""
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo / 'models' / 'emod' / 'v4'))
sys.path.insert(0, str(repo / 'common'))
import inputs # noqa: E402

SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
NOHIV = [0, 100, 200]
EMOD = repo / 'raw_results' / 'emod'
TYPES = inputs.TYPES
KEYS = ['T', 'I', 'M', 'C']
RISKS = {'LOW': 0, 'MEDIUM': 1, 'HIGH': 2}
START = 1960.5
AGEDIFF = np.arange(-20, 41) # Bin edges, years
DURBINS = np.unique(np.round(np.geomspace(1, 1e5, 61))) # Days
ACTBINS = np.array([0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e5])
PERIODS = [1960, 1975, 1990, 2005, 2025]


# %% Reading runs into common tables

def emod_tables(folder, pid):
    """ EMOD's reports as the common tables: starts, ends and yearly snapshots """
    out = EMOD / folder / f'set{pid:03d}' / 'output'
    cols = {'Rel_start_time': 't', 'Rel_scheduled_end_time': 'end', 'Rel_type (0 = TRANSITORY; 1 = INFORMAL; 2 = MARITAL; 3 = COMMERCIAL)': 'type',
            'Original_node_ID': 'node', 'A_ID': 'm', 'A_age': 'age_m', "A_IP='Risk'": 'risk_m', 'B_ID': 'f', 'B_age': 'age_f', "B_IP='Risk'": 'risk_f', 'A_gender': 'a_gender'}
    s = pd.read_csv(out / 'RelationshipStart.csv', usecols=list(cols)).rename(columns=cols)
    assert (s.a_gender == 0).all()
    starts = pd.DataFrame(dict(year=START + s.t/365, type=s.type, node=s.node, age_m=s.age_m, age_f=s.age_f,
                               risk_m=s.risk_m.map(RISKS), risk_f=s.risk_f.map(RISKS), dur_days=s.end - s.t, m=s.m, f=s.f))
    cols = {'Rel_start_time': 't0', 'Rel_actual_end_time': 't1', 'Rel_type (0 = TRANSITORY; 1 = INFORMAL; 2 = MARITAL; 3 = COMMERCIAL)': 'type',
            'num_total_coital_acts': 'acts', 'Termination_Reason': 'reason'}
    e = pd.read_csv(out / 'RelationshipEnd.csv', usecols=list(cols)).rename(columns=cols)
    ends = pd.DataFrame(dict(year=START + e.t0/365, type=e.type, days=e.t1 - e.t0, acts=e.acts, died=e.reason != 'BROKEUP'))
    d = pd.read_csv(out / 'ReportHIVByAgeAndGender.csv')
    d.columns = d.columns.str.strip()
    d = d[np.isclose(d.Year % 1, 0.92)]
    snap = pd.DataFrame(dict(year=d.Year, node=d.NodeId, female=d.Gender, risk=d['IP_Key:Risk'].map(RISKS), age=d.Age, pop=d.Population,
                             debuted=d['Has Debuted'], concurrent=d['Has Concurrent Partners'], partners=d['Current Partners']))
    for k, t in zip(KEYS, TYPES):
        snap[f'cur_{k}'] = d[f'Currently ({t})']
        snap[f'num_{k}'] = d[f'Num_{t}']
        snap[f'que_{k}'] = d[f'Considering ({t})']
        snap[f'avail_{k}'] = d[f'Available for Relationship ({t})']
    snap = snap.groupby(['year', 'node', 'female', 'risk', 'age']).sum().reset_index() # Over HIV and ART status
    return starts, ends, snap


def v4_tables(pid):
    """ Run EMOD′'s network for one parameter set; its records as the common tables """
    import harness
    inp = inputs.load(pid)
    sim = harness.make_sim(pid, seed=inp.run_number % 2**31)
    sim.run()
    net, dt = sim.networks.emodnet, inp.dt
    s, e = net.starts, net.ends
    starts = pd.DataFrame(dict(year=START + (s.ti - 1)/12, type=s.type, node=s.node + 1, age_m=s.age_m, age_f=s.age_f,
                               risk_m=s.risk_m, risk_f=s.risk_f, dur_days=s.dur_days, m=s.p1, f=s.p2))
    ends = pd.DataFrame(dict(year=START + (e.start - 1)/12, type=e.type, days=(e.end - e.start)*dt, acts=e.acts, died=e.died))
    return starts, ends, sim.analyzers.snapshots.df


# %% Summaries, the same for both models

def hist(x, bins):
    return np.histogram(np.clip(x, bins[0], bins[-1] - 1e-9), bins=bins)[0]


def summarize(starts, ends, snap):
    """ Small tables from the common tables of one run """
    out = {}
    out['snap'] = snap
    yr = np.floor(starts.year).astype(int)
    out['formed'] = starts.groupby([yr.rename('year'), 'type']).size().rename('n').reset_index()
    period = pd.cut(starts.year, PERIODS, right=False, labels=PERIODS[:-1]).astype(int).rename('period')
    rows = []
    for (p, t), g in starts.groupby([period, 'type']):
        rows += [dict(period=p, type=t, what='agediff', bin=b, n=n) for b, n in zip(AGEDIFF[:-1], hist(g.age_m - g.age_f, AGEDIFF))]
        rows += [dict(period=p, type=t, what='dur_sched', bin=b, n=n) for b, n in zip(DURBINS[:-1], hist(g.dur_days, DURBINS))]
    eperiod = pd.cut(ends.year, PERIODS, right=False, labels=PERIODS[:-1]).astype(int).rename('period')
    for (p, t), g in ends.groupby([eperiod, 'type']):
        rows += [dict(period=p, type=t, what='dur_actual', bin=b, n=n) for b, n in zip(DURBINS[:-1], hist(g.days, DURBINS))]
        rows += [dict(period=p, type=t, what='acts', bin=b, n=n) for b, n in zip(ACTBINS[:-1], hist(g.acts, ACTBINS))]
        rows += [dict(period=p, type=t, what='died', bin=0, n=int(g.died.sum()))]
    out['dists'] = pd.DataFrame(rows)
    mean = lambda df, col: df.groupby([np.floor(df.year).astype(int).rename('year'), 'type'])[col].mean().rename(col).reset_index()
    out['means'] = mean(starts.assign(agediff=starts.age_m - starts.age_f), 'agediff').merge(mean(starts, 'dur_days'), on=['year', 'type']).merge(
        mean(starts, 'age_f'), on=['year', 'type']).merge(mean(starts, 'age_m'), on=['year', 'type'])
    c = starts[starts.type == 3]
    cy = np.floor(c.year).astype(int).rename('year')
    out['commercial'] = pd.DataFrame(dict(n=c.groupby(cy).size(), fsw_any=c.groupby(cy).f.nunique(), client_any=c.groupby(cy).m.nunique())).reset_index()
    out['mixing'] = starts[starts.year >= 1975].groupby(['type', 'risk_m', 'risk_f']).size().rename('n').reset_index()
    return out


def run_one(job):
    model, pid = job
    tabs = v4_tables(pid) if model == 'v4' else emod_tables(dict(emod='net10', emod_nohiv='net3_nohiv')[model], pid)
    out = summarize(*tabs)
    for df in out.values():
        df.insert(0, 'pid', pid)
        df.insert(0, 'model', model)
    print(f'{model} {pid} done', flush=True)
    return out


def reduce(cores=8):
    jobs = [('v4', p) for p in SETS] + [('emod', p) for p in SETS] + [('emod_nohiv', p) for p in NOHIV]
    with Pool(cores) as pool:
        res = pool.map(run_one, jobs, chunksize=1)
    (here / 'outputs').mkdir(exist_ok=True)
    for name in res[0]:
        pd.concat([r[name] for r in res], ignore_index=True).to_csv(here / 'outputs' / f'{name}.csv.gz', index=False)
    return


def demog_one(seed):
    import harness
    sim = harness.make_sim(0, seed=seed, network=False)
    sim.run()
    return sim.analyzers.snapshots.df.assign(model='v4_demog', pid=seed)


def reduce_demog(cores=8, seeds=20):
    with Pool(cores) as pool:
        res = pool.map(demog_one, range(1, seeds + 1))
    pd.concat(res, ignore_index=True).to_csv(here / 'outputs' / 'population.csv.gz', index=False)
    return


# %% Figures

def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_data, c_muted, c_band
    o = {n: pd.read_csv(here / 'outputs' / f'{n}.csv.gz') for n in ['snap', 'formed', 'dists', 'means', 'commercial', 'mixing']}
    for n, df in o.items(): # EMOD′ on the 3 sets run without HIV in EMOD, for like-for-like lines
        o[n] = pd.concat([df, df[(df.model == 'v4') & df.pid.isin(NOHIV)].assign(model='v4_3')], ignore_index=True)
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)
    looks = dict(emod=dict(color=c_best, label='EMOD with HIV, 10 sets'), v4=dict(color=c_run, label='EMOD′ (no HIV), 10 sets'),
                 emod_nohiv=dict(color=c_best, ls='--', label='EMOD without HIV, 3 sets'), v4_3=dict(color=c_run, ls='--', label='EMOD′, the same 3 sets'))
    tnames = [t.title() for t in TYPES]

    def band(ax, df, x, y, model):
        """ Median over parameter sets, and their range """
        g = df.groupby(x)[y]
        lo, mid, hi = g.min(), g.median(), g.max()
        lk = looks[model]
        if model in ('emod', 'v4'):
            ax.fill_between(mid.index, lo, hi, color=lk['color'], alpha=0.15, lw=0)
        ax.plot(mid.index, mid, lw=1.5, **lk)

    def save(fig, name, note=None):
        if note:
            fig.text(0.01, -0.06, note, fontsize=7, color=c_muted, va='bottom')
        fig.savefig(figdir / name, dpi=150, bbox_inches='tight')
        plt.close(fig)

    snap = o['snap']
    tot = lambda by: snap.groupby(['model', 'pid', 'year'] + by).sum(numeric_only=True).reset_index()
    note = 'Solid: median over the 10 parameter sets, bands their range; dashed: median over sets 0, 100 and 200, the ones EMOD was also run for without HIV.'

    # Population, FSW and clients
    fig, axs = plt.subplots(1, 3, figsize=(11, 3))
    t = tot([])
    hi = tot(['risk', 'female'])
    hi = hi[hi.risk == 2]
    for m in looks:
        band(axs[0], t[t.model == m], 'year', 'pop', m)
        band(axs[1], hi[(hi.model == m) & (hi.female == 1)], 'year', 'pop', m)
        band(axs[2], hi[(hi.model == m) & (hi.female == 0)], 'year', 'pop', m)
    for ax, ti in zip(axs, ['Population (agents)', 'Female sex workers (Risk HIGH women)', 'Clients (Risk HIGH men)']):
        style(ax, ti)
    axs[0].legend(fontsize=7, frameon=False)
    save(fig, 'population.png', note)

    # The demography alone (20 seeds) against EMOD: total, and by age and sex
    pop = pd.concat([pd.read_csv(here / 'outputs' / 'population.csv.gz'), snap[snap.model.isin(['emod', 'emod_nohiv'])]])
    pa = pop.groupby(['model', 'pid', 'year', 'female', 'age'])['pop'].sum().reset_index()
    fig, axs = plt.subplots(1, 3, figsize=(12, 3))
    t = pa.groupby(['model', 'pid', 'year'])['pop'].sum().reset_index()
    ref = t[t.model == 'emod_nohiv'].groupby('year')['pop'].mean()
    for m, lk in [('v4_demog', dict(color=c_run, label='EMOD′ demography (20 seeds)')), ('emod', looks['emod']), ('emod_nohiv', dict(color=c_best, ls='--', marker='o', ms=3, lw=0, label='EMOD, no HIV (3 runs)'))]:
        d = t[t.model == m].pivot(index='year', columns='pid', values='pop').div(ref, axis=0)
        if m == 'emod_nohiv':
            for c in d:
                axs[0].plot(d.index, d[c], color=c_best, marker='o', ms=1.5, lw=0)
            axs[0].plot([], [], **lk)
        else:
            axs[0].fill_between(d.index, d.min(axis=1), d.max(axis=1), color=lk['color'], alpha=0.15, lw=0)
            axs[0].plot(d.index, d.median(axis=1), **{k: v for k, v in lk.items() if k != 'ls'}, ls=lk.get('ls', '-'))
    style(axs[0], 'Population, relative to EMOD without HIV', 'Ratio to the mean of the 3 runs')
    axs[0].legend(fontsize=7, frameon=False)
    for ax, yr, refm in [(axs[1], 1984.92, 'emod'), (axs[2], 2023.92, 'emod_nohiv')]:
        d = pa[np.isclose(pa.year, yr)]
        r = d[d.model == refm].groupby(['female', 'age'])['pop'].mean()
        v = d[d.model == 'v4_demog'].groupby(['pid', 'female', 'age'])['pop'].sum().unstack(['female', 'age']).div(r)
        for f, col in ((0, c_run), (1, c_data)):
            q = v.xs(f, axis=1, level='female')
            ax.fill_between(q.columns + 2.5, q.quantile(0.05), q.quantile(0.95), color=col, alpha=0.12, lw=0)
            ax.plot(q.columns + 2.5, q.median(), color=col, marker='o', ms=3, lw=1.2, label='Women' if f else 'Men')
        ax.axhline(1, color=c_muted, lw=0.6)
        style(ax, f'By age, {int(yr)}: EMOD′ / EMOD{"" if refm == "emod" else " without HIV"}', 'Ratio')
        ax.set_xlabel('Age', fontsize=8, color=c_muted)
    axs[1].legend(fontsize=7, frameon=False)
    save(fig, 'demography.png', 'EMOD′ demography alone, 20 seeds (median, 5-95%), against the mean of EMOD\'s runs: in 1984, all 13 (HIV has hardly acted yet); in 2023, the 3 without HIV.')

    # In a relationship of each type, by sex, among those past debut
    fig, axs = plt.subplots(2, 4, figsize=(13, 5.5), sharex=True)
    t = tot(['female'])
    for f in (0, 1):
        for i, k in enumerate(KEYS):
            ax = axs[f, i]
            for m in looks:
                d = t[(t.model == m) & (t.female == f)].assign(share=lambda d: d[f'cur_{k}']/d.debuted)
                band(ax, d, 'year', 'share', m)
            style(ax, f'{tnames[i]}, {"women" if f else "men"}', 'Share of those past debut' if i == 0 else None)
    axs[0, 0].legend(fontsize=7, frameon=False)
    save(fig, 'in_relationship.png', note)

    # By age, in 1985 (10 sets, before HIV matters) and 2020 (the 3 sets without HIV)
    fig, axs = plt.subplots(2, 4, figsize=(13, 5.5), sharex=True)
    t = tot(['female', 'age'])
    for f in (0, 1):
        for i, k in enumerate(KEYS):
            ax = axs[f, i]
            for yr, alpha, models in ((1985.92, 0.5, ['emod', 'v4']), (2020.92, 1, ['emod_nohiv', 'v4_3'])):
                for m in models:
                    d = t[(t.model == m) & (t.female == f) & np.isclose(t.year, yr)].assign(share=lambda d: d[f'cur_{k}']/d['pop'].clip(lower=1))
                    g = d.groupby('age').share.median()
                    ax.plot(g.index + 2.5, g, color=looks[m]['color'], ls=looks[m].get('ls', '-'), alpha=alpha, lw=1.5, marker='o', ms=2, label=f'{looks[m]["label"]}, {int(yr)}')
            style(ax, f'{tnames[i]}, {"women" if f else "men"}', 'Share of the age group' if i == 0 else None)
    axs[0, 0].legend(fontsize=6, frameon=False)
    for ax in axs[1]:
        ax.set_xlabel('Age', fontsize=8, color=c_muted)
    save(fig, 'by_age.png', 'Median over parameter sets: in 1985, the 10 (EMOD with HIV, which has hardly acted yet); in 2020, the 3 EMOD was run for without HIV.')

    # Concurrency and number of partners
    fig, axs = plt.subplots(2, 2, figsize=(9, 5.5), sharex=True)
    t = tot(['female'])
    for f in (0, 1):
        for m in looks:
            d = t[(t.model == m) & (t.female == f)]
            band(axs[f, 0], d.assign(share=d.concurrent/d.debuted), 'year', 'share', m)
            band(axs[f, 1], d.assign(mean=d.partners/d.debuted), 'year', 'mean', m)
        sex = 'women' if f else 'men'
        style(axs[f, 0], f'Two or more partners at once, {sex}', 'Share of those past debut')
        style(axs[f, 1], f'Current partners (all types), {sex}', 'Mean, those past debut')
    axs[0, 0].legend(fontsize=7, frameon=False)
    save(fig, 'concurrency.png', note)

    # Formation rates
    fig, axs = plt.subplots(1, 4, figsize=(13, 3))
    deb = tot([]).assign(year=lambda d: np.floor(d.year).astype(int))[['model', 'pid', 'year', 'debuted']]
    fr = o['formed'].merge(deb, on=['model', 'pid', 'year'])
    for i in range(4):
        for m in looks:
            d = fr[(fr.model == m) & (fr.type == i)]
            band(axs[i], d.assign(rate=1000*d.n/d.debuted), 'year', 'rate', m)
        style(axs[i], tnames[i], 'Formed per year per 1000 past debut' if i == 0 else None)
    axs[0].legend(fontsize=7, frameon=False)
    save(fig, 'formation.png', note + ' Per 1000 people past debut at the end of the year.')

    # Queues
    fig, axs = plt.subplots(1, 4, figsize=(13, 3))
    t = tot([])
    for i, k in enumerate(KEYS):
        for m in looks:
            d = t[t.model == m]
            band(axs[i], d.assign(q=1000*d[f'que_{k}']/d.debuted), 'year', 'q', m)
        style(axs[i], tnames[i], 'Queue entries per 1000 past debut' if i == 0 else None)
    axs[0].set_yscale('symlog', linthresh=1)
    axs[3].set_yscale('log')
    axs[0].legend(fontsize=7, frameon=False)
    save(fig, 'queues.png', note)

    # Distributions: durations (scheduled and actual), age differences, acts
    di = o['dists']
    def dist_panel(ax, what, i, period=None, cdf=True, logx=True):
        for m in looks:
            d = di[(di.model == m) & (di.what == what) & (di.type == i)]
            if period is not None:
                d = d[d.period == period]
            g = d.groupby(['pid', 'bin']).n.sum().unstack(fill_value=0)
            if g.empty:
                continue
            p = g.div(g.sum(axis=1), axis=0)
            y = p.cumsum(axis=1) if cdf else p
            ax.plot(y.columns, y.median(axis=0), lw=1.5, drawstyle='steps-post' if cdf else 'default', **looks[m])
        if logx:
            ax.set_xscale('log')
    fig, axs = plt.subplots(3, 4, figsize=(13, 8))
    for i in range(4):
        dist_panel(axs[0, i], 'dur_sched', i)
        style(axs[0, i], f'{tnames[i]}: drawn duration', 'Cumulative share' if i == 0 else None)
        dist_panel(axs[1, i], 'dur_actual', i)
        style(axs[1, i], f'{tnames[i]}: actual, ended by 2025', 'Cumulative share' if i == 0 else None)
        dist_panel(axs[2, i], 'acts', i)
        style(axs[2, i], f'{tnames[i]}: coital acts, ended', 'Cumulative share' if i == 0 else None)
        axs[1, i].set_xlabel('Days', fontsize=8, color=c_muted)
        axs[2, i].set_xlabel('Acts in the relationship', fontsize=8, color=c_muted)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.subplots_adjust(hspace=0.55)
    save(fig, 'durations.png', 'Median over parameter sets of the cumulative distribution, all relationships 1960-2025.')

    fig, axs = plt.subplots(2, 4, figsize=(13, 5.5))
    for i in range(4):
        dist_panel(axs[0, i], 'agediff', i, cdf=False, logx=False)
        style(axs[0, i], f'{tnames[i]}', 'Share of relationships' if i == 0 else None)
        axs[0, i].set_xlabel('Man\'s age minus woman\'s (years)', fontsize=8, color=c_muted)
        for m in looks:
            d = o['means'][(o['means'].model == m) & (o['means'].type == i)]
            band(axs[1, i], d, 'year', 'agediff', m)
        style(axs[1, i], f'{tnames[i]}: mean age difference', 'Years' if i == 0 else None)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.subplots_adjust(hspace=0.5)
    save(fig, 'age_difference.png', 'Top: all relationships 1960-2025, median over parameter sets. Bottom: by year formed; ' + note.lower())

    # Commercial contacts per FSW and per client
    hi = tot(['risk', 'female'])
    hi = hi[hi.risk == 2].assign(year=lambda d: np.floor(d.year).astype(int)).pivot_table(index=['model', 'pid', 'year'], columns='female', values='pop').reset_index()
    hi.columns = ['model', 'pid', 'year', 'clients', 'fsw']
    cm = o['commercial'].merge(hi, on=['model', 'pid', 'year'])
    fig, axs = plt.subplots(1, 3, figsize=(11, 3))
    for m in looks:
        d = cm[(cm.model == m) & (cm.year >= 1965)]
        band(axs[0], d.assign(r=d.n/d.fsw), 'year', 'r', m)
        band(axs[1], d.assign(r=d.n/d.clients), 'year', 'r', m)
        band(axs[2], d.assign(r=d.fsw_any/d.fsw), 'year', 'r', m)
    style(axs[0], 'Commercial contacts per FSW', 'Per year')
    style(axs[1], 'Commercial contacts per client', 'Per year')
    style(axs[2], 'FSW with any contact in the year', 'Share of HIGH women at year end')
    axs[0].legend(fontsize=7, frameon=False)
    save(fig, 'commercial.png', note + ' FSW and clients: Risk HIGH at the end of the year.')
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--demog', action='store_true')
    ap.add_argument('--cores', type=int, default=8)
    a = ap.parse_args()
    plot() if a.plot else reduce_demog(a.cores) if a.demog else reduce(a.cores)
