"""
Analysis 03: EMOD′ HIV (models/emod/v4/hiv.py) against EMOD with the campaign reduced to what EMOD′ has. See README.md.

    python run.py          # run EMOD′ (3 sets × 5 seeds), reduce it and EMOD's runs (raw_results/emod/nh3) to outputs/
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

SETS = [0, 100, 200]
SEEDS = [1, 2, 3, 4, 5]
EMOD = repo / 'raw_results' / 'emod' / 'nh3'
START = 1960.5
AGES = np.arange(0, 70, 5)
STAGES = {'ACUTE': 0, 'LATENT': 1, 'AIDS': 2} # EMOD's 1-3
YEARBINS = np.arange(0, 41) # Years from infection to death


def agebin(age):
    return AGES[np.clip(np.searchsorted(AGES, age, side='right') - 1, 0, len(AGES) - 1)]


def emod_tables(pid):
    """ EMOD's reports as common tables: yearly snapshots, infections (sexual), deaths from HIV """
    out = EMOD / f'set{pid:03d}' / 'output'
    d = pd.read_csv(out / 'ReportHIVByAgeAndGender.csv')
    d.columns = d.columns.str.strip()
    d = d[np.isclose(d.Year % 1, 0.92)]
    st = d.HIV_Stage.map(STAGES)
    snap = pd.DataFrame(dict(year=d.Year, node=d.NodeId, female=d.Gender, age=d.Age, pop=d.Population, infected=d.Infected,
                             acute=d.Infected*(st == 0), latent=d.Infected*(st == 1), aids=d.Infected*(st == 2), circumcised=d.Population*d.IsCircumcised*(d.Gender == 0)))
    snap = snap.groupby(['year', 'node', 'female', 'age']).sum().reset_index()
    cols = {'YEAR': 'year', 'NODE_ID': 'node', 'REL_TYPE (0 = TRANSITORY; 1 = INFORMAL; 2 = MARITAL; 3 = COMMERCIAL)': 'type', 'SRC_AGE': 'src_age',
            'SRC_INF_AGE': 'src_inf_age', 'SRC_STAGE': 'src_stage', 'DEST_GENDER': 'female', 'DEST_AGE': 'age'}
    t = pd.read_csv(out / 'TransmissionReport.csv', usecols=list(cols)).rename(columns=cols)
    inf = pd.DataFrame(dict(year=t.year, type=t.type, node=t.node, female=t.female, age=t.age/365, src_stage=t.src_stage - 1, src_days=t.src_age - t.src_inf_age))
    m = pd.read_csv(out / 'HIVMortality.csv', usecols=['Node_ID', 'Death_time', 'Death_was_HIV_cause', 'Gender', 'Age', 'Years_since_infection'])
    m = m[m.Death_was_HIV_cause == 1]
    deaths = pd.DataFrame(dict(year=START + m.Death_time/365, age=m.Age, female=m.Gender, node=m.Node_ID, years_infected=m.Years_since_infection))
    return snap, inf, deaths


def v4_tables(job):
    import harness
    pid, seed = job
    sim = harness.make_sim(pid, seed=seed, hiv=True)
    sim.run()
    h = sim.diseases.hiv
    inf, deaths = h.infections, h.deaths
    inf = inf.assign(node=inf.node + 1, female=inf.female.astype(int))
    deaths = deaths.assign(node=deaths.node + 1, female=deaths.female.astype(int))
    snap = sim.analyzers.snapshots.df.groupby(['year', 'node', 'female', 'age']).sum(numeric_only=True).reset_index()
    return snap, inf, deaths


def summarize(snap, inf, deaths):
    """ Small tables from the common ones """
    out = dict(snap=snap[['year', 'node', 'female', 'age', 'pop', 'infected', 'acute', 'latent', 'aids', 'circumcised']])
    yr = np.floor(inf.year).astype(int).rename('year')
    sex = inf[inf.type >= 0] # Sexual transmissions (EMOD's report has only these)
    syr = np.floor(sex.year).astype(int).rename('year')
    out['incidence'] = inf.groupby([yr, 'type', 'node', 'female', agebin(inf.age)]).size().rename('n').reset_index().rename(columns={'age': 'age'})
    out['by_source'] = sex.groupby([syr, 'type', 'src_stage', (sex.src_days//365).clip(upper=30).rename('src_years')]).size().rename('n').reset_index()
    dyr = np.floor(deaths.year).astype(int).rename('year')
    out['deaths'] = deaths.groupby([dyr, 'node', 'female', agebin(deaths.age), np.floor(deaths.years_infected).clip(upper=40).rename('years_infected')]).size().rename('n').reset_index()
    return out


def run_one(job):
    model, pid, seed = job
    tabs = emod_tables(pid) if model == 'emod' else v4_tables((pid, seed))
    out = summarize(*tabs)
    for df in out.values():
        df.insert(0, 'seed', seed)
        df.insert(0, 'pid', pid)
        df.insert(0, 'model', model)
    print(model, pid, seed, 'done', flush=True)
    return out


def reduce(cores):
    jobs = [('emod', p, 0) for p in SETS] + [('v4', p, s) for p in SETS for s in SEEDS]
    with Pool(cores) as pool:
        res = pool.map(run_one, jobs, chunksize=1)
    (here / 'outputs').mkdir(exist_ok=True)
    for name in res[0]:
        pd.concat([r[name] for r in res], ignore_index=True).to_csv(here / 'outputs' / f'{name}.csv.gz', index=False)
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_muted
    o = {n: pd.read_csv(here / 'outputs' / f'{n}.csv.gz') for n in ['snap', 'incidence', 'by_source', 'deaths']}
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)
    shades = {0: '-', 100: '--', 200: ':'}
    note = 'Each parameter set: EMOD (orange, one run) and EMOD′ (blue, median of 5 seeds, band their range); line style by set.'
    tnames = ['Transitory', 'Informal', 'Marital', 'Commercial']

    def save(fig, name, text=note):
        fig.text(0.01, -0.04, text, fontsize=7, color=c_muted, va='top')
        fig.savefig(figdir / name, dpi=150, bbox_inches='tight')
        plt.close(fig)

    def lines(ax, df, x, y):
        """ EMOD's run and EMOD′'s seeds for each set """
        for pid, ls in shades.items():
            e = df[(df.model == 'emod') & (df.pid == pid)].groupby(x)[y].sum()
            ax.plot(e.index, e, color=c_best, ls=ls, lw=1.3, label=f'EMOD, set {pid}')
            v = df[(df.model == 'v4') & (df.pid == pid)].groupby([x, 'seed'])[y].sum().unstack()
            ax.fill_between(v.index, v.min(axis=1), v.max(axis=1), color=c_run, alpha=0.12, lw=0)
            ax.plot(v.index, v.median(axis=1), color=c_run, ls=ls, lw=1.3, label=f'EMOD′, set {pid}')

    s = o['snap']
    adult = s[(s.age >= 15) & (s.age < 50)]
    # Totals over time
    fig, axs = plt.subplots(2, 3, figsize=(13, 6))
    for f, ax in zip((0, 1), axs[0, :2]):
        a = adult[adult.female == f].groupby(['model', 'pid', 'seed', 'year'])[['infected', 'pop']].sum().reset_index()
        lines(ax, a.assign(prev=a.infected/a['pop']), 'year', 'prev')
        style(ax, f'Prevalence 15-49, {"women" if f else "men"}', 'Share')
    plh = s.groupby(['model', 'pid', 'seed', 'year'])['infected'].sum().reset_index()
    lines(axs[0, 2], plh, 'year', 'infected')
    style(axs[0, 2], 'People with HIV (agents)')
    inc = o['incidence']
    sexual = inc[inc.type >= 0].groupby(['model', 'pid', 'seed', 'year']).n.sum().reset_index()
    lines(axs[1, 0], sexual, 'year', 'n')
    style(axs[1, 0], 'New sexual infections a year (agents)')
    vert = inc[inc.type == -2].groupby(['model', 'pid', 'seed', 'year']).n.sum().reset_index()
    lines(axs[1, 1], vert, 'year', 'n')
    style(axs[1, 1], 'Infections at birth a year (EMOD′ only)')
    de = o['deaths'].groupby(['model', 'pid', 'seed', 'year']).n.sum().reset_index()
    lines(axs[1, 2], de, 'year', 'n')
    style(axs[1, 2], 'Deaths from HIV a year (agents)')
    for ax in axs.flat:
        ax.set_xlim(1978, 2025)
    axs[0, 0].legend(fontsize=6, frameon=False, ncol=2)
    save(fig, 'totals.png')

    # Prevalence by age and sex
    fig, axs = plt.subplots(2, 3, figsize=(13, 6), sharey='row')
    for j, yr in enumerate([1995.92, 2005.92, 2020.92]):
        for f in (0, 1):
            d = s[np.isclose(s.year, yr) & (s.female == f)].groupby(['model', 'pid', 'seed', 'age'])[['infected', 'pop']].sum().reset_index()
            d = d.assign(prev=d.infected/d['pop'].clip(lower=1), age=d.age + 2.5)
            lines(axs[f, j], d, 'age', 'prev')
            style(axs[f, j], f'{int(yr)}, {"women" if f else "men"}', 'Prevalence' if j == 0 else None)
    for ax in axs[1]:
        ax.set_xlabel('Age', fontsize=8, color=c_muted)
    axs[0, 0].legend(fontsize=6, frameon=False, ncol=2)
    save(fig, 'prevalence_by_age.png')

    # Prevalence by node
    fig, axs = plt.subplots(1, 3, figsize=(13, 3))
    for n, ax in zip((1, 2, 3), axs):
        a = adult[adult.node == n].groupby(['model', 'pid', 'seed', 'year'])[['infected', 'pop']].sum().reset_index()
        lines(ax, a.assign(prev=a.infected/a['pop']), 'year', 'prev')
        style(ax, f'Prevalence 15-49, node {n} ({["Northern", "Central", "Southern"][n-1]})', 'Share' if n == 1 else None)
        ax.set_xlim(1978, 2025)
    axs[0].legend(fontsize=6, frameon=False, ncol=2)
    save(fig, 'prevalence_by_node.png')

    # Stage of people with HIV; years from infection to death
    fig, axs = plt.subplots(1, 3, figsize=(13, 3))
    st = s.groupby(['model', 'pid', 'seed', 'year'])[['infected', 'acute', 'aids']].sum().reset_index()
    lines(axs[0], st.assign(r=st.acute/st.infected.clip(lower=1)), 'year', 'r')
    style(axs[0], 'Acute, share of people with HIV')
    lines(axs[1], st.assign(r=st.aids/st.infected.clip(lower=1)), 'year', 'r')
    style(axs[1], 'AIDS, share of people with HIV')
    for ax in axs[:2]:
        ax.set_xlim(1985, 2025)
    d = o['deaths'][o['deaths'].year.between(2000, 2024)].groupby(['model', 'pid', 'seed', 'years_infected']).n.sum().reset_index()
    d['share'] = d.n/d.groupby(['model', 'pid', 'seed']).n.transform('sum')
    lines(axs[2], d, 'years_infected', 'share')
    style(axs[2], 'Years from infection to HIV death, deaths 2000-24', 'Share')
    axs[2].set_xlabel('Whole years since infection', fontsize=8, color=c_muted)
    axs[0].legend(fontsize=6, frameon=False, ncol=2)
    save(fig, 'natural_history.png')

    # Transmissions by relationship type and stage of the source
    b = o['by_source']
    fig, axs = plt.subplots(2, 4, figsize=(13, 5.5))
    tot = b.groupby(['model', 'pid', 'seed', 'year']).n.sum().rename('tot')
    for i in range(4):
        d = b[b.type == i].groupby(['model', 'pid', 'seed', 'year']).n.sum().reset_index().merge(tot.reset_index(), on=['model', 'pid', 'seed', 'year'])
        lines(axs[0, i], d.assign(r=d.n/d.tot), 'year', 'r')
        style(axs[0, i], f'{tnames[i]}: share of sexual infections', 'Share' if i == 0 else None)
        axs[0, i].set_xlim(1982, 2025)
    for i, (k, name) in enumerate([(0, 'Acute'), (1, 'Latent'), (2, 'AIDS')]):
        d = b[b.src_stage == k].groupby(['model', 'pid', 'seed', 'year']).n.sum().reset_index().merge(tot.reset_index(), on=['model', 'pid', 'seed', 'year'])
        lines(axs[1, i], d.assign(r=d.n/d.tot), 'year', 'r')
        style(axs[1, i], f'Partner who transmitted: {name}', 'Share' if i == 0 else None)
        axs[1, i].set_xlim(1982, 2025)
    d = b[b.year.between(2000, 2024) & (b.src_stage >= 0)].groupby(['model', 'pid', 'seed', 'src_years']).n.sum().reset_index()
    d['share'] = d.n/d.groupby(['model', 'pid', 'seed']).n.transform('sum')
    lines(axs[1, 3], d, 'src_years', 'share')
    style(axs[1, 3], 'Partner\'s years since infection, 2000-24', 'Share')
    axs[0, 0].legend(fontsize=6, frameon=False, ncol=2)
    fig.subplots_adjust(hspace=0.45)
    save(fig, 'transmissions.png', note + ' Infections in a relationship only (EMOD\'s TransmissionReport); seeded and at-birth infections excluded.')
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=18)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
