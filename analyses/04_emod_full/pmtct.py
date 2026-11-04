"""
PMTCT and InterventionStatus: EMOD (sets 0, 100, 200 with InterventionStatus in the report, raw_results/emod/ip3)
against EMOD′ as it is and in two variants, to test two behaviours EMOD′ gets from the campaign:
- `unblocked`: the combination-PMTCT decision (Needs_Combination_PMTCT) not disqualified by ART staging, so women
  found at ANC after 2008 get Option B;
- `persist`: PMTCT efficacy kept after its 280 days, instead of reset (PMTCT::Update).
Outputs: outputs/pmtct_status.csv.gz (InterventionStatus by year, sex and HIV), outputs/pmtct_infants.csv.gz
(infections at ages 0-4 by year). The variants are made here, for this test only.

    python pmtct.py          # runs EMOD′ (3 sets × 3 seeds × 3 variants) and reduces EMOD
    python pmtct.py --plot
"""
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
for p in ('models/emod/v4', 'models/emod/reference', 'common'):
    sys.path.insert(0, str(repo / p))
SETS, SEEDS, VARIANTS = [0, 100, 200], [1, 2, 3], ['as_is', 'unblocked', 'persist']
EMOD = repo / 'raw_results' / 'emod' / 'ip3'


def tables(d):
    d = d[d.Year >= 1995]
    s = d[np.isclose(d.Year % 1, 0.92)].assign(year=lambda x: np.floor(x.Year).astype(int))
    status = s.groupby(['year', 'Gender', 'HasHIV', 'InterventionStatus']).Population.sum().rename('n').reset_index()
    kids = d[d.Age == 0]
    infants = kids.groupby(np.floor(kids.Year).astype(int).rename('year'))['Newly Infected'].sum().rename('new_0_4').to_frame()
    infants['infected_0_4'] = kids[np.isclose(kids.Year % 1, 0.92)].groupby(np.floor(kids.Year).astype(int)).Infected.sum()
    infants['pop_0_4'] = kids[np.isclose(kids.Year % 1, 0.92)].groupby(np.floor(kids.Year).astype(int)).Population.sum()
    return status, infants.reset_index()


def run_v4(job):
    import harness, inputs, cascade
    pid, seed, variant = job
    if variant == 'persist':
        fire = cascade.EmodCascade.fire
        cascade.EmodCascade.fire = lambda self, tid, uids, year: [] if self.templates[tid]['class'] == 'PMTCT' else fire(self, tid, uids, year)
    sim = harness.make_sim(pid, seed=inputs.load(pid).run_number*10 + seed, cascade=True, record=False)
    if variant == 'unblocked':
        sim.init()
        cas = sim.interventions.emodcascade
        for li in cas.listeners:
            if 'Needs_Combination_PMTCT' in li['triggers']:
                cas.templates[li['tid']]['disq'] = []
    sim.run()
    status, infants = tables(sim.analyzers.emodreport.df)
    return [t.assign(model=f'v4_{variant}', pid=pid, seed=seed) for t in (status, infants)]


def run_emod(pid):
    d = pd.read_csv(EMOD / f'set{pid:03d}' / 'output' / 'ReportHIVByAgeAndGender.csv', skipinitialspace=True, keep_default_na=False) # Keep the status "None"
    d.columns = d.columns.str.strip()
    d = d.rename(columns={'IP_Key:InterventionStatus': 'InterventionStatus'})
    status, infants = tables(d)
    return [t.assign(model='emod', pid=pid, seed=0) for t in (status, infants)]


def reduce(cores):
    with Pool(cores, maxtasksperchild=1) as pool: # The variants patch classes, so one run per process
        res = pool.map(run_v4, [(p, s, v) for p in SETS for s in SEEDS for v in VARIANTS], chunksize=1) + pool.map(run_emod, SETS)
    for k, name in enumerate(['status', 'infants']):
        pd.concat([r[k] for r in res], ignore_index=True).to_csv(here / 'outputs' / f'pmtct_{name}.csv.gz', index=False)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from targetfit import style, c_run, c_best, c_muted
    looks = {'emod': dict(color=c_best, label='EMOD (3 sets)'), 'v4_as_is': dict(color=c_run, label='EMOD′ (3 sets × 3 seeds)'),
             'v4_unblocked': dict(color=c_run, ls='--', label='EMOD′, Option B not blocked'), 'v4_persist': dict(color=c_run, ls=':', label='EMOD′, PMTCT not reset')}
    st = pd.read_csv(here / 'outputs' / 'pmtct_status.csv.gz', keep_default_na=False)
    inf = pd.read_csv(here / 'outputs' / 'pmtct_infants.csv.gz')
    pos = st[st.HasHIV == 1]
    tot = pos.groupby(['model', 'pid', 'seed', 'year']).n.sum().rename('tot')
    share = pos.groupby(['model', 'pid', 'seed', 'year', 'InterventionStatus']).n.sum().reset_index().merge(tot.reset_index(), on=['model', 'pid', 'seed', 'year'])
    share['s'] = share.n/share.tot
    states = [s for s in share.InterventionStatus.unique() if share[share.InterventionStatus == s].s.max() > 0.005]
    fig, axs = plt.subplots(3, 4, figsize=(14, 8))
    for ax, s in zip(axs.flat, sorted(states)):
        for m in ('emod', 'v4_as_is'):
            g = share[(share.model == m) & (share.InterventionStatus == s)].groupby('year').s
            ax.plot(g.mean().index, g.mean(), lw=1.5, **looks[m])
        style(ax, s, 'Share of PLHIV')
    axs.flat[0].legend(fontsize=7, frameon=False)
    for ax in axs.flat[len(states):]:
        ax.set_visible(False)
    fig.text(0.01, -0.02, 'Mean over runs. InterventionStatus of people with HIV, all ages.', fontsize=7, color=c_muted)
    fig.savefig(here / 'figures' / 'intervention_status_compared.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    for m, lk in looks.items():
        g = inf[inf.model == m].groupby('year')
        axs[0].plot(g.new_0_4.mean().index, g.new_0_4.mean(), lw=1.5, **lk)
        p = g.infected_0_4.sum()/g.pop_0_4.sum()
        axs[1].plot(p.index, p, lw=1.5, **lk)
    style(axs[0], 'New infections at ages 0-4 a year (agents, mean of runs)')
    style(axs[1], 'HIV prevalence, ages 0-4')
    axs[0].legend(fontsize=7, frameon=False)
    fig.savefig(here / 'figures' / 'pmtct_infants.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=27)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
