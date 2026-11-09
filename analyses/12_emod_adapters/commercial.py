"""
Analysis 12 addendum: why hybrid 7's commercial relationships drift above EMOD′'s. See SUMMARY.md (addendum).

    python commercial.py         # outputs/commercial.csv, commercial_rates.csv (the *_before.csv files are from before
                                 # `EmodNetworkAt(burnin_risk=True)`, i.e. the warm-up run with the start's sex work)
    python commercial.py --plot  # figures/commercial.png, outputs/commercial_summary.csv

Records, each mid-year and node: people at Risk HIGH past debut by sex (FSW, clients), commercial queue entries by sex,
commercial relationships formed in the year, and the frozen commercial entry rates per person (by sex).
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
out, figs = here / 'outputs', here / 'figures'
SEEDS = [1, 2, 3, 4, 5, 6]
C = 3 # COMMERCIAL


def make_track():
    import starsim as ss

    class Track(ss.Analyzer):
        def __init__(self, **kwargs):
            super().__init__(name='ctrack', **kwargs)
            self.rows, self.rates, self.formed = [], None, {}
            return

        def step(self):
            net, ppl = self.sim.networks.emodnet, self.sim.people
            year = self.t.now('year')
            if self.rates is None and net.burnin_left <= 0:
                self.rates = net.rates[:, C, :, 0, 1].copy() # Node × sex, per day, flag holders (all HIGH hold it)
                self.rate_year = year
            if round((year % 1)*12) != 6:
                return
            au = ppl.auids
            past = ppl.age.raw[au] >= net.debut.raw[au]
            high = au[past & (ppl.risk.raw[au] == 2)]
            Q = net.queue
            qc = Q.uid[Q.type == C]
            for n in range(3):
                row = dict(year=int(year), node=n)
                for s, lab in ((1, 'f'), (0, 'm')):
                    h = high[(net.node.raw[high] == n) & (ppl.female.raw[high] == s)]
                    row[f'high_{lab}'] = len(h)
                    row[f'age_high_{lab}'] = ppl.age.raw[h].mean() if len(h) else np.nan
                    qq = qc[(net.node.raw[qc] == n) & (ppl.female.raw[qc] == s)]
                    row[f'queue_{lab}'] = len(qq)
                    row[f'qpeople_{lab}'] = len(np.unique(qq))
                row['n1549'] = int(((ppl.age.raw[au] >= 15) & (ppl.age.raw[au] < 50) & (net.node.raw[au] == n)).sum())
                self.rows.append(row)
            return

    return Track()


def run_one(job):
    name, seed = job
    import harness, adapters
    tr = make_track()
    if name == 'home':
        sim = harness.make_sim(pid=0, seed=seed, hiv=False, network=True, stop=2025, analyzers=[tr])
    else:
        sim = adapters.make_sim(pid=0, seed=seed, network='emod', hiv=False, stop=2025, analyzers=[tr])
    sim.run()
    tr = sim.analyzers.ctrack # The sim's copy
    net = sim.networks.emodnet
    df = pd.DataFrame(tr.rows).assign(model=name, seed=seed)
    st = net.starts
    st = st[st.type == C].assign(year=lambda d: sim.t.yearvec[d.ti.values].astype(int))
    formed = st.groupby(['year', 'node']).size().rename('formed')
    df = df.merge(formed.reset_index(), on=['year', 'node'], how='left').fillna({'formed': 0})
    r = net.rates[:, C, :, 0, 1] if tr.rates is None else tr.rates
    rates = pd.DataFrame([dict(model=name, seed=seed, node=n, rate_f=r[n, 1], rate_m=r[n, 0], rate_year=tr.rate_year) for n in range(3)])
    return df, rates


def run():
    jobs = [(m, s) for m in ('home', 'hyb7') for s in SEEDS]
    with Pool(len(jobs)) as pool:
        res = pool.map(run_one, jobs)
    pd.concat([r[0] for r in res]).to_csv(out / 'commercial.csv', index=False)
    pd.concat([r[1] for r in res]).to_csv(out / 'commercial_rates.csv', index=False)
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = []
    for tag, lab in (('_before', 'before'), ('', 'after')):
        d = pd.read_csv(out / f'commercial{tag}.csv')
        d = d.groupby(['model', 'seed', 'year'])[['high_f', 'high_m', 'formed', 'n1549']].sum().reset_index()
        d = d.groupby(['model', 'year']).mean().reset_index()
        d['model'] = np.where(d.model == 'home', 'EMOD′', f'7 ({lab})')
        rows.append(d)
    d = pd.concat(rows).drop_duplicates(['model', 'year'])
    d = d.assign(per_fsw=d.formed/d.high_f, per1549=d.formed/d.n1549, fsw_share=100*d.high_f/d.n1549, clients_per_fsw=d.high_m/d.high_f)
    home = d[d.model == 'EMOD′'].set_index('year')
    for m in ('7 (before)', '7 (after)'):
        h = d[d.model == m].set_index('year')
        d.loc[d.model == m, 'per1549_ratio'] = (h.per1549/home.per1549).reindex(h.index).values
    d.round(4).to_csv(out / 'commercial_summary.csv', index=False)
    col = {'EMOD′': '#333333', '7 (before)': '#f4a582', '7 (after)': '#d62728'}
    fig, axs = plt.subplots(1, 4, figsize=(16, 3.6))
    for m, g in d[d.year >= 1975].groupby('model'):
        for ax, k in zip(axs, ('per_fsw', 'per1549', 'fsw_share', 'clients_per_fsw')):
            ax.plot(g.year, g[k], color=col[m], lw=1.5, label=m)
    for ax, t in zip(axs, ('Commercial relationships formed per FSW per year', 'Commercial formed per 15–49 year old per year', 'FSW (Risk HIGH women) % of 15–49', 'Clients per FSW (Risk HIGH men / women)')):
        ax.set_title(t, fontsize=9)
    axs[0].legend(frameon=False)
    fig.suptitle('Commercial relationships: EMOD′ vs hybrid 7 (no HIV), warm-up with the start\'s sex work (before) or with the sex work at the end of EMOD\'s burn-in, 1975.5 (after); mean of 6 seeds', fontsize=10)
    fig.tight_layout()
    fig.savefig(figs / 'commercial.png', dpi=120)
    return


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--plot', action='store_true')
    a = p.parse_args()
    plot() if a.plot else run()
