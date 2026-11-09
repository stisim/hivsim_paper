"""
Analysis 12 addendum: why hybrid 7's commercial relationships drift above EMOD′'s. See SUMMARY.md (addendum).

    python commercial.py   # outputs/commercial.csv, commercial_rates.csv; figures/commercial.png

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


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--plot', action='store_true')
    a = p.parse_args()
    run()
