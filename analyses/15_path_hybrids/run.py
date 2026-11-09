"""
Analysis 15: the five hybrids of the first path through the runner. See README.md.

    python run.py runs     # 5 hybrids x 5 seeds, 1985-2025, into raw_results/swaps/
    python run.py check    # seeds 1-2 rerun in-process: stored table reproduced, NaNs, outputs' counts against the disease's own
"""
import sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path.insert(0, str(repo))
from common import swaps # noqa: E402

OUT = here / 'outputs'
CONFIG = {'from': 'hivsim', 'to': 'emod', 'slots': ['C', 'V', 'RN', 'D'], 'pairs': {'hivsim': [0], 'emod': [0]}, 'seeds': 5, 'stop': 2025}
SHARES = dict(share_deaths_on_art='hiv_deaths', art_coverage='plhiv', incidence_15_49_f='py_15_49_f', incidence_15_49_m='py_15_49_m') # Output: its denominator
# EmodHIV's new_infections are sexual only, so they are compared with the log's sexual transmissions (tx_n); sti.HIV's are every route
RESULTS = dict(new_infections=('new_infections', 'new_infections'), hiv_deaths=('new_deaths', 'hiv_deaths'), infected=('infected', 'n_infected')) # Outputs' count: (sti.HIV's result, EmodHIV's)


def runs():
    rows = swaps.run(CONFIG, cores=25)
    print(f'{len(rows)} runs, {sum(bool(r["error"]) for r in rows)} errors')
    OUT.mkdir(exist_ok=True)
    man = pd.read_csv(swaps.runner.OUT / 'manifest.csv')
    ids = {swaps.runner.hybrid_id(h, CONFIG, dict(swaps.runner.code_version(), **swaps.runner.pins())) for h in swaps.hybrids(CONFIG)}
    man = man[man.hybrid.isin(ids)].drop_duplicates(['hybrid', 'cell'], keep='last')
    t = man.groupby('name', sort=False).agg(runs=('cell', 'size'), errors=('error', lambda e: int(e.notna().sum() - (e == '').sum())), s_min=('seconds', 'min'), s_max=('seconds', 'max'))
    t.to_csv(OUT / 'run_times.csv')
    print(t.to_string())


def nans(y):
    """ Outputs that are NaN although their denominator is positive (or that are NaN with no denominator to explain them) """
    bad = []
    for c in y.columns:
        if c in ('year',) or not y[c].isna().any():
            continue
        d = SHARES.get(c, 'tx_n' if c.startswith('tx_') else None)
        if d is None or (y.loc[y[c].isna(), d] > 0).any():
            bad.append(c)
    return bad


def check_one(job):
    """ One hybrid and cell, rerun in-process: the stored table reproduced, NaNs, and the outputs' counts by year against the disease's results """
    from common import outputs
    i, choice, cell = job
    folder = swaps.runner.OUT / swaps.runner.hybrid_id(choice, CONFIG, dict(swaps.runner.code_version(), **swaps.runner.pins()))
    sim = swaps.runner.start(choice, cell, CONFIG, folder)
    sim.run()
    o = sim.analyzers.outputs
    y = outputs.yearly(o)
    stored = pd.read_csv(folder / f'{cell["cell"]}.csv.gz')
    same = stored.shape == y.shape and np.allclose(stored.to_numpy(float), y.to_numpy(float), equal_nan=True, rtol=1e-12, atol=0)
    hiv = sim.diseases.hiv
    emod = hasattr(hiv, 'sexual')
    years = np.asarray(sim.t.yearvec, float)
    yr = np.floor(years + 1e-9).astype(int)
    s = o.strata.groupby('year')[['new_infections', 'hiv_deaths', 'infected']].sum()
    rows = []
    for k, (rs, re) in RESULTS.items():
        res = np.asarray(hiv.results[re if emod else rs], float)/float(sim.pop_scale or 1) # In agents
        out = s[k].to_numpy()
        if k == 'infected': # Stock: at the year's last step (the outputs count at the end of the step)
            last = pd.Series(np.arange(len(yr))).groupby(yr).max()
            a, b = out[last.to_numpy()], res[last.to_numpy()]
            b1 = res[np.minimum(last.to_numpy() + 1, len(res) - 1)] # The result recorded one step later
        else:
            a, b = pd.Series(out).groupby(yr).sum().to_numpy(), pd.Series(res).groupby(yr).sum().to_numpy()
            if emod and k == 'new_infections':
                a = y.set_index('year').tx_n.reindex(np.unique(yr)).fillna(0).to_numpy() if 'tx_n' in y else np.zeros(len(b))
            b1 = b
        for j, year in enumerate(np.unique(yr)):
            rows.append(dict(hybrid=i, name=swaps.name(choice), cell=cell['cell'], measure=k, year=year, outputs=a[j], disease=b[j], disease_next_step=b1[j]))
    summary = dict(hybrid=i, name=swaps.name(choice), cell=cell['cell'], stored_reproduced=same, nan_outputs=';'.join(nans(y)),
                   new_infections=int(y.new_infections.sum()), tx_n=int(y.get('tx_n', pd.Series([0])).sum()), hiv_deaths=int(y.hiv_deaths.sum()),
                   prev_2020=float(s.infected.loc[s.index[s.index < 2021].max()]/o.strata.groupby('year')['pop'].sum().loc[s.index[s.index < 2021].max()]))
    return summary, rows, y.assign(hybrid=i, name=swaps.name(choice), cell=cell['cell'])


def check():
    OUT.mkdir(exist_ok=True)
    jobs = [(i, h, c) for i, h in enumerate(swaps.hybrids(CONFIG)) for c in swaps.cells(CONFIG)[:2]]
    with Pool(len(jobs)) as pool:
        res = pool.map(check_one, jobs)
    summ = pd.DataFrame([r[0] for r in res])
    counts = pd.DataFrame([x for r in res for x in r[1]])
    counts['ratio'] = counts.outputs/counts.disease
    pd.concat([r[2] for r in res]).to_csv(OUT / 'yearly.csv', index=False)
    counts.to_csv(OUT / 'counts_vs_disease.csv', index=False)
    agg = counts[counts.disease > 0].groupby(['hybrid', 'measure']).ratio.agg(['min', 'median', 'max']).unstack()
    summ.to_csv(OUT / 'check.csv', index=False)
    print(summ.to_string())
    print(agg.round(4).to_string())
    agg.round(4).to_csv(OUT / 'counts_ratio.csv')
    d = swaps.load(CONFIG)
    took = d[d.year == 2020].assign(plhiv_2020=lambda x: x.plhiv).groupby('name', sort=False).agg(seeds=('cell', 'size'), took_off=('plhiv', lambda v: int((v > 0).sum())),
           plhiv_2020_median=('plhiv', 'median'), new_inf_2020_median=('new_infections', 'median'), art_cov_2020_median=('art_coverage', 'median'), tx_undiag_2020_median=('tx_undiagnosed', 'median'))
    took.to_csv(OUT / 'epidemic_2020.csv')
    print(took.to_string())


if __name__ == '__main__':
    dict(runs=runs, check=check)[sys.argv[1]]()
