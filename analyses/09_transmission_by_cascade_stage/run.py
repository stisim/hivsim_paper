"""
Analysis 09: transmissions by the source's cascade stage, EMOD against EMOD′ (condensed cascade). See README.md.
Needs Starsim at ebd52c1a on PYTHONPATH, and the EMOD runs (models/emod/reference/run_emod.py --transmission
--seeds 1 2 3 4 5 --out tr10x5 for the sets below).

    python run.py          # EMOD′ runs (unless their logs are in raw_results/), and both classified, to outputs/
    python run.py --plot   # figures/ and tables
"""
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
repo = here.parents[1]
for p in ('', 'models/emod/v4'):
    sys.path.insert(0, str(repo / p))
from common import transmission as tx # noqa: E402

SETS = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225]
SEEDS = [1, 2, 3, 4, 5]
EMOD = repo / 'raw_results' / 'emod' / 'tr10x5'
RAW = repo / 'raw_results' / 'emod_prime' / 'tr10x5' # Per-transmission logs (not committed)
YEARS = [2010, 2015, 2020, 2024]


def sexual(sim):
    return getattr(sim.diseases.hiv, 'sexual', (np.zeros(0, int), np.zeros(0, int)))


def run_prime(job):
    import harness, inputs
    from cascade_host import EmodHost
    pid, seed = job
    f = RAW / f'set{pid:03d}_s{seed}.csv.gz'
    if f.exists(): # Already run: classify again only
        return summarize(pd.read_csv(f), 'EMOD′', pid, seed)
    sim = harness.make_sim(pid, seed=inputs.load(pid).run_number*10 + seed, cascade=True, condensed=True, record=False,
                           analyzers=[tx.TransmissionLog(EmodHost, sexual)])
    sim.run()
    log = sim.analyzers.transmissionlog.df
    log.to_csv(f, index=False)
    return summarize(log, 'EMOD′', pid, seed)


def run_emod(job):
    pid, seed = job
    log = tx.from_emod(EMOD / f'set{pid:03d}_s{seed}' / 'output')
    out = summarize(log, 'EMOD', pid, seed)
    ip_on = log.art_status.isin(['OnEffectiveART', 'OnNonSuppressiveART'])
    out['check'] = dict(model='EMOD', pid=pid, seed=seed, n=len(log), art_ip_mismatch=int((log.on_art != ip_on).sum()),
                        on_art_ge6m_not_on_by_ip=int((log.on_art & (log.since_start >= tx.HALF) & ~ip_on).sum()),
                        **daniel_check(log, EMOD / f'set{pid:03d}_s{seed}' / 'output'))
    return out


def summarize(log, model, pid, seed):
    tag = dict(model=model, pid=pid, seed=seed)
    fine = tx.shares(log).reset_index().assign(**tag)
    coarse = tx.shares(log, coarse=True).reset_index().assign(**tag)
    return dict(fine=fine, coarse=coarse)


def daniel_check(log, folder):
    """ Daniel's 'diagnosed, not on ART' (mwi_transmission.Rmd, fun.diag.noart): transmissions after a diagnosis event, minus
    transmissions after a first ART start, each counted separately, with his report settings (events from 2016, diagnoses
    at ages 15+, transmissions from 2016) and with ours (full history). Negative when the second set isn't inside the first """
    ev = pd.read_csv(folder / 'ReportEventRecorder.csv', usecols=['Time', 'Year', 'Event_Name', 'Individual_ID', 'Age'])
    out = {}
    for tag, y0, a0 in (('daniel', 2016, 15), ('full', 0, 0)):
        e = ev[ev.Year >= y0]
        diag = e[(e.Event_Name == 'HIVNewlyDiagnosed') & (e.Age >= a0*365)].groupby('Individual_ID').Time.min()
        art = e[e.Event_Name == 'StartedART'].groupby('Individual_ID').Time.min()
        d = log[log.year >= max(y0, 2016)]
        after_dx = d.time > d.src.map(diag) # NaN compares False
        after_art = d.time > d.src.map(art)
        yr = np.floor(d.year)
        diff = after_dx.groupby(yr).sum() - after_art.groupby(yr).sum()
        out[f'{tag}_diag_noart_min'] = int(diff.min())
        out[f'{tag}_art_without_dx'] = int((after_art & ~after_dx).sum()) # In his n_art but not in his n_diag
        out[f'{tag}_n'] = len(d)
    return out


def reduce(cores):
    RAW.mkdir(parents=True, exist_ok=True)
    jobs = [(p, s) for p in SETS for s in SEEDS]
    with Pool(cores) as pool:
        res = pool.map(run_prime, jobs, chunksize=1) + pool.map(run_emod, jobs, chunksize=1)
    out = here / 'outputs'
    out.mkdir(exist_ok=True)
    for k in ('fine', 'coarse'):
        pd.concat([r[k] for r in res], ignore_index=True).to_csv(out / f'shares_{k}.csv.gz', index=False)
    pd.DataFrame([r['check'] for r in res if 'check' in r]).to_csv(out / 'emod_checks.csv', index=False)
    return


def compare(kind):
    """ Per year and category: each model's median and 5-95% over runs, and the difference of means with its SE """
    s = pd.read_csv(here / 'outputs' / f'shares_{kind}.csv.gz')
    cats = tx.CATEGORIES if kind == 'fine' else tx.COARSE_ORDER
    long = s.melt(id_vars=['year', 'model', 'pid', 'seed', 'n'], value_vars=cats, var_name='category', value_name='share')
    g = long.groupby(['category', 'year', 'model']).share
    t = pd.DataFrame(dict(median=g.median(), lo=g.quantile(0.05), hi=g.quantile(0.95), mean=g.mean(), var=g.var(), runs=g.size())).unstack('model')
    d = pd.DataFrame(index=t.index)
    for m in ('EMOD', 'EMOD′'):
        for c in ('median', 'lo', 'hi'):
            d[f'{m} {c}'] = t[(c, m)]
    d['diff'] = t[('mean', 'EMOD′')] - t[('mean', 'EMOD')]
    d['se'] = np.sqrt(t[('var', 'EMOD′')]/t[('runs', 'EMOD′')] + t[('var', 'EMOD')]/t[('runs', 'EMOD')])
    d['z'] = d['diff']/d.se.where(d.se > 0)
    return d.reset_index()


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(repo / 'common'))
    from targetfit import style, c_run, c_best, c_muted
    figdir = here / 'figures'
    figdir.mkdir(exist_ok=True)
    looks = {'EMOD': dict(color=c_best, label='EMOD'), 'EMOD′': dict(color=c_run, label='EMOD′')}
    labels = dict(undiagnosed='Undiagnosed', diagnosed_naive='Diagnosed, never on ART', art_lt6m='On ART < 6 months', art_supp='On ART ≥ 6 months, suppressive arm',
                  art_nonsupp='On ART ≥ 6 months, non-suppressive arm', interrupted='Interrupted (off ART after starting)',
                  art_first_lt6m='On ART < 6 months, first start', art_reinit_lt6m='On ART < 6 months, restarted', off_first='Off ART, first interruption',
                  off_later_lt6m='Off ART, later interruption, < 6 months', off_later_ge6m='Off ART, later interruption, ≥ 6 months')
    for kind in ('coarse', 'fine'):
        c = compare(kind)
        c.to_csv(here / 'outputs' / f'compare_{kind}.csv', index=False, float_format='%.4f')
        c[c.year.isin(YEARS)].to_csv(here / 'outputs' / f'compare_{kind}_years.csv', index=False, float_format='%.4f')
        cats = tx.COARSE_ORDER if kind == 'coarse' else tx.CATEGORIES
        ncol = 3
        fig, axs = plt.subplots(int(np.ceil(len(cats)/ncol)), ncol, figsize=(13, 3.2*np.ceil(len(cats)/ncol)), squeeze=False)
        for ax, cat in zip(axs.flat, cats):
            d = c[(c.category == cat) & (c.year >= 1995) & (c.year <= 2024)]
            for m, lk in looks.items():
                ax.fill_between(d.year, d[f'{m} lo'], d[f'{m} hi'], color=lk['color'], alpha=0.15, lw=0)
                ax.plot(d.year, d[f'{m} median'], lw=1.5, **lk)
            if cat == 'interrupted':
                ax.axhspan(0.084, 0.201, color=c_muted, alpha=0.12, lw=0)
                ax.text(1996, 0.21, 'MIHPSA Malawi range across models', fontsize=7, color=c_muted)
            style(ax, labels[cat])
        for ax in axs.flat[len(cats):]:
            ax.set_visible(False)
        axs[0, 0].legend(fontsize=7, frameon=False)
        fig.text(0.01, -0.02, 'Share of the year\'s sexual transmissions by the source\'s cascade stage. Lines: median over 50 runs (10 parameter sets × 5 seeds); '
                 'bands: 5-95%. 2025 (part year) left out.', fontsize=7, color=c_muted, va='top')
        fig.tight_layout()
        fig.savefig(figdir / f'shares_{kind}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--cores', type=int, default=50)
    a = ap.parse_args()
    plot() if a.plot else reduce(a.cores)
