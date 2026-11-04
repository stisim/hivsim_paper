"""
Compare HIVsim Malawi runs with the targets in results/targets.csv: measure each target in a run, score runs,
and plot the data against the runs of an experiment.

    python targetfit.py ../experiments/01_coverage_check   # Plots and coverage table from outputs/results.jsonl

Each target row is identified by its row number in targets.csv (`tid`). A run's results are a dict of tid -> value.
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

here = Path(__file__).parent
targetfile = here.parent / 'results' / 'targets.csv'

# Colours (reference palette: categorical slots 1 and 2, blue ramp, text inks)
c_run, c_runs, c_band, c_best, c_data, c_muted = '#2a78d6', '#9ec5f4', '#cde2fb', '#eb6834', '#0b0b0b', '#52514e'

stocks = ['prevalence', 'plhiv', 'on_art', 'population', 'diagnosed', 'on_ART_among_diagnosed', 'suppressed_among_on_ART', 'art_coverage', 'vls_among_plhiv']
flows = ['new_infections', 'aids_deaths', 'incidence']


def load_targets():
    t = pd.read_csv(targetfile)
    t.index.name = 'tid'
    return t


class TargetRecorder(ss.Analyzer):
    """ Measure the stock targets (and the susceptible people that incidence needs) at their times, from the agents """

    def __init__(self, targets=None, **kwargs):
        super().__init__(**kwargs)
        self.targets = load_targets() if targets is None else targets
        rows = self.targets[self.targets.target.isin(stocks + ['incidence'])]
        self.times = np.unique(rows.year.values)
        self.done = set()
        self.values = {}

    def mask(self, row):
        sim, hiv = self.sim, self.sim.diseases.hiv
        ppl = sim.people
        m = ppl.alive & (ppl.age >= row.age_lo) & (ppl.age < row.age_hi)
        if row.sex == 'f':
            m = m & ppl.female
        elif row.sex == 'm':
            m = m & ~ppl.female
        if row.group == 'fsw':
            m = m & sim.networks.structuredsexual.fsw
        return m

    def measure(self, tid, row):
        hiv, art = self.sim.diseases.hiv, self.sim.interventions.art
        m = self.mask(row)
        n = lambda x: float(np.count_nonzero(x)) * self.sim.pars.pop_scale
        share = lambda num, den: float(np.count_nonzero(num & den)) / max(np.count_nonzero(den), 1)
        inf = m & hiv.infected
        if row.target == 'prevalence': return share(hiv.infected, m)
        if row.target == 'plhiv': return n(inf)
        if row.target == 'on_art': return n(m & hiv.on_art)
        if row.target == 'population': return n(m)
        if row.target == 'diagnosed': return share(hiv.diagnosed, inf)
        if row.target == 'on_ART_among_diagnosed': return share(hiv.on_art, inf & hiv.diagnosed)
        if row.target == 'suppressed_among_on_ART': return share(art.on_effective_art, inf & hiv.on_art)
        if row.target == 'art_coverage': return share(hiv.on_art, inf)
        if row.target == 'vls_among_plhiv': return share(art.on_effective_art, inf)
        if row.target == 'incidence': return n(m & ~hiv.infected) # Susceptible people, the denominator
        raise ValueError(row.target)

    def step(self):
        year = self.t.now('year')
        for time in self.times:
            if time not in self.done and year >= time:
                self.done.add(time)
                for tid, row in self.targets[self.targets.year == time].iterrows():
                    if row.target in stocks or row.target == 'incidence':
                        self.values[tid] = self.measure(tid, row)


def measure(sim, targets=None):
    """ The model's value for each target: stocks from TargetRecorder, flows from the results """
    targets = load_targets() if targets is None else targets
    rec = sim.analyzers.targetrecorder
    years = np.asarray(sim.t.yearvec)
    res = sim.results.hiv
    out = {}
    for tid, row in targets.iterrows():
        if row.target in stocks:
            out[tid] = rec.values.get(tid, np.nan)
        elif row.target in ('new_infections', 'aids_deaths'):
            key = 'new_infections' if row.target == 'new_infections' else 'new_deaths'
            y0 = np.floor(row.year)
            out[tid] = float(np.asarray(res[key])[(years >= y0) & (years < y0 + 1)].sum())
        elif row.target == 'incidence': # New infections in the year to the survey, per susceptible person at the survey
            bins = [(lo, lo + 5) for lo in range(int(row.age_lo), int(row.age_hi), 5)]
            inyear = (years > row.year - 1) & (years <= row.year)
            new = sum(np.asarray(res[f'new_infections_{row.sex}_{lo}_{hi}'])[inyear].sum() for lo, hi in bins)
            out[tid] = float(new / rec.values[tid]) if rec.values.get(tid) else np.nan
    return out


def score(values, targets=None, roles=('primary', 'secondary')):
    """
    A summary of how far one run is from the data, for ranking runs (not a likelihood): the mean squared distance
    in units of each target's uncertainty, from its 95% interval where there is one, else 20% of the value
    """
    targets = load_targets() if targets is None else targets
    t = targets[targets.role.isin(roles)]
    model = np.array([values.get(tid, np.nan) for tid in t.index], dtype=float)
    sd = np.where(t.upper.notna() & t.lower.notna(), (t.upper - t.lower) / 3.92, 0.2 * t.value.abs())
    z = (model - t.value.values) / np.maximum(sd, 1e-12)
    return float(np.nanmean(z**2))


def load_runs(expdir):
    """ The runs of an experiment: a DataFrame of parameters and a DataFrame of target values (tid columns), one row per run """
    rows = [json.loads(line) for line in (Path(expdir) / 'outputs' / 'results.jsonl').read_text().splitlines() if line]
    runs = sorted([r for r in rows if r.get('kind', 'sample') == 'sample' and 'values' in r], key=lambda r: r['i'])
    pars = pd.DataFrame([r['pars'] for r in runs])
    vals = pd.DataFrame([{int(k): v for k, v in r['values'].items()} for r in runs])
    return pars, vals, rows


def coverage_table(vals, targets):
    """ For each target and source: how many data points lie inside the range of the runs, and inside their 5-95% range """
    lo, hi = vals.min(), vals.max()
    q5, q95 = vals.quantile(0.05), vals.quantile(0.95)
    t = targets.loc[vals.columns].copy()
    t['in_range'] = (t.value >= lo) & (t.value <= hi)
    t['in_90'] = (t.value >= q5) & (t.value <= q95)
    t['below'] = t.value < lo
    g = t.groupby(['role', 'target', 'source'], sort=False)
    return pd.DataFrame(dict(points=g.size(), in_range=g.in_range.sum(), in_90=g.in_90.sum(), data_below_all_runs=g.below.sum()))


def trajectory_check(vals, targets, factor=1.5, start=2000, art_tol=0.1):
    """
    Whole trajectories, not points, from `start` (Robyn: the later epidemic matters most): for each UNAIDS count,
    the share of runs within `factor` of the data in every year; for UNAIDS ART coverage, within `art_tol` (absolute)
    every year; and the shape of the epidemic curve, new infections in 2020 over their 1990-95 peak
    """
    rows = []
    for target in ['plhiv', 'new_infections', 'aids_deaths']:
        t = targets[(targets.target == target) & (targets.year >= start)].sort_values('year')
        ratio = vals[t.index].div(t.value.values, axis=1)
        within = ((ratio >= 1/factor) & (ratio <= factor)).all(axis=1)
        rows.append(dict(check=f'{target}: within x{factor} of the data every year from {start}', data=np.nan, runs_median=np.nan, share_of_runs=within.mean()))
    t = targets[(targets.target == 'art_coverage') & (targets.source == 'UNAIDS 2025')]
    if len(t) and set(t.index) <= set(vals.columns):
        within = (vals[t.index] - t.value.values).abs().le(art_tol).all(axis=1)
        rows.append(dict(check=f'art_coverage: within {art_tol:.0%} of the data every year', data=np.nan, runs_median=np.nan, share_of_runs=within.mean()))
    t = targets[targets.target == 'new_infections'].sort_values('year')
    early, late = t[t.year < 1996], t[np.floor(t.year) == 2020]
    data = late.value.item() / early.value.max()
    runs = vals[late.index[0]] / vals[early.index].max(axis=1)
    rows.append(dict(check='new infections 2020 / peak 1990-95', data=data, runs_median=runs.median(), share_of_runs=((runs > data/factor) & (runs < data*factor)).mean()))
    return pd.DataFrame(rows).set_index('check')


def style(ax, title=None, ylabel=None):
    ax.spines[['top', 'right']].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color(c_muted)
    ax.tick_params(colors=c_muted, labelsize=8)
    ax.grid(axis='y', color='#e8e7e3', lw=0.6)
    ax.set_axisbelow(True)
    if title: ax.set_title(title, fontsize=9, color=c_data, loc='left')
    if ylabel: ax.set_ylabel(ylabel, fontsize=8, color=c_muted)


def draw(ax, x, runs, best, data, lower=None, upper=None, xdata=None):
    """ Runs (thin lines), their 5-95% band, the best run, and the data with its interval """
    runs = np.asarray(runs, dtype=float)
    ax.fill_between(x, np.nanquantile(runs, 0.05, axis=0), np.nanquantile(runs, 0.95, axis=0), color=c_band, lw=0, label='Runs, 5-95%')
    for r in runs[:200]:
        ax.plot(x, r, color=c_runs, lw=0.4, alpha=0.5)
    ax.plot(x, best, color=c_best, lw=2, label='Best run')
    xd = x if xdata is None else xdata
    yerr = None
    if lower is not None and np.isfinite(lower).any():
        yerr = np.vstack([np.clip(data - lower, 0, None), np.clip(upper - data, 0, None)])
        yerr = np.nan_to_num(yerr)
    ax.errorbar(xd, data, yerr=yerr, fmt='o', ms=4, color=c_data, ecolor=c_data, elinewidth=1, capsize=2, label='Data', zorder=5)


def plot(expdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    expdir = Path(expdir)
    figdir = expdir / 'figures'
    figdir.mkdir(exist_ok=True)
    targets = load_targets()
    pars, vals, _ = load_runs(expdir)
    scores = np.array([score(dict(row.dropna()), targets) for _, row in vals.iterrows()])
    ib = int(np.nanargmin(scores))
    best = vals.iloc[ib]
    nruns = len(vals)

    # Prevalence by age, one panel per survey and sex
    prev = targets[(targets.target == 'prevalence') & (targets.group == 'all')]
    surveys = list(dict.fromkeys(prev.sort_values('year').source))
    fig, axes = plt.subplots(2, len(surveys), figsize=(3.2*len(surveys), 6), sharey=True, constrained_layout=True)
    for j, src in enumerate(surveys):
        for i, sex in enumerate(['f', 'm']):
            t = prev[(prev.source == src) & (prev.sex == sex)].sort_values('age_lo')
            x = (t.age_lo + np.minimum(t.age_hi, 70)) / 2
            draw(axes[i, j], x.values, vals[t.index].values, best[t.index].values, t.value.values, t.lower.values, t.upper.values)
            style(axes[i, j], f'{src}, {"women" if sex == "f" else "men"}', 'HIV prevalence' if j == 0 else None)
            axes[i, j].set_xlabel('Age', fontsize=8, color=c_muted)
            axes[i, j].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    axes[0, 0].legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle(f'Prevalence by age and sex: data and {nruns} runs', fontsize=10, x=0.01, ha='left')
    fig.savefig(figdir / 'coverage_prevalence_by_age.png', dpi=150)
    plt.close(fig)

    # Time series: UNAIDS, incidence, FSW prevalence, numbers on ART by sex
    panels = [('plhiv', 'all', 'all', 'PLHIV (UNAIDS)'), ('new_infections', 'all', 'all', 'New infections (UNAIDS)'),
              ('aids_deaths', 'all', 'all', 'AIDS deaths (UNAIDS)'), ('prevalence', 'f', 'fsw', 'Prevalence, female sex workers'),
              ('incidence', 'f', 'all', 'Incidence 15-49, women (MPHIA)'), ('incidence', 'm', 'all', 'Incidence 15-49, men (MPHIA)'),
              ('on_art', 'f', 'all', 'On ART 15-64, women (check)'), ('on_art', 'm', 'all', 'On ART 15-64, men (check)')]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6), constrained_layout=True)
    for ax, (target, sex, group, title) in zip(axes.flat, panels):
        t = targets[(targets.target == target) & (targets.sex == sex) & (targets.group == group)].sort_values('year')
        if not len(t):
            ax.set_visible(False)
            continue
        draw(ax, t.year.values, vals[t.index].values, best[t.index].values, t.value.values, t.lower.values, t.upper.values)
        style(ax, title)
        if target in ('prevalence', 'incidence'):
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=1 if target == 'incidence' else 0))
        else:
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v/1e3:,.0f}k'))
    axes[0, 0].legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle(f'Epidemic over time: data and {nruns} runs', fontsize=10, x=0.01, ha='left')
    fig.savefig(figdir / 'coverage_time_series.png', dpi=150)
    plt.close(fig)

    # ART and suppression: UNAIDS ART coverage over time, and the MPHIA 2020-21 cascade by sex
    cas = targets[targets.target.isin(['art_coverage', 'vls_among_plhiv', 'diagnosed'])]
    cas = cas[cas.index.isin(vals.columns)]
    if len(cas):
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.3, 1]))
        t = cas[(cas.target == 'art_coverage') & (cas.source == 'UNAIDS 2025')].sort_values('year')
        if len(t):
            draw(axes[0], t.year.values, vals[t.index].values, best[t.index].values, t.value.values, t.lower.values, t.upper.values)
        style(axes[0], 'ART coverage, share of PLHIV on ART (UNAIDS, all ages)')
        axes[0].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
        axes[0].legend(fontsize=7, frameon=False, loc='upper left')
        m = cas[cas.source == 'MPHIA 2020-21'].sort_values(['target', 'sex'])
        names = {'diagnosed': 'Diagnosed', 'art_coverage': 'On ART', 'vls_among_plhiv': 'Suppressed'}
        ax = axes[1]
        for k, (tid, r) in enumerate(m.iterrows()):
            q = vals[tid].quantile([0.05, 0.95]).values
            ax.plot([k, k], q, color=c_band, lw=10, solid_capstyle='butt', label='Runs, 5-95%' if k == 0 else None)
            ax.plot(k, best[tid], 'o', color=c_best, ms=7, label='Best run' if k == 0 else None)
            err = None if np.isnan(r.lower) else [[r.value - r.lower], [r.upper - r.value]]
            ax.errorbar(k, r.value, yerr=err, fmt='o', ms=4, color=c_data, capsize=2, label='Data' if k == 0 else None, zorder=5)
        ax.set_xticks(range(len(m)), [f'{names[r.target]}\n{"women" if r.sex == "f" else "men"}' for _, r in m.iterrows()], fontsize=7)
        style(ax, 'Share of PLHIV 15+ (MPHIA 2020-21)')
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
        fig.suptitle(f'ART and viral suppression: data and {nruns} runs', fontsize=10, x=0.01, ha='left')
        fig.savefig(figdir / 'coverage_art_suppression.png', dpi=150)
        plt.close(fig)

    # Coverage table and best run
    table = coverage_table(vals, targets)
    table.to_csv(expdir / 'outputs' / 'coverage_table.csv')
    summary = dict(n_runs=nruns, best_run=int(ib), best_score=float(scores[ib]), median_score=float(np.nanmedian(scores)),
                   best_pars=pars.iloc[ib].to_dict())
    traj = trajectory_check(vals, targets)
    traj.to_csv(expdir / 'outputs' / 'trajectory_check.csv')
    print('\nWhole trajectories:\n' + traj.round(3).to_string())

    # Noise: spread across seeds of one parameter set, against the spread across parameter sets
    central = [r for r in load_runs(expdir)[2] if r.get('kind') == 'central' and 'values' in r]
    if len(central) > 1:
        cvals = pd.DataFrame([{int(k): v for k, v in r['values'].items()} for r in central])
        ratio = (cvals.std() / vals[cvals.columns].std()).replace([np.inf], np.nan)
        noise = targets.loc[ratio.index].assign(ratio=ratio).groupby(['target', 'source'], sort=False).ratio.median()
        noise.to_csv(expdir / 'outputs' / 'noise_ratio.csv')
        summary['noise_ratio_median'] = float(ratio.median())
        print('\nSpread across seeds / spread across parameter sets (median by target):\n' + noise.round(2).to_string())
    sc.savejson(expdir / 'outputs' / 'coverage_summary.json', summary)
    print(table.to_string())
    print(f'\nBest run {ib}: score {scores[ib]:.2f} (median {np.nanmedian(scores):.2f}); parameters {pars.iloc[ib].round(4).to_dict()}')
    return table, summary


if __name__ == '__main__':
    plot(sys.argv[1])
