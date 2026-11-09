"""
Analysis 12: EMOD′ components on HIVsim's slots, against EMOD′ itself. See README.md.

    python run.py          # 5 models × 6 seeds, to outputs/
    python run.py --plot   # figures/

Needs STIsim v4 (branch rs-v4-art) on PYTHONPATH and Starsim 4.0.0rc3.
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

PID = 0
SEEDS = [1, 2, 3, 4, 5, 6]
MODELS = { # name: (what, HIV)
    'home_hiv': ('EMOD′', True),
    'home_net': ('EMOD′, no HIV', False),
    'hyb6':     ('6: EmodHIV on HIVsim (StructuredSexual, STIsim demography)', True),
    'hyb7':     ('7: EmodHIV, EmodNetwork, EmodRisk on STIsim demography', True),
    'hyb7_net': ('7, no HIV', False),
}
AGEGRP = [0, 15, 25, 35, 200]
AGELAB = ['<15', '15-24', '25-34', '35+']
NETAGE = [15, 25, 35, 50]
NETLAB = ['15-24', '25-34', '35-49']
SNAP_YEARS = [1990, 2000, 2010, 2020]
TYPES = ['T', 'I', 'M', 'C']
STOP = 2025


def make_track():
    import starsim as ss

    class Track(ss.Analyzer):
        """ Infections (start, age, end by HIV death or other death), yearly epidemic counts, and network snapshots """
        def __init__(self, **kwargs):
            super().__init__(name='track', **kwargs)
            self.inf, self.end, self.epi, self.net = [], [], [], []
            return

        def step(self):
            sim, ppl = self.sim, self.sim.people
            year = self.t.now('year')
            hiv = sim.diseases.get('hiv')
            au = ppl.auids
            age = ppl.emod_age.raw[au]/365
            if hiv is not None:
                new = hiv.new_uids
                if len(new):
                    seeded = hiv.year_infected.raw[new] < year - 1e-6 # Fast-forwarded at seeding
                    self.inf.append(pd.DataFrame(dict(uid=new, t0=year, age=ppl.emod_age.raw[new]/365, seeded=seeded)))
                died = au[ppl.ti_dead.raw[au] == self.ti]
                if len(died):
                    hd = np.isin(died, hiv.dead_uids)
                    self.end.append(pd.DataFrame(dict(uid=died, t1=year, hiv=hd)))
            if round((year % 1)*12) == 6: # Mid-year
                female = ppl.female.raw[au]
                a1549 = (age >= 15) & (age < 50)
                row = dict(year=int(year), n1549=a1549.sum(), high1549=(ppl.risk.raw[au][a1549] == 2).sum())
                if hiv is not None:
                    row['inf1549'] = hiv.infected.raw[au][a1549].sum()
                self.epi.append(row)
                net = sim.networks.get('emodnet')
                if net is not None and int(year) in SNAP_YEARS:
                    n = net.n.raw[au]
                    g = np.digitize(age, NETAGE) - 1
                    for s in (0, 1):
                        for k in range(len(NETLAB)):
                            m = (female == s) & (g == k)
                            for t, tk in enumerate(TYPES):
                                self.net.append(dict(year=int(year), female=s, age=NETLAB[k], type=tk, n=m.sum(), rels=n[m, t].sum(), any=(n[m, t] > 0).sum()))
            return

    return Track()


def survival(inf, end, stop=STOP):
    """ Kaplan-Meier survival from infection to HIV death by age group, censored at other deaths and the stop, at whole years """
    if not len(inf):
        return pd.DataFrame()
    df = inf[~inf.seeded].merge(end, on='uid', how='left')
    df['t'] = np.where(df.t1.isna(), stop, df.t1) - df.t0
    df['event'] = df.hiv.fillna(False).astype(bool)
    df['grp'] = pd.cut(df.age, AGEGRP, right=False, labels=AGELAB)
    rows = []
    for g, d in df.groupby('grp', observed=True):
        t, e = d.t.values, d.event.values
        times = np.unique(t[e])
        s, surv = 1.0, []
        for u in times:
            at_risk = (t >= u).sum()
            s *= 1 - (e & (t == u)).sum()/at_risk
            surv.append(s)
        grid = np.arange(0, 31)
        idx = np.searchsorted(times, grid, side='right') - 1
        S = np.where(idx >= 0, np.array(surv + [1.0])[idx], 1.0)
        rows += [dict(grp=g, year=int(y), S=v, n=len(d), at_risk=int((t >= y).sum())) for y, v in zip(grid, S)]
    return pd.DataFrame(rows)


def run_one(job):
    name, seed = job
    import harness
    import adapters
    hiv = MODELS[name][1]
    track = make_track()
    if name.startswith('home'):
        sim = harness.make_sim(pid=PID, seed=seed, hiv=hiv, network=True, stop=STOP, analyzers=[track])
    else:
        sim = adapters.make_sim(pid=PID, seed=seed, network='sti' if name == 'hyb6' else 'emod', hiv=hiv, stop=STOP, analyzers=[track])
    sim.run()
    tr = sim.analyzers.track
    tag = dict(model=name, seed=seed)
    res = dict(epi=pd.DataFrame(tr.epi).assign(**tag), net=pd.DataFrame(tr.net).assign(**tag))
    if hiv:
        inf = pd.concat(tr.inf) if tr.inf else pd.DataFrame(columns=['uid', 't0', 'age', 'seeded'])
        end = pd.concat(tr.end) if tr.end else pd.DataFrame(columns=['uid', 't1', 'hiv'])
        res['surv'] = survival(inf, end).assign(**tag)
        h = sim.diseases.hiv
        years = h.infections.year.astype(int) if len(h.infections) else pd.Series(dtype=int)
        res['epi']['new_inf'] = res['epi'].year.map(years.value_counts()).fillna(0).astype(int)
        res['epi']['n_seeds'] = int((h.infections.type == -1).sum()) if len(h.infections) else 0
    net = sim.networks.get('emodnet')
    if net is not None and len(net.starts):
        st = net.starts.assign(year=lambda d: (sim.t.yearvec[d.ti.values]).astype(int))
        st = st[st.year >= 1986]
        g = st.assign(gap=st.age_m - st.age_f).groupby(['year', 'type'])
        res['form'] = g.size().rename('new').reset_index().merge(g.gap.mean().rename('gap').reset_index()).assign(**tag)
        res['gap'] = st.assign(gap=st.age_m - st.age_f).groupby('type').gap.describe()[['count', 'mean', 'std', '25%', '50%', '75%']].reset_index().assign(**tag)
    return res


def run():
    out.mkdir(exist_ok=True)
    jobs = [(m, s) for m in MODELS for s in SEEDS]
    with Pool(len(jobs)) as pool:
        results = pool.map(run_one, jobs)
    for key in ('epi', 'net', 'surv', 'form', 'gap'):
        dfs = [r[key] for r in results if key in r and len(r[key])]
        if dfs:
            pd.concat(dfs, ignore_index=True).to_csv(out / f'{key}.csv', index=False)
    compare()
    return


def compare():
    """ Tables comparing each hybrid with EMOD′ """
    surv = pd.read_csv(out / 'surv.csv')
    s = surv[surv.year.isin([2, 5, 10, 15])].groupby(['model', 'grp', 'year']).S.agg(['mean', 'std', 'count']).reset_index()
    s['se'] = s['std']/np.sqrt(s['count'])
    home = s[s.model == 'home_hiv'].set_index(['grp', 'year'])
    rows = []
    for m in ('hyb6', 'hyb7'):
        h = s[s.model == m].set_index(['grp', 'year'])
        j = h.join(home, rsuffix='_home', how='inner')
        j['z'] = (j['mean'] - j['mean_home'])/np.sqrt(j.se**2 + j.se_home**2)
        rows.append(j.reset_index().assign(model=m))
    pd.concat(rows)[['model', 'grp', 'year', 'mean', 'se', 'mean_home', 'se_home', 'z']].round(4).to_csv(out / 'survival_vs_home.csv', index=False)

    net = pd.read_csv(out / 'net.csv')
    net = net[net.model.isin(['home_net', 'hyb7_net'])]
    g = net.assign(per=net.rels/net.n, share=net['any']/net.n).groupby(['model', 'year', 'female', 'age', 'type'])[['per', 'share']].agg(['mean', 'std']).reset_index()
    g.columns = ['_'.join(c).strip('_') for c in g.columns]
    w = g.pivot_table(index=['year', 'female', 'age', 'type'], columns='model', values=['per_mean', 'per_std', 'share_mean'])
    w.columns = [f'{a}_{b}' for a, b in w.columns]
    w['ratio'] = w['per_mean_hyb7_net']/w['per_mean_home_net']
    w['z'] = (w['per_mean_hyb7_net'] - w['per_mean_home_net'])/np.sqrt((w['per_std_hyb7_net']**2 + w['per_std_home_net']**2)/len(SEEDS))
    w.round(4).reset_index().to_csv(out / 'network_vs_home.csv', index=False)

    gap = pd.read_csv(out / 'gap.csv')
    gap = gap[gap.model.isin(['home_net', 'hyb7_net'])].groupby(['model', 'type'])[['mean', 'std', '50%']].mean().round(2)
    gap.to_csv(out / 'agegap_vs_home.csv')

    epi = pd.read_csv(out / 'epi.csv')
    form = pd.read_csv(out / 'form.csv').merge(epi[['model', 'seed', 'year', 'n1549']], on=['model', 'seed', 'year'])
    form['per1549'] = form.new/form.n1549
    f = form[form.model.isin(['home_net', 'hyb7_net']) & form.year.isin(SNAP_YEARS)].groupby(['model', 'year', 'type']).per1549.mean().unstack('model').round(4)
    f.to_csv(out / 'formation_vs_home.csv')
    e = epi.assign(prev=epi.inf1549/epi.n1549, high=epi.high1549/epi.n1549)
    e[e.year.isin([1990, 2000, 2010, 2020, 2025])].groupby(['model', 'year'])[['prev', 'high']].agg(['mean', 'std']).round(4).to_csv(out / 'epi_summary.csv')
    if 'n_seeds' in epi:
        epi.groupby(['model', 'seed']).agg(n_seeds=('n_seeds', 'first'), prev2025=('inf1549', 'last')).to_csv(out / 'seeds.csv')
    return


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figs.mkdir(exist_ok=True)
    col = {'home_hiv': '#333333', 'home_net': '#333333', 'hyb6': '#1f77b4', 'hyb7': '#d62728', 'hyb7_net': '#d62728'}

    surv = pd.read_csv(out / 'surv.csv')
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
    for ax, g in zip(axs, AGELAB):
        for m in ('home_hiv', 'hyb6', 'hyb7'):
            d = surv[(surv.model == m) & (surv.grp == g)]
            for _, r in d.groupby('seed'):
                ax.plot(r.year, r.S, color=col[m], alpha=0.25, lw=0.8)
            mean = d.groupby('year').S.mean()
            ax.plot(mean.index, mean.values, color=col[m], lw=2, label=MODELS[m][0].split(':')[0])
        ax.set_title(f'Infected at {g}')
        ax.set_xlabel('Years since infection')
    axs[0].set_ylabel('Not yet dead of HIV (Kaplan–Meier)')
    axs[0].legend(frameon=False)
    fig.suptitle('EmodHIV without ART: survival from infection, EMOD′ vs hybrids (6 seeds each)')
    fig.tight_layout()
    fig.savefig(figs / 'survival.png', dpi=120)

    net = pd.read_csv(out / 'net.csv')
    net = net[net.model.isin(['home_net', 'hyb7_net'])].assign(per=lambda d: d.rels/d.n)
    fig, axs = plt.subplots(2, 4, figsize=(14, 6), sharex=True)
    for i, s in enumerate((0, 1)):
        for j, t in enumerate(TYPES):
            ax = axs[i, j]
            for k, m in enumerate(('home_net', 'hyb7_net')):
                d = net[(net.model == m) & (net.female == s) & (net.type == t)]
                agg = d.groupby(['year', 'age']).per.agg(['mean', 'std']).reset_index()
                for a, mk in zip(NETLAB, 'os^'):
                    r = agg[agg.age == a]
                    ax.errorbar(r.year + (k - 0.5)*1.2, r['mean'], yerr=r['std'], color=col[m], marker=mk, ls='-' if a == '15-24' else ':', ms=4, lw=0.8,
                                label=f'{"EMOD′" if m == "home_net" else "7"} {a}' if (i, j) == (0, 0) else None)
            ax.set_title(f'{["Men", "Women"][s]}, {dict(T="transitory", I="informal", M="marital", C="commercial")[t]}')
    axs[0, 0].set_ylabel('Relationships per person')
    axs[1, 0].set_ylabel('Relationships per person')
    axs[0, 0].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle('EmodNetwork: current relationships per person by type, sex and age, EMOD′ (grey) vs hybrid 7 (red), mean ± SD over 6 seeds')
    fig.tight_layout()
    fig.savefig(figs / 'network.png', dpi=120)

    epi = pd.read_csv(out / 'epi.csv')
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for m in ('home_hiv', 'hyb6', 'hyb7'):
        d = epi[epi.model == m]
        for seed, r in d.groupby('seed'):
            axs[0].plot(r.year, r.inf1549/r.n1549, color=col[m], lw=1, label=MODELS[m][0].split(':')[0] if seed == SEEDS[0] else None)
            axs[1].plot(r.year, r.high1549/r.n1549, color=col[m], lw=1)
    axs[0].set_title('HIV prevalence 15–49 (no ART in any)')
    axs[1].set_title('Share of 15–49 at Risk HIGH')
    axs[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figs / 'epidemic.png', dpi=120)
    return


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--compare', action='store_true')
    args = ap.parse_args()
    if args.plot:
        plot()
    elif args.compare:
        compare()
    else:
        run()
