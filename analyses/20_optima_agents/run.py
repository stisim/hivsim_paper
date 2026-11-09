"""
Analysis 20: Optima′ as agents against its ODE (see README.md)

    python run.py ode                # The original ODE, and the agent formulation as an ODE
    python run.py abm SEED [N]       # One run as agents (N agents, default 100,000)
    python run.py compare            # Tables
"""
import sys
from pathlib import Path
import numpy as np, sciris as sc
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'models' / 'optima' / 'v4'))
import optima_prime as opp, compare as cmp, agents as ag # noqa: E402

OUT = HERE / 'outputs'
YEARS = np.arange(1990, 2022)
END = 2021
SEL = (YEARS >= 2000) & (YEARS <= 2020)
OUTPUTS = dict(popsize=slice(0, 38), plhiv=slice(2, 38), diag=slice(8, 38), treat=slice(26, 38), supp=slice(32, 38))


def extract(sim):
    """ Whole-year stocks by population, and flows per year: in the step at each whole year (Optima's convention) and over each year """
    ref = opp.load()
    people = cmp.prime_people(sim, YEARS) # States × pops × years
    out = sc.objdict(years=YEARS, people=people)
    for k, sl in OUTPUTS.items():
        out[k] = people[sl].sum(axis=0) # Pops × years
    r = sim.results.optima
    t = np.array([d.years for d in r['infection_susreg'].timevec])
    yr = np.floor(t + 1e-9).astype(int)
    at = np.array([np.argmin(abs(t - y)) for y in YEARS])
    dt = t[1] - t[0]
    def flow(keys, weight=None):
        tot = sum(np.asarray(r[k])*(1 if weight is None else weight(k)) for k in keys)
        return tot[at]/dt, np.array([tot[yr == y].sum() for y in YEARS])
    inf = [k for k in r.keys() if k.startswith('infection_') or k in ('birth_undx', 'birth_dx')]
    out.inci_step, out.inci_year = flow(inf)
    out.sexinci_step, out.sexinci_year = flow([k for k in inf if k.startswith('infection_')])
    rr = np.array([np.interp(t, ref.simpars['tvec'], row) for row in ref.simpars['rrcomorbiditydeathtx']]) # Pops × steps: Optima reports deaths on ART scaled
    deaths = 0
    for k in [k for k in r.keys() if k.startswith('hivdeath_')]:
        res = np.asarray(r[k].by('pop')) # Steps (of the group results) × pops
        tg = np.array([d.years for d in r[k].by('pop').timevec])
        w = rr[:, [np.argmin(abs(t - x)) for x in tg]].T if k.split('_')[1] in opp.TREATED else 1
        deaths = deaths + (res*w).sum(axis=1)
    out.deaths_step = deaths[[np.argmin(abs(tg - y)) for y in YEARS]]/(tg[1] - tg[0])
    return out


def group_results():
    names = [s.replace('-', '_') for s in opp.load().settings['statelabels']]
    return names + [f'hivdeath_{n}' for n in names if not n.startswith('sus') and n != 'progcirc']


def run_ode():
    OUT.mkdir(exist_ok=True)
    for label, make in [('ode', lambda: opp.make_sim(components=cmp.COMPONENTS['E'], end=END)), ('ode_agentpieces', lambda: ag.make_sim(solver='ode', end=END))]:
        T = sc.timer()
        sim = make().run()
        el = T.toc(output=True)
        out = extract(sim)
        out.time = el
        sc.save(OUT / f'{label}.obj', out)
        print(f'{label}: {el:.1f} s')


def run_abm(seed, n=100_000):
    OUT.mkdir(exist_ok=True)
    T = sc.timer()
    sim = ag.make_sim(n_agents=n, rand_seed=seed, end=END, group_results=group_results())
    sim.init()
    tinit = T.toc(output=True)
    sim.run()
    out = extract(sim)
    out.time, out.time_init, out.n_agents, out.final_agents = T.toc(output=True), tinit, n, len(sim.people.auids)
    sc.save(OUT / f'abm_{n}_{seed}.obj', out)
    print(f'abm seed {seed}, {n} agents: init {tinit:.0f} s, total {out.time:.0f} s, {out.final_agents} agents alive at the end')


def compare(n=100_000):
    pk = list(opp.load().simpars['popkeys'])
    ode, ode2 = sc.load(OUT / 'ode.obj'), sc.load(OUT / 'ode_agentpieces.obj')
    runs = [sc.load(f) for f in sorted(OUT.glob(f'abm_{n}_*.obj'))]
    lines = [f'Optima′ as agents ({len(runs)} seeds, {n:,} agents) against its ODE, 2000-2020, each whole year',
             f'Run time: ODE {ode.time:.0f} s; agents mean {np.mean([r.time for r in runs]):.0f} s (init {np.mean([r.time_init for r in runs]):.0f} s)', '']
    totals = lambda o: {**{k: o[k].sum(axis=0) for k in OUTPUTS}, 'inci_step': o.inci_step, 'inci_year': o.inci_year, 'sexinci_year': o.sexinci_year, 'deaths_step': o.deaths_step}
    a, b = totals(ode), totals(ode2)
    lines.append('Agent pieces as an ODE vs the original ODE: max |%| difference')
    lines.append('  ' + '  '.join(f'{k} {100*np.max(np.abs(b[k][SEL]/a[k][SEL] - 1)):.2f}' for k in a))
    lines.append('')
    m = {k: np.array([totals(r)[k] for r in runs]) for k in a}
    for refname, ref in [('original ODE', a), ('agent pieces as ODE', b)]:
        lines.append(f'Agents vs {refname}: mean % difference (ABM mean/ODE - 1), max |%|, and z = (mean - ODE)/SE: mean z, share of years |z| > 2')
        lines.append(f'{"output":13s}{"mean %":>9s}{"max |%|":>9s}{"mean z":>9s}{"|z|>2":>8s}{"SE %":>8s}')
        for k in a:
            mean, se = m[k].mean(0)[SEL], m[k].std(0, ddof=1)[SEL]/np.sqrt(len(runs))
            d = 100*(mean/ref[k][SEL] - 1)
            z = (mean - ref[k][SEL])/se
            lines.append(f'{k:13s}{d.mean():9.2f}{np.abs(d).max():9.2f}{z.mean():9.2f}{np.mean(np.abs(z) > 2):8.2f}{100*np.mean(se/ref[k][SEL]):8.2f}')
        lines.append('')
    yrs = [2000, 2005, 2010, 2015, 2020]
    lines.append('By year (ODE | agents mean ± SE):')
    for k in a:
        lines.append(f'  {k:12s} ' + '  '.join(f'{y}: {a[k][YEARS == y][0]:,.0f} | {m[k].mean(0)[YEARS == y][0]:,.0f}±{m[k].std(0, ddof=1)[YEARS == y][0]/np.sqrt(len(runs)):,.0f}' for y in yrs))
    lines.append('')
    lines.append('By population, 2000-2020: |z| > 3 in more than a quarter of the years (z of ABM mean vs original ODE), with the mean % gap and the population\'s share of the output in 2020')
    for k in ['popsize', 'plhiv', 'diag', 'treat']:
        rs = np.array([r[k] for r in runs])[:, :, SEL] # Seeds × pops × years
        o = ode[k][:, SEL]
        mean, se = rs.mean(0), rs.std(0, ddof=1)/np.sqrt(len(runs))
        z = np.where(se > 0, (mean - o)/np.where(se > 0, se, 1), 0)
        share = o[:, -1]/o[:, -1].sum()
        bad = [(pk[j], f'{100*np.nanmean(mean[j]/np.where(o[j] > 0, o[j], np.nan) - 1):+.1f}%', f'{100*share[j]:.1f}%') for j in range(len(pk)) if np.mean(np.abs(z[j]) > 3) > 0.25]
        lines.append(f'  {k:8s} {len(bad):2d}/{len(pk)}: {bad}')
    text = '\n'.join(lines)
    print(text)
    (OUT / 'tables.txt').write_text(text + '\n')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ode': run_ode()
    elif cmd == 'abm': run_abm(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 100_000)
    elif cmd == 'compare': compare(int(sys.argv[2]) if len(sys.argv) > 2 else 100_000)
