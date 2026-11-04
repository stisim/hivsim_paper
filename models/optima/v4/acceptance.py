"""P2 acceptance: Optima′ (the whole model) against the calibrated Optima run, each output and year, 2000-2020, overall and by population"""
import numpy as np, sciris as sc
import optima_prime as opp, compare as cmp

START, END = 2000, 2020
OUTPUTS = dict(plhiv=slice(2, 38), diag=slice(8, 38), treat=slice(26, 38), supp=slice(32, 38))

if __name__ == '__main__':
    refs = sc.load(cmp.REFS)
    o = refs['E']['original']
    years = o['years']
    sel = (years >= START) & (years <= END)
    T = sc.timer()
    sim = opp.make_sim(components=cmp.COMPONENTS['E'], dt=refs['E']['dt'], end=int(years[-1])).run()
    print(f'Optima′ ran in {T.toc(output=True):.1f} s')
    p = cmp.prime_people(sim, years)
    pk = list(opp.load().simpars['popkeys'])
    rows = {}
    for name, sl in OUTPUTS.items():
        a, b = o['people'][sl].sum(axis=0), p[sl].sum(axis=0) # Pops × years
        rows[name] = (100*np.abs(b.sum(0)/a.sum(0) - 1)[sel], 100*np.abs(b/np.where(a > 0, a, np.nan) - 1)[:, sel])
    a, b = o['people'].sum(axis=0), p.sum(axis=0)
    rows['popsize'] = (100*np.abs(b.sum(0)/a.sum(0) - 1)[sel], 100*np.abs(b/a - 1)[:, sel])
    rows['inci'] = (100*np.abs(cmp.prime_inci(sim, years)/o['inci'].sum(0) - 1)[sel], None)
    rows['deaths'] = (100*np.abs(cmp.prime_deaths(sim, years)/o['death'].sum(axis=(0, 1)) - 1)[sel], None)
    print(f'Max % difference from the calibrated Optima run, each year {START}-{END}')
    print(f'{"output":10s}{"overall":>9s}{"worst year":>12s}{"by pop":>9s}{"worst pop":>12s}')
    for name, (tot, bypop) in rows.items():
        yr = years[sel][np.argmax(tot)]
        if bypop is not None:
            j = np.nanargmax(np.nanmax(bypop, axis=1))
            print(f'{name:10s}{tot.max():9.2f}{yr:12d}{np.nanmax(bypop):9.2f}{pk[j]:>12s}')
        else:
            print(f'{name:10s}{tot.max():9.2f}{yr:12d}')
    print('\nPopulations over 2% in any year, by output (worst %, and its share of the output in 2020):')
    for name, (tot, bypop) in rows.items():
        if bypop is None: continue
        sl = OUTPUTS.get(name)
        share = (o['people'][sl] if sl else o['people']).sum(axis=0)[:, -11]
        share = share/share.sum()
        over = [(pk[j], round(float(np.nanmax(bypop[j])), 1), f'{100*share[j]:.1f}%') for j in np.argsort(-np.nanmax(bypop, axis=1)) if np.nanmax(bypop[j]) > 2]
        print(f'  {name:8s} {len(over):2d}: {over}')
    print('\nNew infections, % difference by year:', ' '.join(f'{y}:{d:.2f}' for y, d in zip(years[sel], rows["inci"][0])))
