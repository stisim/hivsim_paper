"""
HIVsim for Malawi: the uncalibrated skeleton for P3.

Data files go in hivsim/data/ (see the README there). Parameter values
marked "to calibrate" are placeholders until the P3 target set is agreed.

    python malawi.py          # short check run
"""
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

datadir = Path(__file__).parent / 'data'
location = 'malawi'
start = 1985 # Optima starts in 1990; the comparison is 2000-2020

required = [f'{location}_age_{start}.csv', f'{location}_asfr.csv', f'{location}_deaths.csv', 'init_prev_hiv.csv', 'n_art.csv']
optional = ['condom_use.csv', 'vmmc.csv', 'vls.csv']


def check_data():
    """ Stop with a list of the required data files that are missing """
    missing = [f for f in required if not (datadir/f).exists()]
    if missing:
        raise FileNotFoundError(f'Missing from {datadir}: {", ".join(missing)} (see the README there)')


def read(filename, **kwargs):
    """ Read an optional data file, or return None if it isn't there """
    path = datadir/filename
    return pd.read_csv(path, **kwargs) if path.exists() else None


def make_testing(test_scale=1, test_rel_m=1, test_start=2004):
    """
    HIV testing: general population, people with low CD4, and ANC. General and low-CD4 testing scale up linearly
    from test_start to 2020, to levels multiplied by test_scale; men's general testing is test_rel_m times women's
    """
    def scaleup(sim):
        return np.clip((sim.t.now('year') - test_start) / (2020 - test_start), 0, 1)

    def general_prob(sim, uids):
        rel = np.where(sim.people.female[uids], 1.0, test_rel_m)
        return np.minimum(0.1 * test_scale * scaleup(sim) * rel, 1)

    def low_cd4_prob(sim, uids):
        return np.full(len(uids), min(0.8 * test_scale * scaleup(sim), 1))

    undiagnosed = lambda sim: ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art
    general = sti.HIVTest(name='general_testing', test_prob_data=general_prob, start=1990, eligibility=undiagnosed)
    low_cd4 = sti.HIVTest(name='low_cd4_testing', test_prob_data=low_cd4_prob, start=1990,
                          eligibility=lambda sim: undiagnosed(sim) & (sim.diseases.hiv.cd4 < 200))
    anc = sti.ANCTest(visit_prob=0.9, start=2004) # Malawi's PMTCT programme; to calibrate
    return [general, low_cd4, anc]


def make_sim(seed=1, n_agents=10e3, stop=2021, verbose=1/12, hiv_pars=None, nw_pars=None, test_scale=1, test_rel_m=1, test_start=2004, condom_scale=1, **kwargs):
    """ Make the Malawi sim; hiv_pars and nw_pars are parameter values for sti.HIV and sti.StructuredSexual; condom_scale multiplies condom use (at most 1); see make_testing for the others """
    check_data()

    condoms = read('condom_use.csv')
    if condoms is not None:
        condoms = condoms.set_index('partnership').mul(condom_scale).clip(upper=1).reset_index()
    nw = sti.StructuredSexual(condom_data=condoms, **(nw_pars or {}))
    hiv = sti.HIV(**{
        'beta_m2f': 0.01, # To calibrate
        'init_prev_data': read('init_prev_hiv.csv'),
        'age_bins': [0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100], # The surveys' age groups
    } | (hiv_pars or {}))

    n_art = read('n_art.csv', index_col='year')
    vmmc = read('vmmc.csv') # Share of men circumcised, by age
    vls = read('vls.csv', index_col='year')
    intvs = make_testing(test_scale, test_rel_m, test_start) + [sti.ART(coverage=n_art, vls_coverage=vls)]
    if vmmc is not None:
        intvs.append(sti.VMMC(coverage=vmmc))

    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        datafolder=datadir,
        demographics=location,
        diseases=hiv,
        networks=[nw, ss.PrenatalNet()],
        interventions=intvs,
        rand_seed=seed,
        verbose=verbose,
        **kwargs,
    )
    return sim


if __name__ == '__main__':
    sim = make_sim(n_agents=2e3, verbose=0)
    with sc.timer():
        sim.run()
    df = sim.to_df(resample='year', use_years=True, sep='.')
    print(df[['time', 'hiv.prevalence_15_49', 'hiv.on_art']].iloc[::5].to_string())
