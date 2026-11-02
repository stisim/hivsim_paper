"""
HIVsim for Malawi: the uncalibrated skeleton for P3.

Data files go in models/hivsim/data/ (see the README there). Parameter values
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


def make_testing(test_scale=1):
    """ HIV testing: general population, people with low CD4, and ANC; test_scale multiplies the first two """
    years = np.arange(1990, 2031)
    scaleup = np.clip((years - 2004) / (2020 - 2004), 0, 1) # Testing scales up from 2004 to 2020; to calibrate
    undiagnosed = lambda sim: ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art
    general = sti.HIVTest(name='general_testing', years=years, test_prob_data=np.minimum(0.1*test_scale*scaleup, 1), eligibility=undiagnosed)
    low_cd4 = sti.HIVTest(name='low_cd4_testing', years=years, test_prob_data=np.minimum(0.8*test_scale*scaleup, 1),
                          eligibility=lambda sim: undiagnosed(sim) & (sim.diseases.hiv.cd4 < 200))
    anc = sti.ANCTest(visit_prob=0.9, start=2004) # Malawi's PMTCT programme; to calibrate
    return [general, low_cd4, anc]


def make_sim(seed=1, n_agents=10e3, stop=2021, verbose=1/12, hiv_pars=None, nw_pars=None, test_scale=1, condom_scale=1, **kwargs):
    """ Make the Malawi sim; hiv_pars and nw_pars are parameter values for sti.HIV and sti.StructuredSexual; condom_scale multiplies condom use (at most 1) """
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
    intvs = make_testing(test_scale) + [sti.ART(coverage=n_art, vls_coverage=vls)]
    if vmmc is not None:
        intvs.append(sti.VMMC(coverage=vmmc))

    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        datafolder=datadir,
        demographics=location,
        diseases=hiv,
        networks=[nw, ss.MaternalNet()],
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
