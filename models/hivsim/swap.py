"""
HIVsim Malawi (malawi/, hivsim_malawi at 803b2f7) with its own care, or with the condensed cascade (common/cascade)
in its place, with EMOD's Malawi table (set 0).

    sim = make_sim('condensed', seed=1); sim.run()
"""
import sys
from pathlib import Path
import pandas as pd
import stisim as sti

here = Path(__file__).resolve().parent
repo = here.parents[1]
sys.path[:0] = [str(here / 'malawi'), str(here), str(repo)]
import malawi # noqa: E402
from common.cascade.params import Params # noqa: E402

TABLES = repo / 'common' / 'cascade' / 'tables'

# hivsim_malawi experiment 03's best run (outputs/coverage_summary.json): calibration is still under way
BEST = dict(hiv_beta_m2f=0.007978243787796735, hiv_rel_init_prev=1.0839243919089498, nw_prop_f0=0.8315692935579891,
            nw_prop_m0=0.612858870804095, nw_f1_conc=0.052841444966013285, nw_m1_conc=0.1988339487797292,
            nw_fsw_shares=0.05474597622860841, nw_client_shares=0.23048969675026132, test_scale=1.1786373930647522,
            condom_scale=1.286413141369469, test_rel_m=1.9641771041623357, test_start=2009.089939951)


def emod_table(set_id=0):
    """ EMOD's Malawi table for HIVsim: flows summed over EMOD's three nodes, and scaled by HIVsim's pop_scale """
    params = pd.read_csv(TABLES / f'malawi_emod_set{set_id:03d}_params.csv')
    params = params[params.pathway != 'flows']
    flows = pd.read_csv(TABLES / f'malawi_emod_set{set_id:03d}_flows.csv', keep_default_na=False)
    keys = [c for c in flows.columns if c not in ('node', 'n')]
    flows = flows.groupby(keys, as_index=False).n.sum().assign(node='all')
    return Params(params, flows)


def traditional_vmmc():
    """ Circumcision before the programmes (from 2010), held at its level then: HIVsim's VMMC data before 2010 """
    v = malawi.read('vmmc.csv')
    pre = v[v.Year <= 2010].groupby('AgeBin', as_index=False).p_vmmc.first()
    return sti.VMMC(coverage=pd.concat([pre.assign(Year=y) for y in (1985, 2100)])[['Year', 'AgeBin', 'p_vmmc']])


def make_sim(variant='own', seed=1, pars=None, params=None, n_agents=10e3, stop=2022, **kwargs):
    """ HIVsim Malawi with its own testing, ART and VMMC (`own`), or with the condensed cascade (`condensed`, with `params`, by default emod_table()); from the slot factories (components.py) """
    import components
    pars = BEST if pars is None else pars
    return components.make_sim(pars, seed=seed, n_agents=n_agents, stop=stop, verbose=0, cascade=variant, params=params, sim_kwargs=kwargs)
