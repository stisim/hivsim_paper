"""
EMOD′ sims: the demography, network and risk modules for one of Daniel's parameter sets, with a yearly analyzer
that counts what EMOD's ReportHIVByAgeAndGender does.

    sim = make_sim(pid=0, seed=1); sim.run()
"""
import numpy as np
import pandas as pd
import starsim as ss

import inputs
from demography import EmodDemography
from network import EmodNetwork, KEYS, available_all
from risk import EmodRisk
from hiv import EmodHIV

AGES = np.arange(0, 70, 5) # Lower edges of the report's age bins, as in the EMOD run


class Snapshots(ss.Analyzer):
    """ Once a year, at the times of EMOD's ReportHIVByAgeAndGender: people and their relationships by node, sex, Risk and age """
    def __init__(self, net='emodnet', name='snapshots', **kwargs):
        super().__init__(name=name, **kwargs)
        self.netname = net
        self.rows = []
        return

    def step(self):
        if (self.ti - 5) % 12: # EMOD reports at 1960.92, 1961.92, ...
            return
        ppl, net = self.sim.people, self.sim.networks.get(self.netname)
        uids = ppl.auids
        uids = uids[~np.isin(uids, self.sim.demographics.emoddemog.born) & ppl.alive.raw[uids]] # EMOD reports before the step's births
        age = AGES[np.searchsorted(AGES*365, ppl.emod_age.raw[uids], side='right') - 1] # EMOD's age in days at its report
        df = pd.DataFrame(dict(node=ppl.node.raw[uids] + 1, female=ppl.female.raw[uids].astype(int), age=age, pop=1))
        if net is None:
            out = df.groupby(['node', 'female', 'age']).sum().reset_index()
            out.insert(0, 'year', round(self.t.now('year') + 1/12, 2))
            self.rows.append(out)
            return
        n = net.stack('n', uids)
        q = net.stack('q', uids)
        mx = net.stack('max', uids)
        debuted = ppl.emod_age.raw[uids]/365 >= net.debut.raw[uids]
        avail = available_all(n, q, mx, net.flags.raw[uids].astype(np.int64)) & debuted[:, None]
        df = df.assign(risk=net.risk.raw[uids], debuted=debuted, concurrent=n.sum(1) > 1, partners=n.sum(1))
        for t, k in enumerate(KEYS):
            df[f'cur_{k}'] = n[:, t] > 0
            df[f'num_{k}'] = n[:, t]
            df[f'que_{k}'] = q[:, t]
            df[f'avail_{k}'] = avail[:, t]
        hiv = self.sim.diseases.get('hiv')
        if hiv is not None:
            st = hiv.stage.raw[uids]
            inf = hiv.infected.raw[uids]
            df = df.assign(infected=inf, acute=inf & (st == 0), latent=inf & (st == 1), aids=inf & (st == 2), circumcised=ppl.circumcised.raw[uids], sti=ppl.sti.raw[uids])
        out = df.groupby(['node', 'female', 'risk', 'age']).sum().reset_index()
        out.insert(0, 'year', round(self.t.now('year') + 1/12, 2))
        self.rows.append(out)
        return

    def finalize(self):
        super().finalize()
        self.df = pd.concat(self.rows, ignore_index=True)
        return


def make_sim(pid=0, seed=1, stop=2025, record=True, network=True, hiv=False, **kwargs):
    """ EMOD′ for parameter set `pid`: its demography; with `network`, its sexual network and Risk; with `hiv`, HIV too """
    inp = inputs.load(pid)
    mods = dict(networks=EmodNetwork(inp, record=record), interventions=EmodRisk(inp)) if network or hiv else {}
    if hiv:
        mods['diseases'] = EmodHIV(inp, record=record)
    return ss.Sim(n_agents=int(inp.init_pop.sum()), start=inp.start, stop=stop, dt=1/12, rand_seed=seed, verbose=0,
                  demographics=EmodDemography(inp), analyzers=Snapshots(), **mods, **kwargs)
