"""
EMOD′ Risk and sex work: the `Risk` individual property and the campaign events that make people female sex workers
(FSW) or clients (docs/emod_inventory.md, section 2).

- `Risk` is drawn when a person is first seen (at birth, or at the start) from their node's initial distribution:
  LOW with probability p_low[node], else MEDIUM (the Risk overlay, by NodeID, as Daniel's inputs have it).
- At sexual debut, a woman will ever be an FSW, and a man a client, with a probability by node. After a uniform delay
  (up to 5 years for women, 10 for men) their Risk becomes HIGH; after a Weibull time (women) or a uniform one (men)
  it becomes MEDIUM, whatever it was before. The listeners that start the delay begin in 1961, so only debuts from 1961
  count (and nobody already past debut at the start).
- The network redraws a person's flags and maxima when their Risk changes.

It runs as an intervention, after the network, as EMOD's campaign runs after its pair formation.
"""
import numpy as np
import starsim as ss

LOW, MEDIUM, HIGH = 0, 1, 2


class EmodRisk(ss.Intervention):
    """
    Args:
        inp (objdict): EMOD's inputs for one parameter set, from `inputs.load()`
        net (str): the name of the network whose debut ages it reads
    """
    def __init__(self, inp, net='emodnet', name='emodrisk', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp, self.netname = inp, net
        self.define_states(
            ss.BoolArr('debuted', default=False),
            ss.FloatArr('t_uptake'),  # Year their Risk becomes HIGH
            ss.FloatArr('t_dropout'), # Year it becomes MEDIUM again
        )
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.risk = sim.people.add_state(ss.IntArr('risk', default=-1))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 3)
        return

    def start_step(self):
        """ Risk for anyone who has none yet """
        super().start_step()
        auids = self.sim.people.auids
        uids = auids[self.risk.raw[auids] < 0]
        low = self.rng.random(len(uids)) < self.inp.p_low[self.node.raw[uids]]
        self.risk[uids] = np.where(low, LOW, MEDIUM)
        return

    def step(self):
        ppl, sw = self.sim.people, self.inp.sexwork
        year = self.t.now('year')
        net = self.sim.networks[self.netname]
        auids = ppl.auids
        female = ppl.female.raw

        # Debut: whether they will ever sell or buy sex, and when they start
        new = auids[~self.debuted.raw[auids] & (ppl.age.raw[auids] >= net.debut.raw[auids])]
        self.debuted[new] = True
        if year >= sw.start_year and len(new):
            sex = female[new].astype(int)
            ever = new[self.rng.random(len(new)) < sw.p_ever[self.node.raw[new], sex]]
            delay = self.rng.random(len(ever))*sw.delay_max[female[ever].astype(int)]
            self.t_uptake[ever] = year + delay/365

        # Uptake and dropout
        due = lambda arr: auids[arr.raw[auids] <= year]
        start = due(self.t_uptake)
        if len(start):
            self.risk[start] = HIGH
            self.t_uptake[start] = np.nan
            f = female[start]
            days = np.where(f, sw.dropout[1]['Lambda']*(-np.log(self.rng.random(len(start))))**(1/sw.dropout[1]['Kappa']),
                            sw.dropout[0]['Min'] + self.rng.random(len(start))*(sw.dropout[0]['Max'] - sw.dropout[0]['Min']))
            self.t_dropout[start] = year + days/365
        stop = due(self.t_dropout)
        if len(stop):
            self.risk[stop] = MEDIUM
            self.t_dropout[stop] = np.nan
        return
