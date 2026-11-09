"""
Optima′ as agents: the declared model run with Starsim's agent solver, with agent versions of the two pieces that write
the ODE engine's compartments (population forcing and Optima's cascade rules). The pieces that read compartments
(Transmission, Births, Regain, VMMC) read optima_prime.state_counts(), which as agents counts the agents here.

The agent versions use only public mechanisms, and so also run as an ODE (make_sim(solver='ode')), which separates
what the agent versions change from what running as agents changes.
"""
import numpy as np, starsim as ss
import optima_prime as opp

_counts = dict(key=None, out=None) # This step's counts, by state: see counts()


def counts(sim, names, fresh=False):
    """ The number of people in each state, by population (states × populations), from the agents (each counts as sim.current_scale people); counted once a step unless `fresh` """
    key = (id(sim), sim.ti)
    if fresh or _counts['key'] != key:
        dis = sim.diseases.optima
        if isinstance(sim.solver, ss.solvers.abm):
            out = {s.name: sim.groups.count(dis[s.name])*sim.current_scale for s in dis.graph.people_states}
        else:
            eng = dis.engine
            out = {n: eng.y[c].sum(axis=0) for n, c in eng.state_comps.items()}
        _counts.update(key=key, out=out)
    return np.array([_counts['out'][n] for n in names])


def stale():
    """ Count again at the next call (after moving people) """
    _counts['key'] = None


class ForcePopsize:
    """
    Optima's population forcing (optima_prime.ForcePopsize) as declared flows, so that it can add people as agents:
    from the counts at the start of each step, the people each population gains or loses so that, as Optima does,
    those with no inflows gain or lose susceptibles to match their size, and then everyone not on ART is scaled to match it.
    Gains are arrivals (ss.immigration) into each state not on ART, in proportion to who is in it; losses, departures
    (ss.emigration) at the rate that removes the share. Each is taken over the step, so the population matches its
    target of the step's start at its end, as the ODE's rescaling after each step does.
    """
    def __init__(self, ref):
        sp = ref.simpars
        self.tvec, self.popsize = sp['tvec'], np.asarray(sp['popsize'], dtype=float)
        inflows = np.asarray(sp['agetransit']).sum(axis=0) + np.asarray(sp['birthtransit']).sum(axis=0)
        self.noinflows = np.nonzero(inflows == 0)[0]
        self.names = [s.replace('-', '_') for s in ref.settings['statelabels']]
        self.art = np.array([n.startswith(tuple(opp.TREATED)) for n in self.names])
        self.sus = np.array([n in ('susreg', 'progcirc') for n in self.names])
        self._ti, self._out = None, None

    def __deepcopy__(self, memo):
        return self # Shared with its transitions' functions

    def flows(self, sim):
        """ Arrivals (people per year) and departures (rate per person per year) by state × population """
        if self._ti != sim.ti:
            y0 = counts(sim, self.names)
            y = y0.copy()
            target = np.array([np.interp(sim.now.years, self.tvec, row) for row in self.popsize])
            j = self.noinflows
            new = target[j] - y[:, j].sum(axis=0)
            y[np.ix_(self.sus, j)] *= 1 + new/np.maximum(y[np.ix_(self.sus, j)].sum(axis=0), 1e-9)
            wanted = target - y[self.art].sum(axis=0)
            y[~self.art] *= wanted/np.maximum(y[~self.art].sum(axis=0), 1e-9)
            delta, dt = y - y0, sim.t.dt_year
            arrive = np.maximum(delta, 0)/dt
            leave = -np.log1p(np.maximum(np.minimum(delta, 0)/np.maximum(y0, 1e-9), -1 + 1e-9))/dt
            self._ti, self._out = sim.ti, (arrive, leave)
        return self._out

    def transitions(self, st, inp):
        """ For optima_prime.make_disease(extra=) """
        trs = {}
        for i, n in enumerate(self.names):
            if self.art[i]: continue
            trs[f'force_in_{n}'] = ss.immigration(st[n], n=ss.by(pop=lambda sim, i=i: ss.peryear(self.flows(sim)[0][i])))
            trs[f'force_out_{n}'] = ss.emigration(st[n], rate=ss.by(pop=lambda sim, i=i: ss.peryear(self.flows(sim)[1][i])))
        return trs


class CascadeTargets(opp.CascadeTargets):
    """
    Optima's cascade rules (optima_prime.CascadeTargets) with ss.move(): the same numbers, worked out from the counts
    after the step, moved as shares of each state (in each population in proportion, as Optima moves them)
    """
    def step(self):
        sim, dis = self.sim, self.sim.diseases.optima
        rows = {c: [f'{c}_{k}' for k in opp.CD4] for c in opp.UNTREATED + opp.TREATED}
        names = [n for c in rows.values() for n in c]
        y = dict(zip(names, counts(sim, names, fresh=True).sum(axis=1))) # Totals by state, kept up to date with each move
        tot = lambda cs: sum(y[n] for c in cs for n in rows[c])
        def move(src, dst, amount=None, p=None):
            """ Move `amount` people from one state to another, as a share of the source (in expectation, as agents); or each person with probability `p` """
            if p is None:
                if amount <= 0 or y[src] <= 0: return
                p = min(amount/y[src], 1)
            moved = ss.move(dis[src], dis[dst], p=p)
            y[src] -= moved; y[dst] += moved
        year = sim.now.years
        structure = dict(propdx=(['undx'], 'dx', ['dx', 'care', 'lost', 'usvl', 'svl'], opp.UNTREATED + opp.TREATED),
                         propcare=(['dx', 'lost'], 'care', ['care', 'usvl', 'svl'], ['dx', 'care', 'lost', 'usvl', 'svl']),
                         proptx=(['care'], None, opp.TREATED, ['care'] + opp.TREATED), propsupp=(['usvl'], 'svl', ['svl'], opp.TREATED))
        for name, (lower, upper, numer, denom) in structure.items():
            if name not in self.held and year >= self.fixyears[name] - 1e-9:
                self.held[name] = tot(numer)/(self.eps + tot(denom))
            actual, available = tot(numer), tot(denom)
            if name in self.held:
                wanted = self.held[name]*available
            elif name == 'proptx':
                wanted = np.interp(year, self.tvec, self.numtx)
            else:
                continue
            diff = wanted - actual
            if diff > self.eps: # Up the cascade
                up = np.array([[y[n] for n in rows[c]] for c in lower]) # Lower states × CD4
                total = up.sum()
                if total <= self.eps: continue
                diff = min(diff, total - self.eps)
                if name == 'proptx' and year < self.eligible: # Lowest CD4 first
                    movers = np.zeros_like(up)
                    for i in reversed(range(len(opp.CD4))):
                        if diff > self.eps and up[:, i].sum() > self.eps:
                            n = min(diff, up[:, i].sum() - self.eps)
                            movers[:, i] = n*up[:, i]/up[:, i].sum()
                            diff -= movers[:, i].sum()
                else:
                    movers = diff*up/total
                for a, c in enumerate(lower):
                    for i, k in enumerate(opp.CD4):
                        src = f'{c}_{k}'
                        if name == 'proptx' and y[src] > 0: # Suppressed at the start or not: each person's chances, the second given not chosen by the first (as shares of those left, the second would be capped at 1 when everyone in the state starts, and too few would be unsuppressed)
                            q = min(movers[a, i]/y[src], 1)
                            move(src, f'svl_{k}', p=q*self.treatvs)
                            move(src, f'usvl_{k}', p=q*(1 - self.treatvs)/(1 - q*self.treatvs))
                        else:
                            move(src, f'{upper}_{k}', movers[a, i])
            elif diff < -self.eps: # Down the cascade
                srcs = opp.TREATED if name == 'proptx' else [upper]
                down = np.array([[y[n] for n in rows[c]] for c in srcs])
                total = down.sum()
                if total <= self.eps: continue
                movers = min(-diff, total - self.eps)*down/total
                dst = 'care' if name == 'proptx' else 'lost' if name == 'propcare' else lower[0]
                for a, c in enumerate(srcs):
                    for i, k in enumerate(opp.CD4):
                        move(f'{c}_{k}', f'{dst}_{k}', movers[a, i])
        stale()
        return


def make_sim(ref=None, components=('demography', 'natural_history', 'transmission', 'cascade', 'treatment', 'mtct'),
             dt=0.02, end=2030, n_agents=100_000, solver='abm', rand_seed=1, interventions=(), cascade=None, **kwargs):
    """
    Optima′ (Optima's own cascade) with the agent versions of its forcing and cascade rules: as agents, or with
    solver='ode' as an ODE; `interventions` run after them. With `cascade` (a common/cascade Params of pathways by rate
    and level), Optima′ takes the table's coverages and leaves its cascade's moves and rules to the agent cascade, which
    `interventions` should then hold (pathways.TableCascade with host.OptimaHost)
    """
    ref = ref or opp.load()
    fp = ForcePopsize(ref) if 'demography' in components and ref.simpars['forcepopsize'] else None
    cx = opp.table_cascade(ref, cascade) if cascade is not None else None
    dis, st = opp.make_disease(ref, components, extra=fp.transitions if fp else None, cascade=cx, table_moves=False)
    interventions = ([CascadeTargets(ref, name='cascade_targets')] if 'treatment' in components and cx is None else []) + list(interventions)
    total = ref.raw['people'][:, :, 0].sum()
    agents = dict(n_agents=n_agents, rand_seed=rand_seed) if solver == 'abm' else {}
    return ss.Sim(dis, interventions=interventions, solver=solver, start=1990, stop=end, dt=dt, total_pop=total, verbose=0, **agents, **kwargs)
