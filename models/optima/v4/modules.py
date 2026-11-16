"""
Optima′ as separate agent modules (docs/swap_design.md §2.3a, step 2): its demography and populations (D, with R),
transmission (T), natural history and ART biology (H) and seeding (S), which talk only through the H interface and
the host; the cascade (C) is `common.cascade.pathways.TableCascade` with Optima's table on `ModulesHost`
(components.py builds them). Each piece does what Optima′ declared as agents does (optima_prime.make_disease,
agents.py), with the same rates: Optima's probability per step of 0.2 years divided by 0.2, with Optima's products.

- `OptimaDR` (D and R, one module, §1.4): each agent's Optima population (`people.optima_pop`), with sex and an age in
  its band (fixed per agent, as `host.OptimaHost`'s); background death; population forcing (agents.ForcePopsize: from
  the counts at the start of the step, departures, and arrivals that are copies of residents of the population, of
  the same kind: without HIV, with HIV not on ART); Optima's births with mother-to-child transmission and PMTCT
  (optima_prime.Births: Optima's births have no mothers, so these read the mothers' places in care by population
  through the host); then ageing and risk transfers, in the declared order.
- `OptimaT`: Optima's force of infection by population pair (optima_prime.Transmission), from the counts at the start
  of the step: each source is the infected agents of a population in one place in care and CD4 band, with the sum of
  their `infectiousness()` (H) times Optima's `effdx` factor (T's, once for the diagnosed and once more on ART, as
  Optima's alltrans). Circumcision is `people.circumcised`. Sources of each step's infections are drawn after the
  step from the parts of the force, as host.OptimaTracker drew them.
- `OptimaH`: each infected agent's CD4 band (acute, gt500, gt350, gt200, gt50, lt50) and place (undx, dx, care, lost,
  svl, usvl: pathways.PLACES; ART's arm is svl or usvl, and the others are the host's record, kept here as EmodHIV
  keeps `diagnosed`). Each step, HIV death by band and arm, progression off ART, and the CD4 moves on ART, as one
  competing draw, among those who survive Optima's step (`hiv_survival()`, times D's `bg_survival()`).
- `OptimaS`: Optima's state in 1990, each population's shares of the 38 states, exactly (as the declared engine).

What crosses between modules: T → H `infect(uids)`; H → T `infectiousness(uids)`; T reads the place in care and CD4
band through the host (for effdx and its sources); C → H `start_art`/`stop_art`/`switch_arm` and the places, through
the host; H → D its deaths, and `hiv_survival()` (Optima's death before other moves: D's background death is among
those who survive HIV); D → H, T and C `bg_survival()`; D → H `infect()` for infected births. Within a step: D, then
T's infections, then H's moves (each among those not dying, and present at the start of the step), then C; the
declared engine instead draws D's, T's and H's moves for each agent in one competing draw, which differs only by
order dt² per step.
"""
import numpy as np
import starsim as ss
import optima_prime as opp
import poplabel as pl
from host import AGE_LO, AGE_W, FEMALE, CHILD, CD4_RANGE, fixed_fraction
from common.cascade.interface import CascadeHost
from common.cascade.pathways import PLACES

UNDX, DX, CARE, LOST, SVL, USVL = range(1, 7)
NB, NP = len(opp.CD4), len(pl.POPS)
NG = (len(PLACES) - 1)*NB # Sources: place × CD4 band
_REF = {}


def ref():
    """ Optima's inputs, loaded once """
    if 'ref' not in _REF:
        _REF['ref'] = opp.load()
    return _REF['ref']


def at(arr, year, tvec):
    """ Linear interpolation of a series over Optima's time points (last axis), as optima_prime.Transmission.at() """
    i = int(np.clip(np.searchsorted(tvec, year, side='right') - 1, 0, len(tvec) - 2))
    w = np.clip((year - tvec[i])/(tvec[i+1] - tvec[i]), 0, 1)
    return arr[..., i]*(1 - w) + arr[..., i+1]*w


def leaving(ppl, uids, ti):
    """ Who of these dies or leaves in this step (already requested) """
    return (ppl.ti_dead.raw[uids] == ti) | (ppl.ti_removed.raw[uids] == ti)


def modules(sim):
    """ The modules the others read: D (Optima's, or None), T, H """
    return (sim.demographics.get('optima_dr') if hasattr(sim.demographics, 'get') else None,
            sim.diseases.get('optima_t'), sim.diseases.get('hiv'))


def pop_of(sim, uids):
    return sim.people.optima_pop.raw[uids]


def bg_survival(sim, uids):
    """ D's chance of surviving background death in Optima's step, for each (1 if D isn't Optima's) """
    d = modules(sim)[0]
    return d.bg_survival(uids) if d is not None and hasattr(d, 'bg_survival') else np.ones(len(uids))


def round_preserving(y):
    return ss.solvers.Compartmental.round_preserving(y)


class OptimaDR(ss.Demographics):
    """
    Optima's demography and populations (D and R): see the module docstring

    Args:
        host (class): the CascadeHost of H's model (who has HIV, who is on ART, and the mothers' places in care)
    """
    def __init__(self, host, name='optima_dr', seed_offset=2401, **kwargs):
        super().__init__(name=name, **kwargs)
        self.host_class, self.seed_offset = host, seed_offset
        self.debuts = self.born = ss.uids()
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.optima_pop = sim.people.add_state(ss.IntArr('optima_pop', default=-1))
        self.rng = np.random.default_rng(sim.pars.rand_seed + self.seed_offset)
        r = ref()
        sp = r.simpars
        pl.check(r)
        self.tvec = sp['tvec']
        self.bg = np.asarray(sp['death'], float)*opp.REF_DT # Probability per 0.2-year step, pops × time
        assert np.all(np.asarray(sp['propemigrate']) == 0), 'Emigration (propemigrate) is not handled'
        self.popsize = np.asarray(sp['popsize'], float)
        inflows = np.asarray(sp['agetransit']).sum(axis=0) + np.asarray(sp['birthtransit']).sum(axis=0)
        self.noinflows = np.nonzero(inflows == 0)[0]
        self.forcing = bool(sp['forcepopsize'])
        b = opp.Births(r) # Optima's births and MTCT arithmetic (rates, effmtct, pmtcteff, numpmtct)
        self.rates, self.effmtct, self.pmtcteff, self.numpmtct = b.rates, b.effmtct, b.pmtcteff, b.numpmtct
        people0 = r.raw['people'][:, :, 0]
        self.fracs = people0.sum(axis=0)/people0.sum()
        self.transfers = self.make_transfers(sp)
        return

    def make_transfers(self, sp):
        """ Ageing (with Optima's old band skip) and risk transitions, as optima_prime.make_disease declares them, in its order: (src, dst, rate per year by time) """
        n = len(self.tvec)
        out = []
        series = lambda x: np.broadcast_to(np.asarray(x, float), (n,)).copy()
        at_, a = np.asarray(sp['agetransit']), sp['agerate']
        onward = {i: j for i, j in zip(*np.nonzero(at_))}
        for i, j in onward.items():
            rate = at_[i, j]*a[i]
            if opp.OPTIONS.age_skip and j in onward:
                k = onward[j]
                skip = rate*at_[j, k]*a[j]*opp.REF_DT
                out.append((i, k, series(opp.OPTIONS.skip_scale*skip)))
                rate = rate - opp.OPTIONS.skip_scale*skip
            out.append((i, j, series(opp.OPTIONS.age_scale*rate)))
        rt = np.asarray(sp['risktransit'], float)
        size = np.asarray(sp['popsize'], float)
        for i, j in zip(*np.nonzero(rt)) if opp.OPTIONS.risk else []:
            r = 1 - np.exp(-opp.REF_DT/rt[i, j])
            out.append((i, j, series(r/opp.REF_DT)))
            out.append((j, i, series(r/opp.REF_DT*size[i]/size[j])))
        return out

    def init_post(self):
        super().init_post()
        uids = self.sim.people.auids
        order = np.argsort(self.rng.random(len(uids)), kind='stable') # Exactly the split, in a random order, as the declared engine's groups
        counts = round_preserving(self.fracs*len(uids)).astype(int)
        self.optima_pop[uids[order]] = np.repeat(np.arange(NP), counts)
        self.set_people(uids)
        return

    def set_people(self, uids):
        """ Sex and age from the population (age fixed within its band for each agent, as host.OptimaHost's) """
        p = self.optima_pop.raw[uids]
        ppl = self.sim.people
        ppl.female[uids] = FEMALE[p]
        ppl.age[uids] = AGE_LO[p] + fixed_fraction(uids, 0)*AGE_W[p]
        return

    @property
    def host(self):
        if getattr(self, '_host', None) is None:
            self._host = self.host_class(self.sim)
        return self._host

    def year(self):
        return self.t.now('year')

    def bg_prob(self):
        """ Optima's probability of background death in its 0.2-year step, by population, now """
        return at(self.bg, self.year(), self.tvec)

    def bg_survival(self, uids):
        return 1 - self.bg_prob()[self.optima_pop.raw[uids]]

    def survival_by_pop(self):
        return 1 - self.bg_prob()

    def places(self, uids):
        """ Each one's class (0 without HIV, 1 with HIV not on ART, 2 on ART), and whether diagnosed """
        h = self.host
        inf = np.asarray(h.infected(uids), bool)
        p = h.care_place(uids)
        if p is not None:
            p = np.asarray(p)
            on, dx = np.isin(p, [SVL, USVL]), inf & (p != UNDX)
        else:
            st = h.cascade_stage(uids)
            on, dx = np.asarray(st['on_art'], bool), inf & np.asarray(st['diagnosed'], bool)
        return np.where(on, 2, np.where(inf, 1, 0)), dx

    def step(self):
        sim, ppl = self.sim, self.sim.people
        ti, dt, year, scale = self.ti, self.t.dt_year, self.year(), sim.current_scale
        au = ppl.auids
        pop = self.optima_pop.raw[au]
        cls, dx = self.places(au)
        y0 = np.bincount(cls*NP + pop, minlength=3*NP).reshape(3, NP)*scale # People, by class × population, at the start of the step

        arrive, leave = self.forcing_flows(y0, year, dt)

        # Background death (among those who survive HIV in Optima's step) and departures, as competing risks
        hb = self.bg_prob()[pop]/opp.REF_DT*self.hiv_survival(au)
        hl = leave[cls, pop]
        tot = hb + hl
        p = -np.expm1(-tot*dt)
        u = self.rng.random(len(au))
        go = u < p
        die = go & (u < p*np.divide(hb, tot, out=np.zeros_like(tot), where=tot > 0))
        ppl.request_death(au[die], cause=self)
        ppl.request_removal(au[go & ~die])

        # Births, with mother-to-child transmission (optima_prime.Births), from the start of the step
        amounts = self.births(pop, cls, dx, year, scale)
        kids = []
        for amount in amounts:
            n = self.rng.poisson(np.maximum(amount, 0)*dt/scale)
            new = self.grow(np.repeat(np.arange(NP), n))
            kids.append(new)
        dis = modules(sim)[2]
        infected = ss.uids(np.concatenate(kids[1:]))
        if len(infected):
            dis.infect(infected, source=None)
            if len(kids[2]):
                self.host.test(kids[2]) # Infected at birth to a diagnosed mother: diagnosed (Optima's birth_dx)
        self.born = ss.uids(np.concatenate(kids))

        # Arrivals (population forcing): copies of residents at the start of the step, of the population and class
        for c in (0, 1):
            n = self.rng.poisson(arrive[c]/scale)
            for j in np.flatnonzero(n):
                members = au[(cls == c) & (pop == j)]
                if not len(members): continue
                donors = members[self.rng.integers(len(members), size=n[j])]
                self.clone(self.grow(np.full(n[j], j)), donors)

        # Ageing and risk transfers, each in turn (as the declared engine takes them), for everyone now
        au = ppl.auids
        pop0 = self.optima_pop.raw[au].copy()
        p = pop0.copy()
        for src, dst, rate in self.transfers:
            idx = np.flatnonzero(p == src)
            if not len(idx): continue
            r = at(rate, year, self.tvec)
            moving = idx[self.rng.random(len(idx)) < -np.expm1(-max(r, 0)*dt)]
            p[moving] = dst
        self.optima_pop[au] = p
        self.debuts = au[CHILD[pop0] & ~CHILD[p] & (AGE_LO[p] >= 15)]
        self.set_people(au)
        return

    def forcing_flows(self, y0, year, dt):
        """ Population forcing (agents.ForcePopsize) by class (people: without HIV, with HIV not on ART, on ART × population, at the start of the step): populations with no inflows gain or lose people without HIV to match their size, then everyone not on ART is scaled to match it. Returns arrivals (people this step) and departures (rate per year) by class × population """
        arrive, leave = np.zeros((3, NP)), np.zeros((3, NP))
        if not self.forcing:
            return arrive, leave
        y = y0.copy()
        target = np.array([np.interp(year, self.tvec, row) for row in self.popsize])
        j = self.noinflows
        new = target[j] - y[:, j].sum(axis=0)
        y[0, j] *= 1 + new/np.maximum(y[0, j], 1e-9)
        wanted = target - y[2]
        y[:2] *= wanted/np.maximum(y[:2].sum(axis=0), 1e-9)
        delta = y - y0
        arrive = np.maximum(delta, 0)
        leave = -np.log1p(np.maximum(np.minimum(delta, 0)/np.maximum(y0, 1e-9), -1 + 1e-9))/dt
        return arrive, leave

    def births(self, pop, cls, dx, year, scale):
        """ Births per year into each population (optima_prime.Births, projects before 2.12.0): without HIV, with HIV to undiagnosed mothers, and with HIV to diagnosed mothers; from each agent's population, class and diagnosis at the start of the step """
        total = np.bincount(pop, minlength=NP)*scale
        undx = np.bincount(pop[(cls == 1) & ~dx], minlength=NP)*scale
        alldx = np.bincount(pop[dx], minlength=NP)*scale
        alltx = np.bincount(pop[cls == 2], minlength=NP)*scale
        rates = at(self.rates, year, self.tvec)
        eps = 1e-3
        pregdx = (alldx*rates.sum(axis=1)).sum()*opp.REF_DT
        proppmtct = min(at(self.numpmtct, year, self.tvec)*opp.REF_DT/(eps*opp.REF_DT + pregdx), 1)
        effmtct, pmtcteff = at(self.effmtct, year, self.tvec), at(self.pmtcteff, year, self.tvec)
        births = rates.T @ total
        mtctundx = rates.T @ undx*effmtct
        mtctdx = rates.T @ alldx*(proppmtct*pmtcteff + (1 - proppmtct)*effmtct) + rates.T @ alltx*pmtcteff # Mothers on ART count as both on treatment and eligible for PMTCT, as in Optima
        return births - mtctundx - mtctdx, mtctundx, mtctdx

    def hiv_survival(self, uids):
        h = modules(self.sim)[2]
        return h.hiv_survival(uids) if h is not None and hasattr(h, 'hiv_survival') else np.ones(len(uids))

    def grow(self, pops):
        """ New agents in these populations """
        if not len(pops):
            return ss.uids()
        new = self.sim.people.grow(len(pops))
        self.optima_pop[new] = pops
        self.set_people(new)
        return new

    SKIP = {'uid', 'slot', 'parent', 'alive', 'ti_dead', 'ti_removed', 'scale', 'cause_of_death', 'age', 'female', 'optima_pop'}

    def clone(self, new, donors):
        """ New agents as copies of others: every array of the people and modules but Starsim's own and the population's """
        for key, arr in self.sim.people.states.items():
            if key in self.SKIP or getattr(arr, 'shared', None) is not None:
                continue
            arr[new] = arr.raw[donors]
        return


class OptimaT(ss.Disease):
    """ Optima's force of infection by population pair (optima_prime.Transmission) on agents: see the module docstring """
    def __init__(self, name='optima_t', seed_offset=2203, **kwargs):
        super().__init__(name=name, **kwargs)
        self.seed_offset = seed_offset
        self.pairs = (ss.uids(), ss.uids())
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.circumcised = sim.people.add_state(ss.BoolArr('circumcised', default=False, result=False))
        self.rng = np.random.default_rng(sim.pars.rand_seed + self.seed_offset)
        r = ref()
        sp = r.simpars
        self.foi = opp.Transmission(r) # Its arrays: partnerships, acts, per-act probabilities, condoms, susceptibility (alleff), force, inhomo
        assert np.all(np.asarray(sp['propemigrate']) == 0), 'Emigration (propemigrate) is not handled'
        self.dxfactor = 1 - sp['effdx']
        self.n_start = 0
        self._new = []
        return

    @property
    def host(self):
        if getattr(self, '_host', None) is None:
            self._host = self.sim.diseases.hiv.host_class(self.sim)
        return self._host

    def sources(self, uids):
        """ Each infected agent's source group (place × CD4 band) and infectiousness (H's, times effdx: once diagnosed, and once more on ART) """
        h = self.host
        place = h.care_place(uids)
        if place is None:
            st = h.cascade_stage(uids)
            place = np.where(st['on_art'], np.where(st['suppressive'], SVL, USVL), np.where(st['diagnosed'], DX, UNDX))
        place = np.asarray(place)
        band = np.asarray(h.cd4_band(uids))
        rel = np.asarray(self.sim.diseases.hiv.infectiousness(uids), float)
        rel = rel*self.dxfactor**((place != UNDX).astype(int) + np.isin(place, [SVL, USVL]))
        return (place - 1)*NB + band, rel

    def start_step(self):
        """ The force of infection, from the agents at the start of the step (as Optima′'s counts) """
        super().start_step()
        sim, f = self.sim, self.foi
        ppl = sim.people
        self.n_start = ppl.uid.len_used
        year = self.t.now('year')
        au = ppl.auids
        pop = pop_of(sim, au)
        allpeople = np.maximum(np.bincount(pop, minlength=NP), 1e-300)
        inf = au[np.asarray(self.host.infected(au), bool)]
        g, rel = self.sources(inf)
        ip = pop_of(sim, inf)
        effprev = np.bincount(g*NP + ip, weights=rel, minlength=NG*NP).reshape(NG, NP)/allpeople # Sources × populations
        prev = np.bincount(ip, minlength=NP)/allpeople
        h = f.inhomopar
        inhomo = (h + f.eps)/(np.exp(h + f.eps) - 1)*np.exp(h*(1 - prev))
        acts = f.at(f.acts, year)*opp.REF_DT
        whole, frac = np.floor(acts), acts - np.floor(acts)
        alleff = f.at(f.alleff, year)[f.p1]
        q = (f.at(f.trans, year)*f.at(f.cond, year))[:, None, None]*alleff[:, :, None]*effprev[:, f.p2].T[:, None, :] # Partnerships × circ × sources
        noinf = (1 - np.minimum(frac[:, None, None]*q, 1))*(1 - np.minimum(q*(whole[:, None, None] != 0), 1))**whole[:, None, None]
        lognoinf = np.zeros((NP, NP) + q.shape[1:])
        np.add.at(lognoinf, (f.p1, f.p2), np.log(noinf))
        self.parts = 1 - np.exp(lognoinf) # Acquiring pop × partner pop × circ × source
        part = self.parts.sum(axis=1)
        bg = 1 - self.survival_by_pop()
        prob = part.sum(axis=2).T*f.force*inhomo*(1 - bg)
        prob = np.minimum(prob, 1 - f.eps - bg.max())
        self.rate = prob/opp.REF_DT # Circ × acquiring pop, per year
        return

    def survival_by_pop(self):
        d = modules(self.sim)[0]
        return d.survival_by_pop() if d is not None and hasattr(d, 'survival_by_pop') else np.ones(NP)

    def step_state(self):
        """ Infections, among those without HIV here at the start of the step and not dying or leaving in it """
        super().step_state()
        sim, ppl = self.sim, self.sim.people
        au = ppl.auids
        au = au[(au < self.n_start) & ~leaving(ppl, au, self.ti)]
        sus = au[~np.asarray(self.host.infected(au), bool)]
        c = (self.circumcised.raw[sus] & ~ppl.female.raw[sus]).astype(int)
        p = -np.expm1(-self.rate[c, pop_of(sim, sus)]*self.t.dt_year)
        hit = self.rng.random(len(sus)) < p
        new, c = sus[hit], c[hit]
        if len(new):
            sim.diseases.hiv.infect(new, source=None)
        self._new = (new, c)
        return

    def step(self):
        """ The sources of this step's infections (after the cascade, as host.OptimaTracker drew them): a source group from the parts of the force on the recipient's population and circumcision, then an agent in it now """
        super().step()
        sim, ppl = self.sim, self.sim.people
        new, circ = self._new if len(self._new) else (ss.uids(), np.zeros(0, int))
        if not len(new):
            self.pairs = (ss.uids(), ss.uids())
            return
        au = ppl.auids
        au = au[~leaving(ppl, au, self.ti)]
        cand = au[np.asarray(self.host.infected(au), bool)]
        g, _ = self.sources(cand)
        key = pop_of(sim, cand)*NG + g
        order = np.argsort(key, kind='stable')
        cand, key = cand[order], key[order]
        cells = np.arange(NP*NG)
        start = np.searchsorted(key, cells)
        count = np.searchsorted(key, cells, side='right') - start
        has = (count > 0).reshape(NP, NG)
        rpop = pop_of(sim, new)
        dst, src = [], []
        for cc in (0, 1):
            for i in np.unique(rpop[circ == cc]):
                who = new[(rpop == i) & (circ == cc)]
                w = self.parts[i, :, cc, :]*has
                if w.sum() <= 0:
                    w = has.astype(float)
                cell = self.rng.choice(w.size, size=len(who), p=(w/w.sum()).ravel())
                pick = start[cell] + (self.rng.random(len(who))*count[cell]).astype(int)
                dst.append(who); src.append(cand[pick])
        self.pairs = (ss.uids(np.concatenate(dst)), ss.uids(np.concatenate(src)))
        return


class OptimaH(ss.Disease):
    """ Optima's natural history and ART biology on agents: see the module docstring """
    def __init__(self, host=None, name='hiv', seed_offset=2501, **kwargs):
        super().__init__(name=name, **kwargs)
        self.host_class = host or ModulesHost
        self.seed_offset = seed_offset
        self.define_states(
            ss.BoolArr('infected', default=False, label='Infected'),
            ss.IntArr('place', default=0), # pathways.PLACES: 0 without HIV, then undx, dx, care, lost, svl, usvl
            ss.IntArr('band', default=-1), # CD4 band (optima_prime.CD4)
            ss.FloatArr('ti_infected'),
            ss.BoolArr('ever_treated', default=False, result=False),
            ss.BoolArr('symptomatic', default=False, result=False),
            ss.FloatArr('lowest_cd4', default=np.inf),
        )
        self.new_uids = self.dead_uids = self.newly_symptomatic = ss.uids()
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.rng = np.random.default_rng(sim.pars.rand_seed + self.seed_offset)
        sp = ref().simpars
        assert np.all(np.asarray(sp['hivdeath']) == 1), 'Relative HIV death rates other than 1 are not handled'
        R = opp.REF_DT
        d = np.array([sp[f'death{k}'] for k in opp.CD4])*R # Probability of HIV death in a step of 0.2 years, by band
        self.death = np.array([d, d*sp['deathsvl'], d*sp['deathusvl']]) # Off ART, on svl, on usvl × band
        tau = np.array([sp['progacute'], sp['proggt500'], sp['proggt350'], sp['proggt200'], sp['proggt50']])
        prog = 1 - np.exp(-R/tau)
        svlrecov = 1 - np.exp(-R/np.array([np.inf, np.inf, sp['svlrecovgt350'], sp['svlrecovgt200'], sp['svlrecovgt50'], sp['svlrecovlt50']]))
        u = lambda key: sp[key]*R
        up, down = np.zeros((3, NB)), np.zeros((3, NB)) # Probabilities per 0.2 years of moving one band down in CD4 (up the list) or up in CD4, as optima_prime's cd4moves
        up[0, :5] = prog # Off ART: progression (none from lt50)
        up[1, 0] = prog[0] # On svl: acute to gt500, recovery from gt350 down
        down[1, 2:] = svlrecov[2:]
        up[2, 0] = prog[0] # On usvl
        up[2, 1] = u('usvlproggt500')
        for i in range(2, 5):
            down[2, i], up[2, i] = u(f'usvlrecov{opp.CD4[i]}'), u(f'usvlprog{opp.CD4[i]}')
        down[2, 5] = u('usvlrecovlt50')
        self.up, self.down = up, down
        self.cd4trans = np.array([sp[f'cd4trans{k}'] for k in opp.CD4])*sp['transnorm']
        self.txfactor = np.array([1, 1 - sp['efftxsupp'], 1 - sp['efftxunsupp']])
        self._new = []
        self.n_start = 0
        return

    def arm(self, uids):
        """ 0 off ART, 1 suppressed (svl), 2 unsuppressed (usvl) """
        p = self.place.raw[uids]
        return np.where(p == SVL, 1, np.where(p == USVL, 2, 0))

    # The H interface
    def infect(self, uids, source=None):
        """ New infections (T, or D's births): into acute, undiagnosed """
        uids = ss.uids(uids)
        if len(uids):
            self.infected[uids] = True
            self.place[uids] = UNDX
            self.band[uids] = 0
            self.ti_infected[uids] = self.ti
            self._new.append(uids)
        return uids

    def seed(self, uids, band):
        """ For S: people with HIV at the start, in a CD4 band, undiagnosed (the host then places them in care) """
        self.infected[uids] = True
        self.place[uids] = UNDX
        self.band[uids] = band
        self.ti_infected[uids] = -1
        return

    def infectiousness(self, uids):
        """ Relative infectiousness per act: Optima's cd4trans by band (times transnorm), and its ART factors (effdx is T's) """
        b = self.band.raw[uids]
        out = np.where(b >= 0, self.cd4trans[np.maximum(b, 0)], 0.0)*self.txfactor[self.arm(uids)]
        return np.where(self.infected.raw[uids], out, 0.0)

    def hiv_survival(self, uids):
        """ The chance of surviving HIV in Optima's step of 0.2 years (1 without HIV): Optima applies death first, and other moves to those who survive """
        b = self.band.raw[uids]
        d = self.death[self.arm(uids), np.maximum(b, 0)]
        return np.where(self.infected.raw[uids], 1 - d, 1.0)

    def cd4(self, uids):
        out = np.full(len(uids), np.inf)
        f = fixed_fraction(uids, 1)
        b = self.band.raw[uids]
        for i, (lo, hi) in enumerate(CD4_RANGE.values()):
            m = self.infected.raw[uids] & (b == i)
            out[m] = lo + f[m]*(hi - lo)
        return out

    def on_art(self, uids):
        return np.isin(self.place.raw[uids], [SVL, USVL])

    def start_art(self, uids, arm):
        """ From any place off ART with HIV into svl ('suppressive') or usvl, in the same band """
        uids = ss.uids(uids)
        u = uids[self.infected.raw[uids] & ~self.on_art(uids)]
        self.place[u] = SVL if arm == 'suppressive' else USVL
        self.ever_treated[u] = True
        return u

    def stop_art(self, uids):
        uids = ss.uids(uids)
        u = uids[self.on_art(uids)]
        self.place[u] = LOST
        return u

    def switch_arm(self, uids, arm):
        uids = ss.uids(uids)
        src, dst = (USVL, SVL) if arm == 'suppressive' else (SVL, USVL)
        u = uids[self.place.raw[uids] == src]
        self.place[u] = dst
        return u

    def start_step(self):
        super().start_step()
        self.n_start = self.sim.people.uid.len_used
        self._new = []
        return

    def step_state(self):
        """ HIV death, progression off ART and CD4 moves on ART (one competing draw), for those with HIV before this step and not dying or leaving """
        super().step_state()
        sim, ppl, ti = self.sim, self.sim.people, self.ti
        au = ppl.auids
        au = au[self.infected.raw[au] & (au < self.n_start)]
        au = au[(self.ti_infected.raw[au] < ti) & ~leaving(ppl, au, ti)]
        arm, b = self.arm(au), self.band.raw[au]
        d = self.death[arm, b]
        surv = (1 - d)*bg_survival(sim, au)
        R = opp.REF_DT
        hd, hu, hn = d/R, self.up[arm, b]*surv/R, self.down[arm, b]*surv/R
        tot = hd + hu + hn
        p = -np.expm1(-tot*self.t.dt_year)
        x = self.rng.random(len(au))*np.divide(tot, p, out=np.zeros_like(tot), where=p > 0) # Where in the hazards each falls, if they leave
        go = x < tot
        dead = au[go & (x < hd)]
        upm = go & (x >= hd) & (x < hd + hu)
        dnm = go & (x >= hd + hu)
        self.band[au[upm]] = b[upm] + 1
        self.band[au[dnm]] = b[dnm] - 1
        if len(dead):
            ppl.request_death(dead, cause=self)
        self.dead_uids = ss.uids(dead)
        au = ppl.auids
        aids = au[self.infected.raw[au] & (self.band.raw[au] >= opp.AIDS) & ~self.symptomatic.raw[au]]
        self.symptomatic[aids] = True
        self.newly_symptomatic = aids
        return

    def step(self):
        super().step()
        self.new_uids = ss.uids(np.concatenate(self._new)) if self._new else ss.uids()
        return


class OptimaS(ss.Intervention):
    """ Optima's state in 1990: each population's shares of the 38 states, exactly in each population (with at least one agent in each state with people, as the declared engine) """
    def __init__(self, host=None, name='optima_s', seed_offset=2601, **kwargs):
        super().__init__(name=name, **kwargs)
        self.host_class, self.seed_offset = host, seed_offset
        return

    def init_post(self):
        super().init_post()
        sim = self.sim
        r = ref()
        names = [s.replace('-', '_') for s in r.settings['statelabels']]
        people0 = r.raw['people'][:, :, 0]
        shares = people0/people0.sum(axis=0) # States × pops
        rng = np.random.default_rng(sim.pars.rand_seed + self.seed_offset)
        uids = sim.people.auids
        pop = pop_of(sim, uids)
        assert (pop >= 0).all(), 'S needs each agent\'s Optima population (D)'
        plan, wanted, placed = {}, np.zeros(len(names)), np.zeros(len(names))
        for g in range(NP):
            inds = uids[pop == g]
            inds = inds[np.argsort(rng.random(len(inds)), kind='stable')]
            want = shares[:, g]*len(inds)
            n = round_preserving(want)
            plan[g] = (inds, want, n)
            wanted += want; placed += n
        for i in np.flatnonzero((wanted > 0) & (placed == 0)): # One agent instead of none, in the population that wanted the most, from its largest state
            g = max(plan, key=lambda g: plan[g][1][i])
            n = plan[g][2]
            if n.max() > 1:
                n[np.argmax(n)] -= 1
                n[i] += 1
        h = self.host_class(sim)
        dis = sim.diseases.hiv
        for g, (inds, _, n) in plan.items():
            edges = np.cumsum(np.append(0, n)).astype(int)
            for i, name in enumerate(names):
                who = ss.uids(inds[edges[i]:edges[i+1]])
                if not len(who) or name == 'susreg':
                    continue
                if name == 'progcirc':
                    h.circumcise(who)
                    continue
                c, k = name.split('_')
                if hasattr(dis, 'seed'):
                    dis.seed(who, opp.CD4.index(k))
                else:
                    dis.infect(who, source=None)
                if c != 'undx':
                    h.test(who)
                if c in ('care', 'lost'):
                    h.link_to_care(who)
                if c == 'lost':
                    h.lose_from_care(who)
                if c in opp.TREATED:
                    h.start_art(who, 'suppressive' if c == 'svl' else 'nonsuppressive')
        return

    def step(self):
        return


class ModulesHost(CascadeHost):
    """ `CascadeHost` on Optima′'s modules: as host.OptimaHost (its docstring), reading H's places and bands, D's populations, T's sources """
    def __init__(self, sim):
        super().__init__(sim)
        self.ppl = sim.people
        self.dis = sim.diseases.hiv
        return

    def pop(self, uids):
        return self.ppl.optima_pop.raw[uids]

    def place(self, uids):
        return self.dis.place.raw[uids]

    # People
    def age_years(self, uids):
        return self.ppl.age.raw[uids]

    def female(self, uids):
        return self.ppl.female.raw[uids]

    def pregnant(self, uids):
        return np.zeros(len(uids), bool)

    def post_debut(self, uids):
        return self.ppl.age.raw[uids] >= 15

    def dying(self, uids):
        return leaving(self.ppl, uids, self.sim.ti)

    def risk(self, uids):
        return np.where(self.pop(uids) <= pl.CLIENTS, 2, 1)

    # Events
    def debuts(self):
        d = modules(self.sim)[0]
        return d.debuts if d is not None else ss.uids()

    def newly_symptomatic(self):
        return self.dis.newly_symptomatic

    def twelve_weeks_pregnant(self):
        return ss.uids()

    def six_weeks_old(self):
        return ss.uids()

    # HIV
    def infected(self, uids):
        return self.dis.infected.raw[uids]

    def test(self, uids):
        uids = ss.uids(uids)
        pos = self.infected(uids)
        u = uids[self.place(uids) == UNDX]
        self.dis.place[u] = DX
        return pos

    def cd4(self, uids):
        return self.dis.cd4(uids)

    def record_cd4(self, uids):
        lc = self.dis.lowest_cd4
        lc[uids] = np.minimum(lc.raw[uids], self.cd4(uids))
        return lc.raw[uids]

    def who_stage(self, uids):
        c = self.cd4(uids)
        return np.select([c >= 500, c >= 350, c >= 200], [1, 2, 3], 4).astype(float)

    def on_art(self, uids):
        return self.dis.on_art(uids)

    def start_art(self, uids, arm):
        self.dis.start_art(uids, arm)

    def stop_art(self, uids):
        self.dis.stop_art(uids)

    def start_pmtct(self, uids, efficacy, days):
        if len(uids):
            raise NotImplementedError('Optima′ has no PMTCT per mother (its births have no mothers)')

    def circumcise(self, uids):
        uids = ss.uids(uids)
        men = uids[~self.female(uids) & ~self.infected(uids) & ~self.ppl.circumcised.raw[uids]]
        self.ppl.circumcised[men] = True
        return men

    def start_prep(self, uids, efficacy, days):
        if len(uids):
            raise NotImplementedError('Optima′ has no PrEP per person (a coverage by population in its force of infection)')

    # Places in care (common/cascade/pathways.py)
    def care_place(self, uids):
        return self.place(uids)

    def cd4_band(self, uids):
        return np.maximum(self.dis.band.raw[uids], 0)

    def population(self, uids):
        return np.array(pl.POPS)[self.pop(uids)]

    def rate_scale(self, uids):
        """ Optima′ moves people in its cascade among those who survive its step of 0.2 years: H's survival, times D's """
        return self.dis.hiv_survival(uids)*bg_survival(self.sim, uids)

    def link_to_care(self, uids):
        uids = ss.uids(uids)
        u = uids[np.isin(self.place(uids), [DX, LOST])]
        self.dis.place[u] = CARE

    def lose_from_care(self, uids):
        uids = ss.uids(uids)
        u = uids[self.place(uids) == CARE]
        self.dis.place[u] = LOST

    def switch_arm(self, uids, arm):
        self.dis.switch_arm(uids, arm)

    # Measurement
    def cascade_stage(self, uids):
        p = self.place(uids)
        inf = self.infected(uids)
        on = np.isin(p, [SVL, USVL])
        diagnosed = inf & (p != UNDX)
        ever = self.dis.ever_treated.raw[uids] | on
        category = np.select([~diagnosed, on, ever], ['undiagnosed', 'on_art', 'interrupted'], 'diagnosed_never_treated')
        category = np.where(inf, category, '')
        return dict(diagnosed=diagnosed, on_art=on, suppressive=p == SVL, ever_treated=ever, category=category)

    # For common/outputs.py
    def sexual_pairs(self):
        t = modules(self.sim)[1]
        return t.pairs if t is not None else (ss.uids(), ss.uids())

    def new_infections(self):
        return self.dis.new_uids

    def hiv_deaths(self):
        return self.dis.dead_uids
