"""
EMOD′ sexual network: EMOD-HIV's pair-formation algorithm (PFA) and relationship types, as a Starsim v4 network.

What EMOD does (docs/emod_inventory.md, section 3; EMOD source in brackets), and so what this does, on each step:

1. Relationships whose timer has run out end (RelationshipManager::Update). A relationship lasts a Weibull duration,
   drawn at formation, rounded up to whole steps; one whose partner dies ends then.
2. Everyone past their debut is available for a type if, for every type they already have (or are queued for), they
   hold its extra-relational flag, and their relationships plus queue entries of that type are fewer than their maximum
   (IndividualHumanSTI::AvailableForRelationship). Flags are drawn along a chain (COMMERCIAL, TRANSITORY, INFORMAL,
   MARITAL), stopping at the first failure; maxima are randomly rounded. Both are drawn again when Risk changes
   (ConcurrencyConfiguration).
3. During the burn-in (PFA_Burnin_Duration_In_Days), the rate at which each available person joins a type's queue is
   set from the number available: the formation rate times all those available, half to each sex, split across age
   bins by the marginals of the joint age table, and divided among those in the bin, with 10 times the rate for those
   holding the type's flag (FlowControllerImpl). After the burn-in these rates are frozen.
4. Each available person joins queues by exponential waiting times at the sum of their rates over types, until the
   end of the step or until they are no longer available (IndividualHumanSTI::ConsiderRelationships). A woman stays in
   the age bin she joined in; queue entries last until matched or death.
5. Each man in the queue, in order of joining, weights the age bins that have women queued by his row of the joint age
   table (his current age); if the total is above PFA_Cum_Prob_Selection_Threshold, he picks a bin, then a woman in it
   by the Risk assortativity weights (BehaviorPfa::Update, Assortivity::FindPartnerIP). No match if no woman has a
   positive weight, or they are already partners; he stays queued.
6. Each relationship has Poisson coital acts in each step, at its type's rate times a dilution factor for the partner
   with more relationships, and at least one act in its first step (Relationship::Consummate).

Nodes are separate: people pair only within their node (`people.node`). `people.risk` (0 LOW, 1 MEDIUM, 2 HIGH) is
set by the risk module (risk.py). Edges: p1 is the man, p2 the woman.
"""
import numpy as np
import numba as nb
import pandas as pd
import sciris as sc
import starsim as ss

TYPES = ['TRANSITORY', 'INFORMAL', 'MARITAL', 'COMMERCIAL']
KEYS = ['T', 'I', 'M', 'C']
NT = len(TYPES)
MAX_SLOTS = 63 # The most relationships plus queue entries anyone can have (IRelationship.h)
BIG = np.int64(1) << 32 # For pair keys, man*BIG + woman


@nb.njit(cache=True)
def nb_seed(seed):
    np.random.seed(seed)


def available_all(a, q, mx, flags):
    """ Availability of each person (rows) for each type (columns) """
    eff = a + q
    bits = ((flags[:, None] >> np.arange(eff.shape[1])) & 1).astype(bool)
    ok = np.all((eff == 0) | bits, axis=1) & (eff.sum(axis=1) <= MAX_SLOTS)
    return ok[:, None] & (eff < mx)


@nb.njit(cache=True)
def to_array(lst):
    out = np.empty(len(lst), np.int64)
    for i in range(len(lst)):
        out[i] = lst[i]
    return out


@nb.njit(cache=True)
def available(a, q, mx, flags, t):
    """ Whether one person is available for type t, given their relationships (a) and queue entries (q) by type """
    total = 0
    for k in range(a.shape[0]):
        e = a[k] + q[k]
        total += e
        if e > 0 and not (flags >> k) & 1:
            return False
    return (a[t] + q[t] < mx[t]) and (total <= MAX_SLOTS)


@nb.njit(cache=True)
def nb_enter(rates, a, q, mx, flags, dt):
    """ Queue entries in one step: exponential waiting times at the summed rates of the types each person is available for; returns person index and type of each, in order """
    n, nt = rates.shape
    out_i, out_t = nb.typed.List.empty_list(nb.int64), nb.typed.List.empty_list(nb.int64)
    r = np.zeros(nt)
    for i in range(n):
        cum = 0.0
        for t in range(nt):
            r[t] = rates[i, t] if available(a[i], q[i], mx[i], flags[i], t) else 0.0
            cum += r[t]
        elapsed = 0.0
        while cum > 0 and elapsed < dt:
            elapsed += -np.log(1 - np.random.random())/cum
            if elapsed <= dt:
                u = np.random.random()*cum
                acc = 0.0
                t = 0
                for t in range(nt):
                    acc += r[t]
                    if acc > u:
                        break
                out_i.append(i)
                out_t.append(t)
                q[i, t] += 1
            cum = 0.0
            for t in range(nt):
                if r[t] > 0 and not available(a[i], q[i], mx[i], flags[i], t):
                    r[t] = 0.0
                cum += r[t]
    return to_array(out_i), to_array(out_t)


@nb.njit(cache=True)
def nb_match(m_uid, m_bin, m_risk, f_uid, f_bin, f_risk, joint, assort, threshold, keys):
    """
    Match one queue (one type, one node). Men in queue order; each picks a female age bin among those with women
    queued, by his row of the (row-normalized) joint table, if its total is above the threshold, then a woman in the
    bin with probability proportional to the assortativity weight of her Risk (drawn as a Risk level, weighted by the
    number queued at it, then a woman at random within it). `keys` are the existing pairs (man*BIG + woman), sorted.
    Returns the indices of the matched men and women entries.
    """
    nfb = joint.shape[1]
    nf, nm = len(f_uid), len(m_uid)
    cls = f_bin*3 + f_risk
    order = np.argsort(cls, kind='mergesort')
    count = np.zeros(nfb*3, np.int64)
    for c in cls:
        count[c] += 1
    start = np.zeros(nfb*3, np.int64)
    for c in range(1, nfb*3):
        start[c] = start[c-1] + count[c-1]
    slots = order.copy() # Live entries of class c are slots[start[c]:start[c] + count[c]]
    binsize = np.zeros(nfb, np.int64)
    for c in range(nfb*3):
        binsize[c//3] += count[c]
    pref = np.zeros(nfb)
    w = np.zeros(3)
    new = set()
    new.add(np.int64(-1))
    mi, fi = nb.typed.List.empty_list(nb.int64), nb.typed.List.empty_list(nb.int64)
    for i in range(nm):
        cum = 0.0
        for b in range(nfb):
            pref[b] = joint[m_bin[i], b] if binsize[b] > 0 else 0.0
            cum += pref[b]
        if cum <= threshold:
            continue
        u = np.random.random()*cum
        b, acc = 0, pref[0]
        while u > acc and b < nfb - 1:
            b += 1
            acc += pref[b]
        tot = 0.0
        for r in range(3):
            w[r] = count[b*3 + r]*assort[m_risk[i], r]
            tot += w[r]
        if tot <= 0:
            continue
        u = np.random.random()*tot
        r, acc = 0, w[0]
        while u > acc and r < 2:
            r += 1
            acc += w[r]
        c = b*3 + r
        k = start[c] + int(np.random.random()*count[c])
        j = slots[k]
        key = m_uid[i]*BIG + f_uid[j]
        pos = np.searchsorted(keys, key)
        if (pos < len(keys) and keys[pos] == key) or key in new:
            continue
        new.add(key)
        slots[k] = slots[start[c] + count[c] - 1]
        count[c] -= 1
        binsize[b] -= 1
        mi.append(i)
        fi.append(j)
    return to_array(mi), to_array(fi)


def sigmoid(c, year):
    """ EMOD's variable width and height sigmoid (condom use by year) """
    return c['Min'] + (c['Max'] - c['Min'])/(1 + np.exp(-c['Rate']*(year - c['Mid'])))


class EmodNetwork(ss.SexualNet):
    """
    EMOD-HIV's pair-formation algorithm, with its four relationship types (see the module docstring)

    Args:
        inp (objdict): EMOD's inputs for one parameter set, from `inputs.load()`
        record (bool): whether to keep a record of every relationship's start and end (`starts`, `ends`)
    """
    def __init__(self, inp, record=True, name='emodnet', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp = inp
        self.record = record
        self.define_states(
            *[ss.IntArr(f'n_{k}', default=0) for k in KEYS],   # Relationships of each type
            *[ss.IntArr(f'q_{k}', default=0) for k in KEYS],   # Queue entries of each type
            *[ss.IntArr(f'max_{k}', default=0) for k in KEYS], # Most relationships of each type at a time
            ss.IntArr('flags', default=0),       # Extra-relational flags, one bit per type
            ss.IntArr('risk_drawn', default=-1), # The Risk the flags and maxima were drawn for
        )
        self.meta.type = ss.dtypes.int
        self.meta.start = ss.dtypes.int  # Step formed
        self.meta.dur_days = ss.dtypes.float
        self.meta.total_acts = ss.dtypes.int
        self.edge_defaults.update(beta=1.0, acts=0.0, total_acts=0)
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.node = sim.people.add_state(ss.IntArr('node', default=0))
        self.risk = sim.people.add_state(ss.IntArr('risk', default=-1))
        inp = self.inp
        self.dt_days = inp.dt
        self.nbins = [len(t.edges[0]) for t in inp.types]
        assert all(len(t.edges[0]) == len(t.edges[1]) for t in inp.types)
        self.rates = np.zeros((3, NT, 2, max(self.nbins), 2)) # Per day: node, type, sex, age bin, has the type's flag
        self.burnin_left = inp.burnin
        self.queue = sc.objdict(uid=np.zeros(0, np.int64), type=np.zeros(0, np.int64), bin=np.zeros(0, np.int64)) # In order of joining
        self.starts, self.ends = [], []
        self.rng = np.random.default_rng(sim.pars.rand_seed + 1)
        nb_seed(sim.pars.rand_seed + 2)
        return

    def init_results(self):
        super().init_results()
        self.define_results(*[ss.Result(f'new_{k}', dtype=int, scale=True, label=f'New {t.lower()} relationships') for k, t in zip(KEYS, TYPES)],
                            *[ss.Result(f'n_{k}', dtype=int, scale=True, label=f'{t.title()} relationships') for k, t in zip(KEYS, TYPES)])
        return

    # %% People

    def stack(self, prefix, uids):
        """ Per-type states for some people, as an array of people × type """
        return np.stack([getattr(self, f'{prefix}_{k}').raw[uids] for k in KEYS], axis=1)

    def set_people(self):
        """ Debut ages for those with none; flags and maxima for those whose Risk has changed (or who have none) """
        inp, ppl = self.inp, self.sim.people
        auids = ppl.auids
        uids = auids[np.isnan(self.debut.raw[auids])]
        if len(uids):
            sex = ppl.female.raw[uids].astype(int)
            draw = np.array(inp.debut.scale)[sex]*(-np.log(self.rng.random(len(uids))))**np.array(inp.debut.het)[sex]
            self.debut[uids] = np.maximum(inp.debut.min, draw)
        risk = self.risk.raw[auids]
        uids = auids[(risk >= 0) & (risk != self.risk_drawn.raw[auids])]
        if len(uids):
            risk, sex = self.risk.raw[uids], ppl.female.raw[uids].astype(int)
            flags, going = np.zeros(len(uids), int), np.ones(len(uids), bool)
            p_extra = np.array([typ.p_extra for typ in inp.types]) # Type × risk × sex
            for pos in range(NT): # Along the chain, stopping at the first failure
                t = np.array(inp.flag_order)[risk, pos]
                p = p_extra[t, risk, sex]
                got = going & (self.rng.random(len(uids)) < p)
                flags[got] |= 1 << t[got]
                going &= got
            self.flags[uids] = flags
            for k, typ in zip(KEYS, inp.types):
                m = typ.max_rels[risk, sex]
                getattr(self, f'max_{k}')[uids] = np.floor(m).astype(int) + (self.rng.random(len(uids)) < m - np.floor(m))
            self.risk_drawn[uids] = risk
        return

    def bins(self, uids):
        """ Each person's age bin for each type (people × type), as EMOD's BinIndexForAgeAndSex """
        ppl = self.sim.people
        age, female = ppl.age.raw[uids], ppl.female.raw[uids]
        out = np.zeros((len(uids), NT), np.int64)
        for t, typ in enumerate(self.inp.types):
            for s in (0, 1):
                m = female == s
                out[m, t] = np.minimum(np.searchsorted(typ.edges[s], age[m], side='right'), len(typ.edges[s]) - 1)
        return out

    # %% Steps

    def step(self):
        self.end_pairs()
        self.set_people()
        self.form()
        self.consummate()
        return

    def end_pairs(self):
        """ End relationships whose time is up """
        if not len(self):
            return
        self.edges.dur = self.edges.dur - 1
        self.drop(self.edges.dur <= 0, died=False)
        return

    def drop(self, gone, died):
        """ Remove the edges in a boolean mask, updating counts and the record """
        if not gone.any():
            return
        e = self.edges
        for t, k in enumerate(KEYS):
            m = gone & (e.type == t)
            for p in (e.p1[m], e.p2[m]):
                np.add.at(getattr(self, f'n_{k}').raw, p, -1)
        if self.record:
            self.ends.append(dict(type=e.type[gone], start=e.start[gone], end=np.full(gone.sum(), self.ti), acts=e.total_acts[gone], died=np.full(gone.sum(), died)))
        self.remove_edges(gone)
        return

    def remove_uids(self, uids):
        """ Relationships of people who die end """
        if len(self):
            is_gone = np.zeros(len(self.sim.people.alive.raw), bool)
            is_gone[uids] = True
            self.drop(is_gone[self.edges.p1] | is_gone[self.edges.p2], died=True)
        q = self.queue
        keep = ~np.isin(q.uid, uids)
        for k in q.keys():
            q[k] = q[k][keep]
        return

    def form(self):
        """ Eligibility, queue entry (rates set during the burn-in) and matching """
        inp, ppl = self.inp, self.sim.people
        auids = ppl.auids
        uids = auids[(self.risk.raw[auids] >= 0) & (ppl.age.raw[auids] >= self.debut.raw[auids])]
        a, q, mx = self.stack('n', uids), self.stack('q', uids), self.stack('max', uids)
        flags = self.flags.raw[uids].astype(np.int64)
        node, female = self.node.raw[uids], ppl.female.raw[uids].astype(int)
        bins = self.bins(uids)
        has_flag = (flags[:, None] >> np.arange(NT)) & 1

        if self.burnin_left > 0:
            self.burnin_left -= self.dt_days
            avail = available_all(a, q, mx, flags)
            for n in range(3):
                for t, typ in enumerate(inp.types):
                    m = avail[:, t] & (node == n)
                    nb_ = self.nbins[t]
                    counts = np.zeros((2, nb_, 2)) # Sex, bin, flag
                    np.add.at(counts, (female[m], bins[m, t], has_flag[m, t]), 1)
                    desired = 0.5*typ.rate*counts.sum()*np.array(typ.marginal) # Sex × bin
                    eff = np.array(typ.rate_ratio)[:, None]*counts[:, :, 1] + counts[:, :, 0]
                    low = np.divide(desired, eff, out=np.zeros_like(desired), where=eff > 0)
                    self.rates[n, t, :, :nb_, 0] = low
                    self.rates[n, t, :, :nb_, 1] = np.array(typ.rate_ratio)[:, None]*low

        rates = self.rates[node[:, None], np.arange(NT)[None, :], female[:, None], bins, has_flag]
        ei, et = nb_enter(rates, a, q, mx, flags, self.dt_days)
        for t, k in enumerate(KEYS):
            getattr(self, f'q_{k}').raw[uids] = q[:, t]
        Q = self.queue
        Q.uid = np.concatenate([Q.uid, uids[ei]])
        Q.type = np.concatenate([Q.type, et])
        Q.bin = np.concatenate([Q.bin, bins[ei, et]])

        for t, typ in enumerate(inp.types):
            self.match(t, typ)
        return

    def match(self, t, typ):
        """ Match the queue of one type, node by node, and start the relationships """
        ppl, Q = self.sim.people, self.queue
        female, node, risk = ppl.female.raw, self.node.raw, self.risk.raw
        keys = np.sort(self.edges.p1.astype(np.int64)*BIG + self.edges.p2) if len(self) else np.zeros(0, np.int64)
        taken = np.zeros(len(Q.uid), bool)
        p1, p2 = [], []
        for n in range(3):
            inq = (Q.type == t) & (node[Q.uid] == n)
            men, women = np.flatnonzero(inq & ~female[Q.uid]), np.flatnonzero(inq & female[Q.uid])
            if not len(men) or not len(women):
                continue
            mu = Q.uid[men]
            mbin = np.minimum(np.searchsorted(typ.edges[0], ppl.age.raw[mu], side='right'), len(typ.edges[0]) - 1) # His current age
            mi, fi = nb_match(mu, mbin, risk[mu], Q.uid[women], Q.bin[women], risk[Q.uid[women]], typ.joint, typ.assort, self.inp.threshold, keys)
            taken[men[mi]] = taken[women[fi]] = True
            p1.append(mu[mi])
            p2.append(Q.uid[women[fi]])
        if not p1:
            return
        p1, p2 = ss.uids(np.concatenate(p1)), ss.uids(np.concatenate(p2))
        for k in Q.keys():
            Q[k] = Q[k][~taken]
        k = KEYS[t]
        for arr in (getattr(self, f'q_{k}'), getattr(self, f'n_{k}')):
            for p in (p1, p2):
                np.add.at(arr.raw, p, -1 if arr.name.startswith('q') else 1)
        days = 365*typ.dur_scale*(-np.log(self.rng.random(len(p1))))**typ.dur_het
        dur = np.maximum(1, np.ceil(days/self.dt_days)) # Steps, including this one
        self.append(p1=p1, p2=p2, dur=dur, dur_days=days, type=np.full(len(p1), t), start=np.full(len(p1), self.ti))
        self.results[f'new_{k}'][self.ti] = len(p1)
        if self.record:
            age = ppl.age.raw
            self.starts.append(dict(type=np.full(len(p1), t), ti=np.full(len(p1), self.ti), node=node[p1], age_m=age[p1], age_f=age[p2],
                                    risk_m=risk[p1], risk_f=risk[p2], dur_days=days, p1=p1, p2=p2))
        return

    def consummate(self):
        """ Coital acts in each relationship in this step, and the share of them with a condom """
        if not len(self):
            return
        e, inp = self.edges, self.inp
        nrels = sum(getattr(self, f'n_{k}').raw for k in KEYS)
        most = np.minimum(np.maximum(nrels[e.p1], nrels[e.p2]), 4)
        rate = np.array([typ.coital_rate for typ in inp.types])[e.type]
        lam = self.dt_days*np.array(inp.dilution)[most]*rate
        acts = self.rng.poisson(lam)
        acts[(e.total_acts == 0) & (acts == 0) & (lam > 0)] = 1 # The first act is forced
        e.total_acts = e.total_acts + acts
        e.acts = acts.astype(float)
        year = self.t.now('year')
        cond = np.array([[sigmoid(c, year) for c in typ.condoms] for typ in inp.types]) # Type × node
        e.condoms = cond[e.type, self.node.raw[e.p1]]
        return

    def update_results(self):
        super().update_results()
        for t, k in enumerate(KEYS):
            self.results[f'n_{k}'][self.ti] = np.count_nonzero(self.edges.type == t) if len(self) else 0
        return

    def finalize(self):
        super().finalize()
        cat = lambda recs: pd.DataFrame({k: np.concatenate([r[k] for r in recs]) for k in recs[0]}) if recs else pd.DataFrame()
        self.starts, self.ends = cat(self.starts), cat(self.ends)
        return
