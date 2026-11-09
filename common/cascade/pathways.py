"""
The cascade's pathways by rate and level (params.py: `rate`, `target`, `hold`; e.g. Optima's own cascade) for
agent-based models, as one Starsim intervention acting through a `CascadeHost` (interface.py). rates.py has their
compartmental form.

Each person has a place in care: undx, dx, care, lost, svl or usvl (Optima's). A host that keeps it gives it
(`care_place()`); otherwise it comes from the host's HIV (`infected`, and `cascade_stage()`'s diagnosed, on_art and
suppressive) and, for the diagnosed not on ART, this cascade's own record of who it put into care or lost.

Each step, after the host's own step:
- `rate` rows: everyone in a row's source places and group moves at its rate (or, given as people a year, that number
  over everyone in its `per` places), times the host's `rate_scale()`; a person's rows are competing risks (one move
  at most a step);
- `target` and `hold` rows: the number on ART filled from those in care (lowest CD4 band first before
  `lowest_cd4_before`), starting suppressed with `p_suppressive`, and removed back to care when over; shares held
  from their years, in the order diagnosed, in care, on ART, suppressed. As Optima (models/optima/v4/agents.py's
  `CascadeTargets`): the numbers from the counts (agents × people per agent), moved as each person's chance.

Moves through the host: undx → dx `test`; dx or lost → care `link_to_care`; care → lost `lose_from_care`; on ART →
lost `stop_art`; on ART → care `stop_art` then `link_to_care`; into svl or usvl `start_art` (arm suppressive or
nonsuppressive), or `switch_arm` if on ART. Coverages (pmtct, prep, pep rows) are not applied here: a model applies them
itself, except PMTCT with `pmtct=True`: at 12 weeks of pregnancy, a diagnosed woman not on ART starts PMTCT (`start_pmtct`,
the table's efficacy, to the end of pregnancy) with Optima's share, the table's number a year over the diagnosed women
giving birth a year (on ART or not, as Optima counts them; those pregnant now over the length of a pregnancy), at most 1.
"""
import numpy as np
import starsim as ss
from .params import parse_places
from .rates import BANDS, HOLDS

PLACES = ['', 'undx', 'dx', 'care', 'lost', 'svl', 'usvl']
UNDX, DX, CARE, LOST, SVL, USVL = range(1, 7)
ON_ART = [SVL, USVL]
ARM = {SVL: 'suppressive', USVL: 'nonsuppressive'}
EPS = 1e-3 # People (Optima's eps)
PREGNANCY = 280 # Days


class TableCascade(ss.Intervention):
    """
    Args:
        params (Params): the table
        host (class): the model's `CascadeHost`
        scale (float): people per agent, for the table's numbers (default: `sim.current_scale`)
        pmtct (bool): apply the table's PMTCT number per person (see above; default False, PMTCT left to the model)
    """
    def __init__(self, params, host, name='table_cascade', seed_offset=2311, scale=None, pmtct=False):
        super().__init__(name=name)
        self.params, self.host_class, self.seed_offset, self.scale = params, host, seed_offset, scale
        self.pmtct = pmtct and params.has('pmtct', 'number')
        P = params
        self.rates = []
        for r in P.rows('rate'):
            srcs, dst = parse_places(r['places'])
            self.rates.append(dict(name=r['parameter'], srcs=[PLACES.index(x) for x in srcs], dst=PLACES.index(dst), group=r['group'], fn=r['fn'],
                                   per=[PLACES.index(x) for x in r['per'].split('|')] if r['per'] else None))
        v = lambda pw, k, default: P.value(pw, k) if P.has(pw, k) else default
        self.target = P.has('target', 'on_art')
        self.lowest_before, self.p_supp = v('target', 'lowest_cd4_before', -np.inf), v('target', 'p_suppressive', 1.0)
        self.fix = {h: v('hold', h, np.inf) for h in HOLDS}
        self.held = {}
        self.needs = {k for r in self.rates for k in r['group']}
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.record = sim.people.add_state(ss.IntArr('cx_place', default=0)) # This cascade's record: care or lost, for hosts without care_place()
        self.host = self.host_class(sim)
        self.rng = np.random.default_rng(sim.pars.rand_seed + self.seed_offset)
        return

    # Places
    def places(self, uids):
        p = self.host.care_place(uids)
        if p is not None:
            return np.asarray(p)
        h = self.host
        inf = np.asarray(h.infected(uids), bool)
        st = h.cascade_stage(uids)
        rec = self.record.raw[uids]
        out = np.where(st['diagnosed'], np.where(np.isin(rec, [CARE, LOST]), rec, DX), UNDX)
        out = np.where(st['on_art'], np.where(st['suppressive'], SVL, USVL), out)
        return np.where(inf, out, 0)

    def groups(self, uids):
        h, g = self.host, {}
        if 'sex' in self.needs: g['sex'] = np.where(h.female(uids), 'f', 'm')
        if 'age' in self.needs: g['age'] = h.age_years(uids)
        if 'cd4' in self.needs: g['cd4'] = np.array(BANDS)[h.cd4_band(uids)]
        if 'pop' in self.needs: g['pop'] = np.asarray(h.population(uids))
        return g

    @staticmethod
    def in_group(group, g, n):
        m = np.ones(n, bool)
        for key, vals in group.items():
            if key == 'age':
                m &= np.any([(g['age'] >= float(v.split('-')[0])) & (g['age'] < float(v.split('-')[1])) for v in vals], axis=0)
            else:
                m &= np.isin(g[key], vals)
        return m

    # Moves through the host
    def act(self, uids, src, dst):
        if not len(uids):
            return
        h = self.host
        if dst == DX and src == UNDX:
            h.test(uids)
        elif dst == CARE and src in (DX, LOST):
            h.link_to_care(uids)
        elif dst == LOST and src == CARE:
            h.lose_from_care(uids)
        elif dst in (LOST, CARE) and src in ON_ART:
            h.stop_art(uids)
            if dst == CARE:
                h.link_to_care(uids)
        elif dst in ON_ART:
            h.switch_arm(uids, ARM[dst]) if src in ON_ART else h.start_art(uids, ARM[dst])
        else:
            raise NotImplementedError(f'No host move from {PLACES[src]} to {PLACES[dst]}')
        self.record[uids] = dst if dst in (CARE, LOST) else 0
        return

    def step(self):
        sim = self.sim
        self.year, dt, scale = self.t.now('year'), self.t.dt_year, self.scale or sim.current_scale
        au = sim.people.auids
        au = au[~np.asarray(self.host.dying(au), bool)]
        place = self.places(au)
        self.step_rates(au, place, dt, scale)
        self.step_levels(au, scale)
        if self.pmtct:
            self.step_pmtct(au, scale)
        return

    def step_pmtct(self, au, scale):
        """ Optima's PMTCT share for the diagnosed women at 12 weeks of pregnancy who are not on ART """
        h = self.host
        new = ss.uids(np.asarray(h.twelve_weeks_pregnant(), dtype=np.int64))
        if not len(new):
            return
        f = au[np.asarray(h.female(au), bool)]
        preg = f[np.asarray(h.pregnant(f), bool)]
        dx = preg[np.asarray(h.infected(preg), bool) & np.asarray(h.cascade_stage(preg)['diagnosed'], bool)]
        births = len(dx)*scale/(PREGNANCY/365) # Diagnosed women giving birth a year
        p = min(self.params.value('pmtct', 'number', self.year)/(EPS + births), 1)
        new = new[np.isin(new, au)]
        new = new[np.asarray(h.infected(new), bool)]
        st = h.cascade_stage(new)
        new = new[np.asarray(st['diagnosed'], bool) & ~np.asarray(st['on_art'], bool)]
        h.start_pmtct(new[self.rng.random(len(new)) < p], self.params.value('pmtct', 'efficacy'), PREGNANCY)
        return

    def step_rates(self, au, place, dt, scale):
        """ The rate rows, as competing risks """
        n = len(au)
        g = self.groups(au)
        rate = np.zeros((len(PLACES), n)) # By destination
        for r in self.rates:
            who = np.isin(place, r['srcs'])
            if r['group']:
                who &= self.in_group(r['group'], g, n)
            value = r['fn'](self.year)
            if r['per'] is not None:
                value = value/(EPS + np.isin(place, r['per']).sum()*scale)
            rate[r['dst'], who] += value
        rate *= self.host.rate_scale(au)[None, :]
        total = rate.sum(axis=0)
        u = self.rng.random((2, n))
        go = u[0] < -np.expm1(-total*dt)
        cum = np.cumsum(rate[:, go], axis=0)/total[go]
        dst = (u[1, go] > cum).sum(axis=0)
        moving, src = au[go], place[go]
        for s in np.unique(src):
            for d in np.unique(dst[src == s]):
                self.act(moving[(src == s) & (dst == d)], s, d)
        return

    def step_levels(self, au, scale):
        """ The target and hold rows (Optima's rules on levels) """
        if not self.target and not np.isfinite(list(self.fix.values())).any():
            return
        year, host = self.year, self.host
        band = None
        structure = dict(diagnosed=([UNDX], DX, [DX, CARE, LOST, SVL, USVL], [UNDX, DX, CARE, LOST, SVL, USVL]),
                         in_care=([DX, LOST], CARE, [CARE, SVL, USVL], [DX, CARE, LOST, SVL, USVL]),
                         on_art=([CARE], None, ON_ART, [CARE] + ON_ART), suppressed=([USVL], SVL, [SVL], ON_ART))
        for name, (lower, upper, numer, denom) in structure.items():
            if name not in self.held and year < self.fix[name] - 1e-9 and not (name == 'on_art' and self.target):
                continue
            place = self.places(au) # Counted again after each rule's moves
            count = lambda cs: np.isin(place, cs).sum()*scale
            if name not in self.held and year >= self.fix[name] - 1e-9:
                self.held[name] = count(numer)/(EPS + count(denom))
            actual, available = count(numer), count(denom)
            wanted = self.held[name]*available if name in self.held else self.params.value('target', 'on_art', year)
            diff = wanted - actual
            if diff > EPS: # Up: a share of each lower place (and CD4 band, lowest first for the target before its year)
                pool = np.isin(place, lower)
                total = pool.sum()*scale
                if total <= EPS: continue
                diff = min(diff, total - EPS)
                p = np.zeros(len(au))
                if name == 'on_art' and year < self.lowest_before:
                    band = host.cd4_band(au) if band is None else band
                    for k in reversed(range(len(BANDS))):
                        inb = pool & (band == k)
                        nk = inb.sum()*scale
                        if diff > EPS and nk > EPS:
                            m = min(diff, nk - EPS)
                            p[inb] = m/nk
                            diff -= m
                else:
                    p[pool] = diff/total
                u = self.rng.random(len(au))
                chosen = u < p
                if name == 'on_art':
                    supp = u < p*self.p_supp
                    for d, mask in ((SVL, chosen & supp), (USVL, chosen & ~supp)):
                        for s in lower:
                            self.act(au[mask & (place == s)], s, d)
                else:
                    for s in lower:
                        self.act(au[chosen & (place == s)], s, upper)
            elif diff < -EPS: # Down: a share of the upper place(s)
                srcs = ON_ART if name == 'on_art' else [upper]
                pool = np.isin(place, srcs)
                total = pool.sum()*scale
                if total <= EPS: continue
                chosen = pool & (self.rng.random(len(au)) < min(-diff, total - EPS)/total)
                dst = CARE if name == 'on_art' else LOST if name == 'in_care' else lower[0]
                for s in srcs:
                    self.act(au[chosen & (place == s)], s, dst)
            band = None # CD4 bands don't change with these moves, but recount anyway after a rule's moves
        return
