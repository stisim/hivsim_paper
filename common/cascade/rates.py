"""
The condensed cascade (a parameter table, params.py) as rates between the care states of a compartmental model, by
CD4 band, population group and year. See README.md for the mapping and where it is exact or approximate.

Care states: undx (undiagnosed), dx (diagnosed, not in care: the cascade's none, uptake or testing), care (staging,
pre-ART), lost (dropped out of ART: the cascade's out and lost), svl and usvl (on ART, in the suppressive and the
non-suppressive arm). Groups are dicts with `sex` ('f' or 'm'), `age_lo` and `age_hi` (years).

    cx = CompartmentalCascade(Params(table, flows), groups, years)
    cx.rates[('undx', 'dx')]               # per year, bands × groups × years
    cx.on_progression[('undx', 'care')]    # probability, on entering the symptomatic band, groups × years
    cx.births['care']                      # probability for an infected birth (the infant test), groups × years
    cx.mtct_averted                        # share of mother-to-child transmission averted by ANC testing and PMTCT, for mothers in undx or dx, years

Pathways by rate and level (params.py: `rate`, `target`, `hold`; e.g. Optima's own cascade), when the table has them:

    cx.path_rates[('undx', 'dx')]          # per year, bands × groups × years: the `rate` rows (rows on the same move add)
    cx.path_numbers                        # `rate` rows given as people a year: (sources, destination, per places, people a year by year)
    cx.targets                             # the `target` and `hold` rows (see TARGETS), for a rule applied to the levels each step
    cx.coverage('prep')                    # a coverage by group, groups × years (prep, pep), and cx.value('pmtct', 'number', year) etc.
"""
import numpy as np

STATES = ['undx', 'dx', 'care', 'lost', 'svl', 'usvl']
ARMS = dict(svl='suppressive', usvl='nonsuppressive')
OPTIMA_BANDS = [(np.inf, np.inf), (500, np.inf), (350, 500), (200, 350), (50, 200), (0, 50)] # CD4 bounds: acute (at the post-infection CD4), >500, ..., <50
EMOD_NH = dict(cd4_post=540.55, cd4_death=31.63, who_kappa=[0.9664, 0.9916, 0.9356], who_lambda=[0.26596, 0.19729, 0.34721]) # EMOD-HIV's CD4 and WHO stages (models/emod/v4/hiv.py)
GESTATION = 280/365
BANDS = ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50'] # Names of OPTIMA_BANDS, as a table's `cd4` group
PLACES = ['undx', 'dx', 'care', 'lost', 'svl', 'usvl'] # Places in care of rate pathways (STATES)
HOLDS = ['diagnosed', 'in_care', 'on_art', 'suppressed'] # Shares a `hold` row keeps, in the order applied


def group_weight(group, g):
    """ The share of population group `g` (a dict: sex, age_lo, age_hi, and optionally pop) in a row's group (params.parse_group), ages spread evenly """
    w = 1.0
    for key, vals in group.items():
        if key == 'sex':
            w *= float(g['sex'] in vals)
        elif key == 'pop':
            w *= float(g.get('pop') in vals)
        elif key == 'age':
            share = 0
            for v in vals:
                lo, hi = (float(x) for x in v.split('-'))
                share += max(0, min(hi, g['age_hi']) - max(lo, g['age_lo']))/(g['age_hi'] - g['age_lo'])
            w *= min(share, 1)
        elif key != 'cd4':
            raise ValueError(f'Unknown group key {key}')
    return w


def band_mask(group, bands=BANDS):
    """ Which CD4 bands a row's group covers """
    return np.array([('cd4' not in group) or (b in group['cd4']) for b in bands], dtype=float)


class Eligibility:
    """
    The share of each CD4 band meeting a WHO-stage and CD4 rule, from a natural history in which √CD4 falls linearly
    over the prognosis and WHO stages are shares of it (EMOD's): people spread evenly in time over the band, with
    WHO-stage shares drawn as EMOD draws them. Monte Carlo, fixed seed.
    """
    def __init__(self, bands, nh=EMOD_NH, n=20000, seed=1):
        rng = np.random.default_rng(seed)
        sp, sd = np.sqrt(nh['cd4_post']), np.sqrt(nh['cd4_death'])
        frac = lambda cd4: np.clip((sp - np.sqrt(np.minimum(cd4, nh['cd4_post'])))/(sp - sd), 0, 1)
        self.cd4, self.who = [], []
        for lo, hi in bands:
            f = rng.uniform(frac(hi), frac(lo), n)
            w = np.array([lam*(-np.log(rng.random(n)))**(1/k) for k, lam in zip(nh['who_kappa'], nh['who_lambda'])]).T
            self.cd4.append((sp - f*(sp - sd))**2)
            self.who.append(1 + (f[:, None] > np.cumsum(w, axis=1)).sum(axis=1)) # Whole WHO stage, 1-4
        self.cache = {}

    def share(self, band, who_min, cd4_max):
        key = (band, who_min, cd4_max)
        if key not in self.cache:
            self.cache[key] = float(np.mean((self.who[band] >= who_min) | (self.cd4[band] <= cd4_max)))
        return self.cache[key]


def integrate(deriv, x0, years, start, step=0.01):
    """ A scalar ODE from `start` (x0 then, and before), at the given years """
    out = np.full(len(years), float(x0))
    x, t = float(x0), float(start)
    for i in np.argsort(years):
        while t < years[i] - 1e-12:
            h = min(step, years[i] - t)
            x = float(np.clip(x + h*deriv(t, x), 0, 1))
            t += h
        if years[i] >= start:
            out[i] = x
    return out


class CompartmentalCascade:
    """
    Args:
        params (Params): the cascade's table and flows
        groups (list): one dict per group: `sex`, `age_lo`, `age_hi`
        years (array): the years to give rates at
        bands (list): CD4 bounds of each band (OPTIMA_BANDS: Optima's); `symptom_band` is the band whose entry is the
            first symptoms (Optima's CD4 50-200, its AIDS)
        pregnancy (array): pregnancies per woman per year, groups × years (e.g. the birth rate); none if not given
        post_debut (array): the share of each group past sexual debut (reached by routine testing); default 1 from 15
        other_exit (float): rate per year at which the undiagnosed and diagnosed leave their pool other than by routine
            testing (death, other testing, ART), for the share of them in the testing loop
        debut_turnover (float): new debuts a year per person past debut, for the share of everyone in the testing loop
        lost_exit (float): death rate per year of those dropped out of ART, for the share of them still out (not lost for good)
    """
    def __init__(self, params, groups, years, bands=OPTIMA_BANDS, symptom_band=4, pregnancy=None, post_debut=None,
                 nh=EMOD_NH, other_exit=0.1, debut_turnover=0.05, lost_exit=0.1):
        self.P, self.groups, self.years = params, groups, np.asarray(years, dtype=float)
        self.bands, self.symptom_band = bands, symptom_band
        G, T = len(groups), len(self.years)
        self.sex = np.array([g['sex'] for g in groups])
        self.lo, self.hi = (np.array([float(g[k]) for g in groups]) for k in ('age_lo', 'age_hi'))
        self.pregnancy = np.zeros((G, T)) if pregnancy is None else np.asarray(pregnancy, dtype=float)*(self.sex == 'f')[:, None]
        self.post_debut = (self.lo >= 15).astype(float) if post_debut is None else np.asarray(post_debut, dtype=float)
        self.rates, self.on_progression, self.births, self.mtct_averted = {}, {}, {}, np.zeros(T)
        if params.has('hct'): # The condensed cascade's pathways
            self.elig = Eligibility(bands, nh)
            self.hidden_states(other_exit, debut_turnover, lost_exit)
            self.build()
        self.build_pathways()

    # Table values
    def p(self, pathway, parameter, year, sex='all'):
        return self.P.p(pathway, parameter, year, sex)

    def v(self, pathway, parameter, year=None, default=None):
        return self.P.value(pathway, parameter, year) if self.P.has(pathway, parameter) else default

    def started(self, pathway, year):
        return year >= self.v(pathway, 'start', default=-np.inf)

    def mean_days(self, pathway, parameter):
        form, v = self.P.delay(pathway, parameter)
        return (v[0] + v[1])/2 if form == 'uniform_days' else v[0]

    def hidden_states(self, other_exit, turnover, lost_exit):
        """
        The cascade's states that the care states don't hold, as shares over time: of everyone past debut, of the
        undiagnosed and of the diagnosed in the routine testing loop (the rest wait in uptake), and of those dropped
        out of ART who are still out rather than lost for good. The testing pools are taken as steady in size (inflow = outflow).
        """
        P, yrs = self.P, self.years
        tt, tu, td = self.mean_days('hct', 'test_interval')/365, self.mean_days('hct', 'uptake_interval')/365, self.mean_days('art', 'dropout_delay')/365
        a = lambda y: self.p('hct', 'p_uptake', y)/tu # From uptake into the loop
        start = self.v('hct', 'start', default=yrs.min())
        fine = np.arange(start, yrs.max() + 1e-9, 0.05)
        s = integrate(lambda y, s: a(y)*(1 - s) + turnover*(self.p('hct', 'p_debut', y) - s), 0, fine, start)
        sgen = lambda y: np.interp(y, fine, s)
        pool = lambda sin: (lambda y, x: (x/tt + other_exit)*(sin(y) - x) + a(y)*(1 - x) - x*(1 - x)/tt)
        self.loop_share = dict(everyone=np.interp(yrs, fine, s, left=0),
                               undx=integrate(pool(sgen), 0, yrs, start), dx=integrate(pool(lambda y: 0), 0, yrs, start))
        # Dropped out of ART: out (offered dropout and return again) or lost for good. Of each dropout, i_out go out
        # and i_lost are lost; those out leave at `leave` (to lost for good at p_lost/td). Dropouts at a steady rate
        # from the first year anyone can start ART; both die at lost_exit
        def drop(y):
            pl, po, pr = (self.p('dropout', k, y) for k in ('p_lost', 'p_out', 'p_restart'))
            return pl, po, pr, (pl + (1 - pl)*(1 - po) + (1 - pl)*po*pr)/td
        self.drop = drop
        grid = np.arange(min(yrs.min(), start), yrs.max() + 1e-9, 0.01)
        art = next((y for y in grid if self.v('eligibility', 'cd4', y, -np.inf) > 0 or max(self.p('hct', 'p_immediate_art', y, s) for s in 'fm') > 0), np.inf)
        o = f = 0.0
        phi = np.zeros(len(grid))
        for i, y in enumerate(grid):
            pl, po, pr, leave = drop(y)
            i_out, i_lost = (1 - pl)*po*(1 - pr), pl
            if y >= art:
                o, f = o + 0.01*(i_out - (leave + lost_exit)*o), f + 0.01*(i_lost + pl/td*o - lost_exit*f)
            phi[i] = o/(o + f) if o + f > 0 else i_out/max(i_out + i_lost, 1e-12)
        self.out_share = np.interp(yrs, grid, phi)

    def eligible(self, band, g, year, pregnant):
        """ The share of band `band` of group `g` eligible for ART at staging: the eligibility rule of `year` """
        v = lambda name, default: self.v('eligibility', name, year, default)
        cd4 = v('cd4', -np.inf)
        if pregnant:
            cd4 = max(cd4, v('cd4_pregnant', -np.inf))
        if self.lo[g] >= v('adult_age', 5):
            if pregnant and v('pregnant', 0) > 0:
                return 1.0
            return self.elig.share(band, v('who_adult', 6), cd4)
        child = np.clip((v('child_age', -1) - self.lo[g])/(self.hi[g] - self.lo[g]), 0, 1) # Ages spread evenly over the group
        return child + (1 - child)*self.elig.share(band, v('who_child', 6), cd4)

    def eligible_share(self, band, g, year):
        """ As eligible(), with the group's share pregnant (from its pregnancy rate) """
        pi = min(self.pregnancy[g, self.ti]*GESTATION, 1)
        return pi*self.eligible(band, g, year, True) + (1 - pi)*self.eligible(band, g, year, False)

    def arms(self, g, year, age=None):
        """ Probabilities of the suppressive (svl) and non-suppressive (usvl) arm at an ART start """
        sex, split = self.sex[g], self.v('art', 'age_split')
        under = float(age < split) if age is not None else np.clip((split - self.lo[g])/(self.hi[g] - self.lo[g]), 0, 1)
        ps = under*self.p('art', 'p_suppressive_under', year, sex) + (1 - under)*self.p('art', 'p_suppressive_over', year, sex)
        return dict(svl=ps, usvl=1 - ps)

    def staging(self, band, g, year, pregnant=False, elig=None):
        """ Where someone sent to staging goes: {state: probability}, the rest back to where they were (dx, or undx if not tested) """
        if elig is None:
            elig = self.eligible(band, g, year, True) if pregnant else self.eligible_share(band, g, year)
        kept = 1 - self.p('staging', 'p_lost', year)
        link = elig*self.p('staging', 'p_link_art', year, self.sex[g]) # ART after the start delay, taken as at once
        out = {arm: kept*link*q for arm, q in self.arms(g, year).items()}
        out['care'] = kept*(1 - elig)*self.p('staging', 'p_link_preart', year)
        return out

    def build(self):
        B, G, T = len(self.bands), len(self.groups), len(self.years)
        R = {}
        def add(src, dst, k, g, value):
            R.setdefault((src, dst), np.zeros((B, G, T)))[k, g, self.ti] += value
        prog = {}
        births = {s: np.zeros((G, T)) for s in ['dx', 'care', 'svl', 'usvl']}
        self.mtct_averted = np.zeros(T)
        tt = self.mean_days('hct', 'test_interval')/365
        interval, td = self.mean_days('preart', 'interval')/365, self.mean_days('art', 'dropout_delay')/365
        for self.ti, y in enumerate(self.years):
            hct, anc, sym = self.started('hct', y), self.started('anc', y), self.started('symptomatic', y)
            pl, po, pr, _ = self.drop(y)
            for g in range(G):
                sex = self.sex[g]
                arms = self.arms(g, y)
                rate_anc = self.pregnancy[g, self.ti]*self.p('anc', 'p_test', y) if anc else 0
                p_imm = self.p('hct', 'p_immediate_art', y, sex)
                p_stage = self.p('anc', 'p_stage', y)
                for k in range(B):
                    st = self.staging(k, g, y)
                    st_anc = self.staging(k, g, y, pregnant=True)
                    for src in ('undx', 'dx'):
                        rate_hct = self.post_debut[g]*self.loop_share[src][self.ti]/tt if hct else 0
                        for dst in ('svl', 'usvl', 'care'):
                            add(src, dst, k, g, rate_hct*((p_imm*arms[dst] if dst != 'care' else 0) + (1 - p_imm)*st[dst]) + rate_anc*p_stage*st_anc[dst])
                        if src == 'undx': # A positive test diagnoses
                            add(src, 'dx', k, g, rate_hct*(1 - p_imm)*(1 - sum(st.values())) + rate_anc*(1 - p_stage*sum(st_anc.values())))
                    # Pre-ART: a check each interval, lost (back to dx) or else the rule again
                    add('care', 'dx', k, g, self.p('preart', 'p_lost', y)/interval)
                    el = self.eligible_share(k, g, y)
                    for arm, q in arms.items():
                        add('care', arm, k, g, (1 - self.p('preart', 'p_lost', y))*el*q/interval)
                    # ART: dropout, then lost for good, out (and restarting, with the arm drawn again), or back to dx
                    for arm, other in (('svl', 'usvl'), ('usvl', 'svl')):
                        add(arm, 'lost', k, g, (pl + (1 - pl)*po*(1 - pr))/td)
                        add(arm, 'dx', k, g, (1 - pl)*(1 - po)/td)
                        add(arm, other, k, g, (1 - pl)*po*pr*arms[other]/td)
                    phi = self.out_share[self.ti]/td # Only those out (not lost for good) are offered return
                    add('lost', 'dx', k, g, phi*(1 - pl)*(1 - po))
                    for arm, q in arms.items():
                        add('lost', arm, k, g, phi*(1 - pl)*po*pr*q)
                # First symptoms, on entering the symptomatic band: a test (to staging), or else staging untested
                k = self.symptom_band
                st = self.staging(k, g, y)
                tested, again = (self.p('symptomatic', 'p_present', y), (1 - self.p('symptomatic', 'p_present', y))*self.p('symptomatic', 'p_present_again', y)) if sym else (0, 0)
                for src in ('undx', 'dx'):
                    for dst, q in st.items():
                        prog.setdefault((src, dst), np.zeros((G, T)))[g, self.ti] = (tested + again)*q
                    if src == 'undx':
                        prog.setdefault(('undx', 'dx'), np.zeros((G, T)))[g, self.ti] = tested*(1 - sum(st.values()))
                # Infected births: the test at 6 weeks (a child, under the age split); positive: staging
                if self.started('infant', y):
                    v = lambda name, default: self.v('eligibility', name, y, default)
                    infant = 1.0 if v('child_age', -1) >= 6/52 else self.elig.share(0, v('who_child', 6), v('cd4', -np.inf))
                    p_test, kept = self.p('infant', 'p_test', y), 1 - self.p('staging', 'p_lost', y)
                    link = infant*self.p('staging', 'p_link_art', y, sex)
                    for arm, q in self.arms(g, y, age=0).items():
                        births[arm][g, self.ti] = p_test*kept*link*q
                    births['care'][g, self.ti] = p_test*kept*(1 - infant)*self.p('staging', 'p_link_preart', y)
                    births['dx'][g, self.ti] = p_test*(1 - kept*(link + (1 - infant)*self.p('staging', 'p_link_preart', y)))
            # PMTCT after a positive ANC test (for mothers not in care): combination (option B or A), or single-dose
            # nevirapine; with blocked_if_staged, those sent to staging get none of the combination
            if anc:
                eff = lambda name: self.v('pmtct', name)
                comb, b, staged = self.p('pmtct', 'p_combination', y), self.p('pmtct', 'p_option_b', y), self.p('anc', 'p_stage', y)
                e_comb = b*eff('eff_option_b') + (1 - b)*eff('eff_option_a')
                blocked = staged*(eff('blocked_if_staged') > 0)
                self.mtct_averted[self.ti] = self.p('anc', 'p_test', y)*((1 - blocked)*comb*e_comb + (1 - comb)*eff('eff_sdnvp'))
        self.rates, self.on_progression, self.births = R, prog, births
        del self.ti

    @property
    def has_pathways(self):
        """ Whether the table has pathways by rate and level (Optima's kind) """
        return bool(self.path_rates or self.path_numbers or self.targets)

    def build_pathways(self):
        """ The `rate` rows as rates by band × group × year (or numbers a year), and the `target` and `hold` rows """
        B, G, T = len(self.bands), len(self.groups), len(self.years)
        self.path_rates, self.path_numbers = {}, []
        for r in self.P.rows('rate'):
            srcs, dst = r['places'].split('>')
            values = np.array([r['fn'](y) for y in self.years])
            if r['per']:
                assert not r['group'], 'A rate given as a number is for everyone in its places'
                self.path_numbers.append((srcs.split('|'), dst, r['per'].split('|'), values))
                continue
            w = np.array([group_weight(r['group'], g) for g in self.groups])
            for src in srcs.split('|'):
                arr = self.path_rates.setdefault((src, dst), np.zeros((B, G, T)))
                arr += band_mask(r['group'], BANDS[:B])[:, None, None]*w[None, :, None]*values[None, None, :]
        t = {}
        if self.P.has('target', 'on_art'):
            t.update(on_art=lambda y: self.P.value('target', 'on_art', y), lowest_cd4_before=self.v('target', 'lowest_cd4_before', default=-np.inf),
                     p_suppressive=self.v('target', 'p_suppressive', default=1.0))
        hold = {h: self.v('hold', h) for h in HOLDS if self.P.has('hold', h)}
        if hold:
            t['hold'] = hold
        self.targets = t

    def coverage(self, pathway, parameter='coverage'):
        """ A coverage by group from the rows by group (missing groups 0), groups × years """
        out = np.zeros((len(self.groups), len(self.years)))
        if self.P.has(pathway, parameter) and not any(r['parameter'] == parameter for r in self.P.rows(pathway)): # One row for everyone
            return out + np.array([self.P.value(pathway, parameter, y) for y in self.years])[None, :]
        for r in self.P.rows(pathway):
            if r['parameter'] == parameter:
                out += np.array([group_weight(r['group'], g) for g in self.groups])[:, None]*np.array([r['fn'](y) for y in self.years])[None, :]
        return out

    def flow_rows(self, intervention, year):
        """ The flows (VMMC, PrEP) of `year`: (sex, age_lo, age_hi, people a year at full scale, restrict) """
        f = self.P.flows
        if f is None:
            return []
        f = f[(f.intervention == intervention) & (f.start <= year) & (year < f.end)]
        return [(r.sex, r.age_lo, r.age_hi, r.n/(r.end - r.start), r.restrict) for r in f.itertuples()]
