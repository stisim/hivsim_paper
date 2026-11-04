"""
Optima′: Optima HIV (the Malawi project) as a declared Starsim model.

Each rate is Optima's probability per step of 0.2 years divided by 0.2, with Optima's products, so that Optima's
handling of time is kept (see docs/optima_inventory.md). Inputs are Optima's own series at dt=0.2.
"""
from pathlib import Path
import numpy as np, sciris as sc, starsim as ss

REF_DT = 0.2 # The step Optima's Malawi project was calibrated at
REFFILE = Path(__file__).resolve().parents[1] / 'reference' / 'optima_malawi_ref.obj'
CD4 = ['acute', 'gt500', 'gt350', 'gt200', 'gt50', 'lt50']
UNTREATED = ['undx', 'dx', 'care', 'lost'] # Cascade states not on ART
OPTIONS = sc.objdict(age_skip=True, age_scale=1.0, risk=True, skip_scale=1.0) # For checks: Optima's ageing band-skip, and a factor on ageing rates
TREATED = ['usvl', 'svl']
NEXT = dict(undx='dx', dx='care', care='lost', lost='care') # Where each untreated cascade state moves to
AIDS = CD4.index('gt50') # CD4 stages from here on are AIDS (Optima's settings.aidsind)


def load(reffile=REFFILE):
    """ Optima's inputs (simpars at dt=0.2) and its calibrated run """
    return sc.objdict(sc.load(reffile))


class Regain:
    """ Optima's chance, in a step of 0.2 years, that someone on unsuppressive ART regains suppression: viral load tests (numvlmon) found failing and switched (regainvs), spread over everyone on ART, while the suppressed proportion is not held """
    def __init__(self, ref):
        sp, st = ref.simpars, ref.settings
        self.tvec, self.eps = sp['tvec'], 1e-3
        names = [s.replace('-', '_') for s in st['statelabels']]
        self.alltx = [n for n in names if n.startswith(tuple(TREATED))]
        self.series = np.where(self.tvec <= sp['fixpropsupp'], sp['regainvs']*sp['numvlmon']*REF_DT, 0.0)
        self._ti, self._prob = None, None

    def prob(self, sim):
        if self._ti != sim.ti:
            eng = sim.diseases.optima.engine
            ontx = sum(eng.y[eng.state_comps[n]].sum() for n in self.alltx)
            self._ti, self._prob = sim.ti, min(Transmission.at(self, self.series, sim.now.years)/(self.eps + ontx), 1)
        return self._prob


class Births:
    """
    Optima's births and mother-to-child transmission (model.py, deprecated_births(): projects before 2.12.0, which do
    not diagnose mothers for PMTCT), worked out once a step from the compartments at its start. Returns births per year
    into each population: uninfected, infected to undiagnosed mothers, and infected to diagnosed mothers.
    """
    def __init__(self, ref, mtct=True):
        sp, st = ref.simpars, ref.settings
        self.tvec, self.eps = sp['tvec'], 1e-3
        self.names = [s.replace('-', '_') for s in st['statelabels']]
        rows = lambda cs: [self.names.index(f'{c}_{k}') for c in cs for k in CD4]
        self.undx, self.alldx, self.alltx = rows(['undx']), rows(['dx', 'care', 'lost', 'usvl', 'svl']), rows(TREATED)
        assert np.all(np.asarray(sp['relhivbirth']) == 1) and np.all(np.isnan(sp['proppmtct'])), 'relhivbirth other than 1, or proppmtct, are not handled yet'
        self.rates = np.einsum('ij,ik->ijk', np.asarray(sp['birthtransit']), sp['birth']) # Mother pop × child pop × time, births per woman per year
        effmtct = sp['mtctbreast']*sp['breast'] + sp['mtctnobreast']*(1 - sp['breast'])
        self.effmtct = effmtct if mtct else 0*effmtct
        self.pmtcteff = (1 - sp['effpmtct'])*self.effmtct
        self.numpmtct = sp['numpmtct']
        self.transmission = None # To read the compartments
        self._ti, self._out = None, None

    def at(self, arr, year):
        return Transmission.at(self, arr, year)

    def inflows(self, sim):
        if self._ti == sim.ti:
            return self._out
        eng = sim.diseases.optima.engine
        people = np.array([eng.y[eng.state_comps[n]].sum(axis=0) for n in self.names])
        year = sim.now.years
        rates = self.at(self.rates, year) # Per year
        total, undx, alldx, alltx = people.sum(axis=0), people[self.undx].sum(axis=0), people[self.alldx].sum(axis=0), people[self.alltx].sum(axis=0)
        pregdx = (alldx*rates.sum(axis=1)).sum()*REF_DT # Diagnosed women giving birth in a step
        proppmtct = min(self.at(self.numpmtct, year)*REF_DT/(self.eps*REF_DT + pregdx), 1) # Optima: numpmtct over the steps of a year, over those eligible
        effmtct, pmtcteff = self.at(self.effmtct, year), self.at(self.pmtcteff, year)
        births = rates.T @ total
        mtctundx = rates.T @ undx*effmtct
        mtctdx = rates.T @ alldx*(proppmtct*pmtcteff + (1 - proppmtct)*effmtct) + rates.T @ alltx*pmtcteff # Mothers on ART count as both on treatment and eligible for PMTCT, as in Optima
        self._ti, self._out = sim.ti, (births - mtctundx - mtctdx, mtctundx, mtctdx)
        return self._out


def cascade_probs(sp):
    """
    Optima's probabilities, in a step of 0.2 years, of moving along the cascade for those not on ART (model.py, "Calculate
    probabilities of shifting along cascade"): for each untreated cascade state, an array of CD4 stage × population × time
    """
    shape = sp['death'].shape
    full = lambda x: np.broadcast_to(np.asarray(x, dtype=float), shape)
    aids = lambda base, aidsval: np.array([base if i < AIDS else aidsval for i in range(len(CD4))])
    test, aidstest = full(sp['hivtest']*REF_DT), full(sp['aidstest']*REF_DT)
    link, aidslink = 1 - np.exp(-REF_DT/full(np.maximum(1e-3, sp['linktocare']))), 1 - np.exp(-REF_DT/full(np.maximum(1e-3, sp['aidslinktocare'])))
    ret, loss, aidsloss = full(sp['returntocare']*REF_DT), full(sp['leavecare']*REF_DT), full(sp['aidsleavecare']*REF_DT)
    return dict(
        undx = aids(test, np.maximum(aidstest, test)),
        dx   = aids(link, np.maximum(aidslink, link)),
        care = aids(loss, np.minimum(aidsloss, loss)),
        lost = aids(ret, 1 - (1 - aidstest)*(1 - ret)),
    )


class Inputs:
    """ Optima's input series, as Starsim values by population group and year """
    def __init__(self, ref):
        self.sp = ref.simpars
        self.tvec = self.sp['tvec']
        self.pops = list(self.sp['popkeys'])
        self.npops = len(self.pops)

    def by(self, arr):
        """ A (pops × times) array, or one value per pop, as ss.by(pop=...) """
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 1:
            return ss.by(pop=dict(zip(self.pops, arr)))
        return ss.by(pop={p: ss.timeseries(self.tvec, row, mode='linear') for p, row in zip(self.pops, arr)})

    def series(self, row):
        """ One time series """
        return ss.timeseries(self.tvec, np.asarray(row, dtype=float), mode='linear')


class Transmission:
    """
    Optima's force of infection (model.py, "Calculate probability of getting infected"), worked out once a step from the
    compartments at its start, as Optima does. Returns the rate of infection (per year) for each population, for the
    uncircumcised and circumcised susceptibles: Optima's probability of infection in a step of 0.2 years, divided by 0.2.
    """
    def __init__(self, ref):
        sp, st = ref.simpars, ref.settings
        self.tvec = sp['tvec']
        pops = list(sp['popkeys'])
        self.names = [s.replace('-', '_') for s in st['statelabels']]
        self.eps = 1e-3 # Optima's settings.eps
        male, female = np.asarray(sp['male']), np.asarray(sp['female'])

        # Infectiousness of each state (Optima's alltrans)
        cd4trans = np.array([sp[f'cd4trans{k}'] for k in CD4])*sp['transnorm']
        dxfactor = 1 - sp['effdx']
        alltrans = np.zeros(len(self.names))
        for c, factor in dict(undx=1, dx=dxfactor, care=dxfactor, lost=dxfactor, usvl=dxfactor*(1 - sp['efftxunsupp'])*dxfactor, svl=dxfactor*(1 - sp['efftxsupp'])*dxfactor).items():
            alltrans[[self.names.index(f'{c}_{k}') for k in CD4]] = cd4trans*factor
        self.alltrans = alltrans
        self.hiv = alltrans > 0

        # Effects on susceptibility by population, time, and circumcision status (Optima's alleff)
        prepeff = 1 - sp['effprep']*sp['prep'] - sp['effpep']*np.minimum(sp['pep'], 1 - sp['prep'])
        stieff = 1 + sp['effsti']*sp['stiprev']
        circeff = 1 - sp['propcirc']*sp['effcirc']
        allcirceff = np.outer([1, 1 - sp['effcirc']], male) + np.outer([1, 1], female) # Circumcision status × population
        self.alleff = np.einsum('ab,ab,ab,ca->acb', prepeff, stieff, circeff, allcirceff) # Population × circumcision status × time

        # Partnerships: acts per year, per-act probability, condom effect (Optima's version < 2.12.0 behaviour)
        p1, p2, acts, trans, cond = [], [], [], [], []
        for act in ['reg', 'cas', 'com']:
            for key in sorted(sp[f'acts{act}']):
                i, j = pops.index(key[0]), pops.index(key[1])
                total = np.asarray(sp[f'acts{act}'][key], dtype=float)*np.ones(len(self.tvec))
                if   male[i] and male[j]:   t = (sp['transmmi'] + sp['transmmr'])/2 # Insertive and receptive acts half each
                elif male[i] and female[j]: t = sp['transmfi']
                elif female[i] and male[j]: t = sp['transmfr']
                else: raise ValueError(f'Cannot tell the sexes of {key}')
                c = sp[f'cond{act}'].get(key, sp[f'cond{act}'].get((key[1], key[0]), 0.0))
                p1.append(i); p2.append(j); acts.append(total); trans.append(t*np.ones(len(self.tvec))); cond.append(1 - np.asarray(c)*sp['effcondom']*np.ones(len(self.tvec)))
        self.p1, self.p2 = np.array(p1), np.array(p2)
        self.acts, self.trans, self.cond = np.array(acts), np.array(trans), np.array(cond)
        self.force, self.inhomopar = np.asarray(sp['force'], dtype=float), np.asarray(sp['inhomo'], dtype=float)
        self.bg = sp['death']*REF_DT + sp['propemigrate']*REF_DT
        self.npops = len(pops)
        self._ti, self._rate = None, None

    def at(self, arr, year):
        """ Linear interpolation of a series over Optima's dt=0.2 time points (last axis) """
        i = int(np.clip(np.searchsorted(self.tvec, year, side='right') - 1, 0, len(self.tvec) - 2))
        w = np.clip((year - self.tvec[i])/(self.tvec[i+1] - self.tvec[i]), 0, 1)
        return arr[..., i]*(1 - w) + arr[..., i+1]*w

    def people(self, sim):
        eng = sim.diseases.optima.engine
        return np.array([eng.y[eng.state_comps[n]].sum(axis=0) for n in self.names]) # States × populations

    def rates(self, sim):
        if self._ti == sim.ti:
            return self._rate
        year = sim.now.years
        people = self.people(sim)
        allpeople = people.sum(axis=0)
        effprev = self.alltrans[:, None]*people/allpeople # States × populations
        prev = people[2:].sum(axis=0)/allpeople # Everyone with HIV (Optima's allplhiv), including those suppressed
        h = self.inhomopar
        inhomo = (h + self.eps)/(np.exp(h + self.eps) - 1)*np.exp(h*(1 - prev))

        acts = self.at(self.acts, year)*REF_DT
        whole, frac = np.floor(acts), acts - np.floor(acts)
        alleff = self.at(self.alleff, year)[self.p1] # Partnerships × circumcision status
        q = (self.at(self.trans, year)*self.at(self.cond, year))[:, None, None]*alleff[:, :, None]*effprev[self.hiv][:, self.p2].T[:, None, :] # Partnerships × circ × infected states
        noinf = (1 - np.minimum(frac[:, None, None]*q, 1))*(1 - np.minimum(q*(whole[:, None, None] != 0), 1))**whole[:, None, None]
        lognoinf = np.zeros((self.npops, self.npops) + q.shape[1:]) # Acquiring pop × partner pop × circ × infected state
        np.add.at(lognoinf, (self.p1, self.p2), np.log(noinf))
        bg = self.at(self.bg, year)
        prob = (1 - np.exp(lognoinf)).sum(axis=(1, 3)).T*self.force*inhomo*(1 - bg) # Circ × acquiring pop
        prob = np.minimum(prob, 1 - self.eps - bg.max())
        self._ti, self._rate = sim.ti, prob/REF_DT
        return self._rate


def make_disease(ref, components=('demography', 'natural_history')):
    """ The Optima′ disease: Optima's 38 states, stratified by its 39 populations """
    inp = Inputs(ref)
    sp, pops = inp.sp, inp.pops
    names = [s.replace('-', '_') for s in ref.settings['statelabels']]
    st = sc.objdict({s.name: s for s in ss.states(' '.join(names))})
    hiv = [st[f'{c}_{k}'] for c in UNTREATED + TREATED for k in CD4]
    trs = sc.objdict()

    # Background death: Optima removes it from everyone (1-bg), after HIV death (1-d)
    bg = sp['death']*REF_DT # Probability per step, by pop and time
    d = np.array([sp[f'death{k}'] for k in CD4])*REF_DT # HIV death probability per step by CD4 (relhivdeath is 1 for every pop here)
    assert np.all(np.asarray(sp['hivdeath']) == 1), 'Relative HIV death rates other than 1 are not handled yet'
    tau = np.array([sp['progacute'], sp['proggt500'], sp['proggt350'], sp['proggt200'], sp['proggt50']])
    prog = 1 - np.exp(-REF_DT/tau) # Probability per step of moving to the next CD4 stage

    if 'transmission' in components:
        foi = Transmission(ref)
        for i, s in enumerate([st.susreg, st.progcirc]):
            trs[f'infection_{s.name}'] = s >> st.undx_acute @ ss.by(pop=lambda sim, i=i: ss.peryear(foi.rates(sim)[i]))

    if 'demography' in components:
        for s in st.values(): # Uninfected and on ART: background death only (ART deaths come with treatment)
            if s.name in ('susreg', 'progcirc'):
                trs[f'bgdeath_{s.name}'] = ss.death(s, rate=ss.peryear(1)*inp.by(bg/REF_DT))

    if 'natural_history' in components:
        q = cascade_probs(sp) if 'cascade' in components else {c: np.zeros((len(CD4),) + bg.shape) for c in UNTREATED}
        for c in UNTREATED:
            for i, k in enumerate(CD4):
                s = st[f'{c}_{k}']
                trs[f'hivdeath_{s.name}'] = ss.death(s, rate=ss.peryear(d[i]/REF_DT))
                trs[f'bgdeath_{s.name}'] = ss.death(s, rate=ss.peryear(1)*inp.by(bg*(1 - d[i])/REF_DT))
                survive = (1 - d[i])*(1 - bg) # Optima applies death first, and every other move to those who survive it
                p = prog[i] if i < len(CD4) - 1 else 0 # No progression from CD4<50
                nxt = CD4[min(i+1, len(CD4) - 1)]
                moves = {}
                for dst, prob in [(f'{c}_{nxt}', p*(1 - q[c][i])), (f'{NEXT[c]}_{k}', (1 - p)*q[c][i]), (f'{NEXT[c]}_{nxt}', p*q[c][i])]:
                    moves[dst] = moves.get(dst, 0) + prob # From CD4<50, the last two are the same move
                for dst, prob in moves.items():
                    if dst != s.name and np.any(prob):
                        trs[f'{s.name}_to_{dst}'] = s >> st[dst] @ (ss.peryear(1)*inp.by(prob*survive/REF_DT))

    if 'treatment' in components:
        q = cascade_probs(sp)
        fixsupp = sp['fixpropsupp'] # Treatment failure runs only while the suppressed proportion is not held (Optima's userate(propsupp))
        fail = np.where(inp.tvec <= fixsupp, sp['treatfail']*REF_DT, 0.0)*np.ones(bg.shape)
        regain = Regain(ref) # Regaining suppression through viral load monitoring: a probability per step from the number on ART, so worked out each step
        svlrecov = 1 - np.exp(-REF_DT/np.array([np.inf, np.inf, sp['svlrecovgt350'], sp['svlrecovgt200'], sp['svlrecovgt50'], sp['svlrecovlt50']]))
        u = lambda key: sp[key]*REF_DT # Rates on unsuppressive ART, used as probabilities per step
        cd4moves = dict( # Probability of each CD4 stage next step, from each stage (model.py lines 303-351)
            svl  = [{0: 1 - prog[0], 1: prog[0]}, {1: 1}] + [{i: 1 - svlrecov[i], i-1: svlrecov[i]} for i in range(2, 6)],
            usvl = [{0: 1 - prog[0], 1: prog[0]}, {1: 1 - u('usvlproggt500'), 2: u('usvlproggt500')}] +
                   [{i-1: u(f'usvlrecov{CD4[i]}'), i: 1 - u(f'usvlrecov{CD4[i]}') - u(f'usvlprog{CD4[i]}'), i+1: u(f'usvlprog{CD4[i]}')} for i in range(2, 5)] +
                   [{4: u('usvlrecovlt50'), 5: 1 - u('usvlrecovlt50')}],
        )
        for c, other, switch in [('svl', 'usvl', fail), ('usvl', 'svl', regain)]:
            for i, k in enumerate(CD4):
                s = st[f'{c}_{k}']
                dtx = d[i]*sp[f'death{c}'] # HIV death on ART
                trs[f'hivdeath_{s.name}'] = ss.death(s, rate=ss.peryear(dtx/REF_DT))
                trs[f'bgdeath_{s.name}'] = ss.death(s, rate=ss.peryear(1)*inp.by(bg*(1 - dtx)/REF_DT))
                survive = (1 - dtx)*(1 - bg)
                loss = q['care'][i]
                cascade = {c: (1 - loss)*(1 - switch), other: (1 - loss)*switch, 'lost': loss} if c == 'svl' else {'lost': loss}
                for j, pcd4 in cd4moves[c][i].items():
                    for dc in [c, other, 'lost']:
                        dst = f'{dc}_{CD4[j]}'
                        if dst == s.name: continue
                        if c == 'usvl' and dc != 'lost': # The chance of regaining suppression changes each step
                            base = pcd4*(1 - loss)*survive*np.ones(bg.shape)/REF_DT
                            factor = (lambda r: r) if dc == 'svl' else (lambda r: 1 - r)
                            rate = lambda sim, base=base, factor=factor: ss.peryear(Transmission.at(regain, base, sim.now.years)*factor(regain.prob(sim)))
                            trs[f'{s.name}_to_{dst}'] = s >> st[dst] @ ss.by(pop=rate)
                        else:
                            prob = pcd4*cascade[dc]*survive*np.ones(bg.shape)
                            if np.any(prob):
                                trs[f'{s.name}_to_{dst}'] = s >> st[dst] @ (ss.peryear(1)*inp.by(prob/REF_DT))

    # Births into F0-4 and M0-4, from all women in the mother populations, at Optima's birth rate per woman, with mother-to-child transmission
    if 'demography' in components:
        births = Births(ref, mtct='mtct' in components)
        for i, (name, dst) in enumerate([('birth', st.susreg), ('birth_undx', st.undx_acute), ('birth_dx', st.dx_acute)]):
            if i == 0 or 'mtct' in components:
                trs[name] = ss.null >> dst @ ss.by(pop=lambda sim, i=i: ss.peryear(births.inflows(sim)[i]))

    # Ageing between bands (Optima's agerate per year, by source population), and risk transitions (FSW, Clients <-> general)
    if 'demography' in components:
        at, a = np.asarray(sp['agetransit']), sp['agerate']
        onward = {i: j for i, j in zip(*np.nonzero(at))} # The band each band ages into
        for i, j in onward.items():
            rate = at[i, j]*a[i]
            if OPTIONS.age_skip and j in onward: # Optima's old ageing code (projects before 2.12.0) ages people who arrive in a band again in the same step, so some of those leaving skip it
                k = onward[j]
                skip = rate*at[j, k]*a[j]*REF_DT
                trs[f'skip_{pops[i]}'] = ss.transfer(pop=(pops[i], pops[k]), rate=ss.peryear(1)*inp.series(OPTIONS.skip_scale*skip))
                rate = rate - OPTIONS.skip_scale*skip
            trs[f'age_{pops[i]}'] = ss.transfer(pop=(pops[i], pops[j]), rate=ss.peryear(1)*inp.series(OPTIONS.age_scale*rate))
        rt = np.asarray(sp['risktransit'], dtype=float)
        size = sp['popsize']
        for i, j in zip(*np.nonzero(rt)) if OPTIONS.risk else []:
            r = 1 - np.exp(-REF_DT/rt[i, j]) # Probability per step
            trs[f'risk_{pops[i]}_{pops[j]}'] = ss.transfer(pop=(pops[i], pops[j]), rate=ss.peryear(r/REF_DT))
            trs[f'risk_{pops[j]}_{pops[i]}'] = ss.transfer(pop=(pops[j], pops[i]), rate=ss.peryear(1)*inp.series(r/REF_DT*size[i]/size[j])) # Balanced: as many move back (population sizes are forced)

    # Initial state: Optima's own, as the share of each population in each state
    people0 = ref.raw['people'][:, :, 0]
    shares0 = people0/people0.sum(axis=0)
    init = {name: inp.by(shares0[i]) for i, name in enumerate(names) if shares0[i].any()}

    dis = ss.Disease(*st.values(), **trs, init=init, name='optima')
    dis.stratify(pop=dict(zip(pops, people0.sum(axis=0)/people0.sum())))
    return dis, st


class ForcePopsize(ss.Intervention):
    """
    Optima's population forcing, applied to the compartments each step: populations with no inflows from ageing or births
    gain or lose susceptibles to match their size; every population then has everyone not on ART scaled to match it.

    Rescaling compartments has no public Starsim interface yet, so this writes the engine's compartments directly.
    """
    def __init__(self, ref, **kwargs):
        super().__init__(**kwargs)
        sp = ref.simpars
        self.tvec, self.popsize = sp['tvec'], sp['popsize']
        inflows = np.asarray(sp['agetransit']).sum(axis=0) + np.asarray(sp['birthtransit']).sum(axis=0)
        self.noinflows = np.nonzero(inflows == 0)[0]
        self.statelabels = [s.replace('-', '_') for s in ref.settings['statelabels']]

    def step(self):
        eng = self.sim.diseases.optima.engine
        y, comps = eng.y, eng.state_comps
        rows = lambda names: np.concatenate([comps[n] for n in names])
        sus = rows(['susreg', 'progcirc'])
        notart = rows([n for n in self.statelabels if not n.startswith(tuple(TREATED))])
        art = rows([n for n in self.statelabels if n.startswith(tuple(TREATED))])
        year = self.sim.now.years if hasattr(self.sim.now, 'years') else float(self.sim.now)
        target = np.array([np.interp(year, self.tvec, row) for row in self.popsize])
        j = self.noinflows
        allsus = y[np.ix_(sus, j)].sum(axis=0)
        new = target[j] - y[:, j].sum(axis=0)
        y[np.ix_(sus, j)] *= 1 + new/allsus # Add (or remove) susceptibles, in proportion
        wanted = target - y[art].sum(axis=0)
        y[notart] *= wanted/y[notart].sum(axis=0)
        return


class CascadeTargets(ss.Intervention):
    """
    Optima's rules on the cascade's levels, applied each step after the population forcing (model.py, "Proportions"):
    the number on ART from numtx (people in care start it, lowest CD4 first before allcd4eligibletx; Optima's chance of
    viral suppression within a step of 0.2 years, at the start), and proportions held from their fix years.

    Moving a computed number of people has no interface that records it, so this writes the engine's compartments directly.
    """
    def __init__(self, ref, **kwargs):
        super().__init__(**kwargs)
        sp = ref.simpars
        self.tvec, self.numtx = sp['tvec'], sp['numtx']
        self.eligible = sp['allcd4eligibletx']
        self.treatvs = 1 - np.exp(-REF_DT/max(1e-3, sp['treatvs'])) # Share suppressed at the start of ART
        self.fixyears = dict(propdx=sp['fixpropdx'], propcare=sp['fixpropcare'], proptx=sp['fixproptx'], propsupp=sp['fixpropsupp'])
        for name in self.fixyears:
            assert np.all(np.isnan(sp[name])), f'Given values of {name} are not handled yet'
        self.held = {} # Proportions held from their fix years
        self.eps = 1e-3

    def step(self):
        eng = self.sim.diseases.optima.engine
        y, comps = eng.y, eng.state_comps
        idx = lambda c: np.array([comps[f'{c}_{k}'][0] for k in CD4]) # One row per CD4 stage (each state has one compartment)
        undx, dx, care, lost, usvl, svl = [idx(c) for c in UNTREATED + TREATED]
        alltx, allcare = np.concatenate([usvl, svl]), np.concatenate([care, usvl, svl])
        alldx, allplhiv = np.concatenate([dx, care, lost, usvl, svl]), np.concatenate([undx, dx, care, lost, usvl, svl])
        year = self.sim.now.years
        structure = dict(propdx=(undx, dx, alldx, allplhiv), propcare=(np.concatenate([dx, lost]), care, allcare, alldx),
                         proptx=(care, alltx, alltx, allcare), propsupp=(usvl, svl, svl, alltx))
        for name, (lower, upper, numer, denom) in structure.items():
            if name not in self.held and year >= self.fixyears[name] - 1e-9: # Hold the proportion from its fix year, at its value then
                self.held[name] = y[numer].sum()/(self.eps + y[denom].sum())
            actual, available = y[numer].sum(), y[denom].sum()
            if name in self.held:
                wanted = self.held[name]*available
            elif name == 'proptx':
                wanted = np.interp(year, self.tvec, self.numtx)
            else:
                continue
            diff = wanted - actual
            if diff > self.eps: # Move people up the cascade
                up = y[lower]
                total = up.sum()
                if total <= self.eps: continue
                diff = min(diff, total - self.eps)
                if name == 'proptx' and year < self.eligible: # Lowest CD4 first
                    movers = np.zeros_like(up)
                    for i in reversed(range(len(CD4))):
                        if diff > self.eps and up[i].sum() > self.eps:
                            n = min(diff, up[i].sum() - self.eps)
                            movers[i] = n*up[i]/up[i].sum()
                            diff -= movers[i].sum()
                else:
                    movers = diff*up/total
                y[lower] -= movers
                if name == 'proptx':
                    y[usvl] += movers*(1 - self.treatvs)
                    y[svl] += movers*self.treatvs
                elif name == 'propcare':
                    y[care] += movers[:len(CD4)] + movers[len(CD4):] # From dx and from lost
                else:
                    y[upper] += movers
            elif diff < -self.eps: # Move people down the cascade
                down = y[upper] if name != 'proptx' else y[alltx]
                total = down.sum()
                if total <= self.eps: continue
                movers = min(-diff, total - self.eps)*down/total
                if name == 'proptx':
                    y[alltx] -= movers
                    y[care] += movers[:len(CD4)] + movers[len(CD4):]
                else:
                    y[upper] -= movers
                    y[lost if name == 'propcare' else lower] += movers
        return


def make_sim(ref=None, components=('demography', 'natural_history'), dt=0.02, end=2030, **kwargs):
    ref = ref or load()
    dis, st = make_disease(ref, components)
    total = ref.raw['people'][:, :, 0].sum()
    interventions = [ForcePopsize(ref, name='force_popsize')] if 'demography' in components and ref.simpars['forcepopsize'] else []
    if 'treatment' in components:
        interventions.append(CascadeTargets(ref, name='cascade_targets'))
    sim = ss.Sim(dis, interventions=interventions, solver='ode', start=1990, stop=end, dt=dt, total_pop=total, verbose=0, **kwargs)
    return sim


if __name__ == '__main__':
    T = sc.timer()
    sim = make_sim().run()
    T.toc('Optima′ run')
    r = sim.results.optima
    print({k: v for k, v in list(r.items())[:5]})
