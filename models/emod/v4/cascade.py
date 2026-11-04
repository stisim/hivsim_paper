"""
EMOD′ care cascade: Daniel's campaign file, interpreted as EMOD runs it (docs/emod_inventory.md, sections 7-9, 11, 13).

The campaign is a state machine of events. A listener (NodeLevelHealthTriggeredIV) gives an intervention to whoever
broadcasts one of its events, if they qualify (node, sex, age, properties, coverage, from its start year); the
intervention may broadcast further events. This module reads the events from the campaign JSON and runs them with
EMOD's semantics (campaign/Interventions.cpp, InterventionsContainer::Update, interventions/*.cpp):

- An intervention is refused if the person's `InterventionStatus` is one of its `Disqualifying_Properties`, when it
  is given and again whenever it is updated; on its first update it sets its `New_Property_Value`.
- Interventions given during a person's update are updated in the same step, so chains of events with no delay, and
  delays shorter than what is left of the step, resolve within the step. Here they run in waves within the step.
- Decisions (HIVSigmoidByYearAndSexDiagnostic, HIVPiecewiseByYearAndSexDiagnostic, HIVRapidHIVDiagnostic,
  STIIsPostDebut, HIVDrawBlood, HIVARTStagingCD4AgnosticDiagnostic, HIVARTStagingByCD4Diagnostic, HIVRandomChoice)
  broadcast a positive or negative event at once. A sigmoid reads only its `Ramp_*` keys, so a `Time_Value_Map` given
  to one is ignored (the ANC uptake), as in EMOD.
- Delays (HIVDelayedIntervention, HIVMuxer) count down by the step from the step they are given in; a muxer is refused
  if the person already holds `Max_Entries` of that name.
- PropertyValueChanger, BroadcastEvent, ARTMortalityTable and ARTDropout (through the HIV module), PMTCT (its efficacy
  for 280 days), MaleCircumcision and ControlledVaccine (PrEP) act at once.
- NChooserEventCoordinatorHIV (VMMC, PrEP) gives a number of interventions each year (see `NChooser`).

Events from elsewhere: STIDebut (risk module), NewlySymptomatic (HIV module), SixWeeksOld (age), TwelveWeeksPregnant
(demography); the last two happen after EMOD updates interventions, so what they give is first updated next step.
The parts of the campaign other modules already run (sex work, STI co-infection, traditional circumcision, seeding)
and the Accessibility tracker (everyone is Accessibility:Yes) are skipped. People: `people.istatus`
(InterventionStatus), `people.artstatus` (ARTStatus), `people.onprep`.
"""
import json
import numpy as np
import starsim as ss

SKIP = ['Commercial_', 'ModifyStiCoInfectionStatus', 'TraditionalMaleCircumcision', 'OutbreakIndividual', 'ReferenceTrackingEventCoordinator']
DELAYS = ('HIVDelayedIntervention', 'HIVMuxer')


def sigmoid(year, c):
    lo, hi, mid, rate = c.get('Ramp_Min', 0), c.get('Ramp_Max', 1), c.get('Ramp_MidYear', 2000), c.get('Ramp_Rate', 1)
    return lo + (hi - lo)/(1 + np.exp(-rate*(year - mid)))


def piecewise(tvm, year, default, order=0):
    """ InterpolatedValueMap: getValuePiecewiseConstant (order 0) or getValueLinearInterpolation """
    t, v = np.array(tvm['Times'], float), np.array(tvm['Values'], float)
    if year < t[0]:
        return default
    return v[np.searchsorted(t, year, side='right') - 1] if order == 0 else float(np.interp(year, t, v))


class EmodCascade(ss.Intervention):
    """
    Args:
        inp (objdict): EMOD's inputs for one parameter set, from `inputs.load()`
    """
    def __init__(self, inp, name='emodcascade', **kwargs):
        super().__init__(name=name, **kwargs)
        self.inp = inp
        self.statuses = inp.ip_values['InterventionStatus']
        self.artstatuses = inp.ip_values['ARTStatus']
        self.listeners, self.oneoffs, self.choosers, self.templates = [], [], [], []
        for i, ev in enumerate(inp.events):
            if any(k in json.dumps(ev) for k in SKIP):
                continue
            ec = ev['Event_Coordinator_Config']
            start = ev.get('Start_Year', inp.start + ev.get('Start_Day', 0)/365)
            nodes = ev['Nodeset_Config'].get('Node_List', [1, 2, 3])
            if ec['class'] == 'NChooserEventCoordinatorHIV':
                self.choosers.append(dict(ec=ec, nodes=nodes, tid=self.add(ec['Intervention_Config'])))
                continue
            ic = ec['Intervention_Config']
            if ic['class'] == 'NodeLevelHealthTriggeredIV':
                self.listeners.append(dict(triggers=ic['Trigger_Condition_List'], start=start, nodes=nodes, target=ic, tid=self.add(ic['Actual_IndividualIntervention_Config'])))
            else:
                self.oneoffs.append(dict(start=start, nodes=nodes, target=ec, tid=self.add(ic), done=False))
        self.by_trigger = {}
        for k, li in enumerate(self.listeners):
            for t in li['triggers']:
                self.by_trigger.setdefault(t, []).append(k)
        self.define_states(ss.BoolArr('vmmc_given', default=False)) # Holds the programme's MaleCircumcision (named Any_MC)
        return

    def add(self, cfg):
        cfg = dict(cfg)
        cfg['disq'] = [self.statuses.index(d.split(':')[1]) for d in cfg.get('Disqualifying_Properties', []) if d.startswith('InterventionStatus:')]
        npv = cfg.get('New_Property_Value', '')
        cfg['npv'] = self.statuses.index(npv.split(':')[1]) if npv.startswith('InterventionStatus:') else None
        self.templates.append(cfg)
        return len(self.templates) - 1

    def init_pre(self, sim):
        super().init_pre(sim)
        self.istatus = sim.people.add_state(ss.IntArr('istatus', default=0))
        self.artstatus = sim.people.add_state(ss.IntArr('artstatus', default=0))
        self.onprep = sim.people.add_state(ss.BoolArr('onprep', default=False))
        self.rng = np.random.default_rng(sim.pars.rand_seed + 6)
        self.items = dict(uid=np.zeros(0, np.int64), tid=np.zeros(0, np.int64), left=np.zeros(0)) # Interventions waiting (delays, PMTCT)
        self.muxers = {}
        self.carry = []
        self.choosers_state = [NChooser(c, self) for c in self.choosers]
        return

    def step(self):
        """ Nothing here: the HIV module calls `update()` at EMOD's point in the step """
        return

    # %% Who qualifies

    def qualifies(self, target, uids, nodes=None):
        ppl = self.sim.people
        ok = ppl.alive.raw[uids] & ~self.dying[uids]
        if nodes is not None:
            ok &= np.isin(self.sim.people.node.raw[uids] + 1, nodes)
        td = target.get('Target_Demographic', 'Everyone')
        if td in ('ExplicitGender', 'ExplicitAgeRangesAndGender') and target.get('Target_Gender', 'All') != 'All':
            ok &= ppl.female.raw[uids] == (target['Target_Gender'] == 'Female')
        if td == 'ExplicitAgeRangesAndGender':
            age = ppl.emod_age.raw[uids]/365
            ok &= (age >= target.get('Target_Age_Min', 0)) & (age < target.get('Target_Age_Max', 125))
        ok &= self.has_props(uids, target.get('Property_Restrictions_Within_Node'), target.get('Property_Restrictions'))
        cov = target.get('Demographic_Coverage', 1)
        if cov < 1:
            ok &= self.rng.random(len(uids)) < cov
        return uids[ok]

    def has_props(self, uids, within=None, plain=None):
        """ Property restrictions: within-node ones are alternatives (any), plain ones all hold """
        def holds(kv):
            key, val = kv
            if key == 'Accessibility':
                return np.full(len(uids), val == 'Yes')
            if key == 'Risk':
                return self.sim.people.risk.raw[uids] == ['LOW', 'MEDIUM', 'HIGH'].index(val)
            if key == 'InterventionStatus':
                return self.istatus.raw[uids] == self.statuses.index(val)
            raise NotImplementedError(f'property restriction {key}')
        ok = np.ones(len(uids), bool)
        if within:
            ok &= np.any([np.all([holds(kv) for kv in d.items()], axis=0) for d in within], axis=0)
        for r in plain or []:
            ok &= holds(r.split(':'))
        return ok

    def allowed(self, tid, uids):
        disq = self.templates[tid]['disq']
        return uids[~np.isin(self.istatus.raw[uids], disq)] if disq else uids

    # %% Events and interventions

    def broadcast(self, name, uids, year):
        """ Give the interventions of every listener for this event; returns (template, people) to update """
        out = []
        if not len(uids) or name not in self.by_trigger:
            return out
        for k in self.by_trigger[name]:
            li = self.listeners[k]
            if year < li['start']:
                continue
            got = self.allowed(li['tid'], self.qualifies(li['target'], uids, li['nodes']))
            if len(got):
                out.extend(self.give(li['tid'], got))
        return out

    def give(self, tid, uids):
        """ Distribution: what acts at once (ART, dropout, PMTCT, circumcision, PrEP) does so now; all are then updated """
        c, hiv = self.templates[tid], self.sim.diseases.hiv
        cls = c['class']
        if cls == 'ARTMortalityTable':
            hiv.go_on_art(uids, c, c['ART_Multiplier_On_Transmission_Prob_Per_Act'], c['Days_To_Achieve_Viral_Suppression'])
        elif cls == 'ARTDropout':
            hiv.go_off_art(uids)
        elif cls == 'PMTCT':
            hiv.pmtct[uids] = np.maximum(hiv.pmtct.raw[uids], c['Efficacy']) # ApplyProbMaternalTransmissionModifier keeps the larger
        elif cls == 'MaleCircumcision':
            men = uids[~self.sim.people.female.raw[uids] & ~self.sim.people.circumcised.raw[uids]] # Refused to those already circumcised
            self.sim.people.circumcised[men] = True
            self.vmmc_given[men] = True
        elif cls == 'ControlledVaccine':
            w = c['Waning_Config']
            hiv.prep_eff[uids] = c.get('Initial_Effect', w.get('Initial_Effect', 1))
            hiv.prep_until[uids] = self.t.now('year') + w['Durability_Map']['Times'][-2]/365
        return [(tid, uids)]

    def set_status(self, tid, uids):
        npv = self.templates[tid]['npv']
        if npv is not None:
            self.istatus[uids] = npv

    def run(self, tid, uids, year):
        """ First update of an intervention: returns what it gives in turn """
        c, ppl, hiv = self.templates[tid], self.sim.people, self.sim.diseases.hiv
        cls = c['class']
        if cls not in ('ARTMortalityTable', 'ARTDropout', 'PropertyValueChanger', 'BroadcastEvent', 'PMTCT', 'MaleCircumcision', 'ControlledVaccine') or c['disq']:
            uids = self.allowed(tid, uids)
        if not len(uids):
            return []
        self.set_status(tid, uids)
        if cls in DELAYS:
            return self.start_delay(tid, uids, year)
        if cls == 'PMTCT':
            self.store(tid, uids, np.full(len(uids), 280.0 - self.dt_days))
            return []
        if cls == 'PropertyValueChanger':
            key, val = c['Target_Property_Key'], c['Target_Property_Value']
            if key == 'InterventionStatus':
                self.istatus[uids] = self.statuses.index(val)
            elif key == 'ARTStatus':
                self.artstatus[uids] = self.artstatuses.index(val)
            elif key == 'OnPrEP':
                self.onprep[uids] = val == 'True'
            return []
        if cls == 'BroadcastEvent':
            return self.broadcast(c['Broadcast_Event'], uids, year)
        if cls == 'HIVRandomChoice':
            u = self.rng.random(len(uids))
            cum = np.cumsum(c['Choice_Probabilities'])
            pick = np.searchsorted(cum, u, side='left') # First choice whose cumulative probability is at least u
            out = []
            for k, name in enumerate(c['Choice_Names']):
                out += self.broadcast(name, uids[pick == k], year)
            return out
        if cls in ('ARTMortalityTable', 'ARTDropout', 'MaleCircumcision', 'ControlledVaccine'):
            return []
        pos = self.decide(c, uids, year)
        out = self.broadcast(c['Positive_Diagnosis_Event'], uids[pos], year) if c.get('Positive_Diagnosis_Event') else []
        if c.get('Negative_Diagnosis_Event'):
            out += self.broadcast(c['Negative_Diagnosis_Event'], uids[~pos], year)
        return out

    def decide(self, c, uids, year):
        ppl, hiv = self.sim.people, self.sim.diseases.hiv
        cls = c['class']
        female = ppl.female.raw[uids]
        if cls == 'HIVSigmoidByYearAndSexDiagnostic':
            p = max(0.0, sigmoid(year, c))*np.where(female, c.get('Female_Multiplier', 1), 1)
            return self.rng.random(len(uids)) < p
        if cls == 'HIVPiecewiseByYearAndSexDiagnostic':
            p = piecewise(c['Time_Value_Map'], year, c.get('Default_Value', 0), c.get('Interpolation_Order', 0))*np.where(female, c.get('Female_Multiplier', 1), 1)
            return self.rng.random(len(uids)) < p
        if cls == 'HIVRapidHIVDiagnostic':
            pos = hiv.infected.raw[uids].copy()
            hiv.ever_tested[uids] = True
            hiv.last_test[uids] = year
            hiv.diagnosed[uids[pos]] = True
            return pos
        if cls == 'STIIsPostDebut':
            return ppl.emod_age.raw[uids]/365 >= self.sim.networks.emodnet.debut.raw[uids]
        if cls == 'HIVDrawBlood':
            hiv.lowest_cd4[uids] = np.minimum(hiv.lowest_cd4.raw[uids], np.where(hiv.infected.raw[uids], hiv.cd4(uids), np.inf))
            return np.ones(len(uids), bool)
        inf = hiv.infected.raw[uids]
        pregnant = ~np.isnan(self.sim.demographics.emoddemog.preg_days.raw[uids])
        if cls == 'HIVARTStagingCD4AgnosticDiagnostic':
            who = np.where(inf, hiv.who_stage(uids), 1.0)
            age = ppl.emod_age.raw[uids]
            adult = age >= c['Adult_Treatment_Age']*365
            res_a = (piecewise(c['Adult_By_WHO_Stage'], year, 6) <= who) | (pregnant & (piecewise(c['Adult_By_Pregnant'], year, 0) > 0))
            res_c = (age <= piecewise(c['Child_Treat_Under_Age_In_Years_Threshold'], year, -1)*365) | (piecewise(c['Child_By_WHO_Stage'], year, 6) <= who)
            return inf & np.where(adult, res_a, res_c)
        if cls == 'HIVARTStagingByCD4Diagnostic':
            cd4 = hiv.lowest_cd4.raw[uids]
            return inf & ((cd4 <= piecewise(c['Threshold'], year, -np.inf)) | (pregnant & (cd4 <= piecewise(c['If_Pregnant'], year, -np.inf))))
        raise NotImplementedError(cls)

    # %% Delays

    def draw_delay(self, c, n):
        d = c.get('Delay_Period_Distribution', 'CONSTANT_DISTRIBUTION')
        u = self.rng.random(n)
        if d == 'CONSTANT_DISTRIBUTION':
            return np.full(n, float(c['Delay_Period_Constant']))
        if d == 'UNIFORM_DISTRIBUTION':
            return c['Delay_Period_Min'] + u*(c['Delay_Period_Max'] - c['Delay_Period_Min'])
        if d == 'EXPONENTIAL_DISTRIBUTION':
            return -np.log(1 - u)*c['Delay_Period_Exponential']
        if d == 'WEIBULL_DISTRIBUTION':
            return c['Delay_Period_Lambda']*(-np.log(1 - u))**(1/c['Delay_Period_Kappa'])
        raise NotImplementedError(d)

    def muxer(self, name):
        n = len(self.sim.people.alive.raw)
        arr = self.muxers.get(name)
        if arr is None or len(arr) < n:
            new = np.zeros(n, np.int64)
            if arr is not None:
                new[:len(arr)] = arr
            self.muxers[name] = arr = new
        return arr

    def start_delay(self, tid, uids, year):
        c = self.templates[tid]
        if c['class'] == 'HIVMuxer':
            m = self.muxer(c['Muxer_Name'])
            np.add.at(m, uids, 1)
            over = m[uids] > c.get('Max_Entries', 1)
            np.add.at(m, uids[over], -1)
            uids = uids[~over]
        left = self.draw_delay(c, len(uids)) - self.dt_days
        due = left <= 0
        self.store(tid, uids[~due], left[~due])
        return self.fire(tid, uids[due], year)

    def store(self, tid, uids, left):
        it = self.items
        it['uid'] = np.concatenate([it['uid'], uids])
        it['tid'] = np.concatenate([it['tid'], np.full(len(uids), tid)])
        it['left'] = np.concatenate([it['left'], left])

    def fire(self, tid, uids, year):
        """ A delay or PMTCT that has run out """
        c = self.templates[tid]
        if c['class'] == 'HIVMuxer':
            np.add.at(self.muxer(c['Muxer_Name']), uids, -1)
        if c['class'] == 'PMTCT':
            self.sim.diseases.hiv.pmtct[uids] = 0.0
            return []
        return self.broadcast(c['Broadcast_Event'], uids, year) if c.get('Broadcast_Event') else []

    def tick(self, year):
        """ Waiting interventions given in earlier steps: disqualified ones end; the others count down """
        it, ppl = self.items, self.sim.people
        if not len(it['uid']):
            return []
        alive = ppl.alive.raw[it['uid']] & ~self.dying[it['uid']]
        disq = np.zeros(len(it['uid']), bool)
        for tid in np.unique(it['tid']):
            m = it['tid'] == tid
            d = self.templates[tid]['disq']
            if d:
                disq[m] = np.isin(self.istatus.raw[it['uid'][m]], d)
        it['left'] = it['left'] - self.dt_days
        end = ~alive | disq
        due = ~end & (it['left'] <= 0)
        out = []
        for tid in np.unique(it['tid'][end | due]):
            m = it['tid'] == tid
            c = self.templates[tid]
            if c['class'] == 'HIVMuxer':
                np.add.at(self.muxer(c['Muxer_Name']), it['uid'][m & end], -1)
            out += self.fire(tid, it['uid'][m & due], year)
        keep = ~(end | due)
        for k in it:
            it[k] = it[k][keep]
        return out

    # %% The step

    def update(self, dying):
        """ One step of the campaign, called by the HIV module after progression and before exposure """
        ppl = self.sim.people
        self.dt_days = self.inp.dt
        year = self.t.now('year')
        self.dying = np.zeros(len(ppl.alive.raw), bool)
        self.dying[dying] = True
        pend, self.carry = self.carry, []
        pend = [(tid, uids[ppl.alive.raw[uids] & ~self.dying[uids]]) for tid, uids in pend]
        for o in self.oneoffs:
            if not o['done'] and year >= o['start']:
                o['done'] = True
                got = self.allowed(o['tid'], self.qualifies(o['target'], ppl.auids, o['nodes']))
                pend += self.give(o['tid'], got)
        for ch in self.choosers_state:
            pend += ch.step(year)
        risk = self.sim.interventions.emodrisk
        pend += self.broadcast('STIDebut', risk.new_debuts, year)
        pend += self.broadcast('NewlySymptomatic', self.sim.diseases.hiv.newly_symptomatic, year)
        pend += self.tick(year)
        while pend:
            nxt = []
            for tid, uids in pend:
                nxt += self.run(tid, uids, year)
            pend = nxt
        # Events after EMOD's update of interventions: what they give is first updated next step
        age = ppl.emod_age.raw[ppl.auids]
        six = ppl.auids[(age - self.dt_days < 42) & (age >= 42)]
        self.carry = self.broadcast('SixWeeksOld', six, year) + self.broadcast('TwelveWeeksPregnant', self.sim.demographics.emoddemog.twelve_weeks, year)
        return


class NChooser:
    """
    NChooserEventCoordinatorHIV: in each period of `Distributions`, a target number of people per node, sex and age
    range, chosen at random among those who qualify (HIV status, not having the intervention, properties), spread
    evenly over the period's steps, with the remainder in the earliest steps. Targets are numbers of people, scaled by
    x_Base_Population and truncated (NChooserEventCoordinator.cpp).
    """
    def __init__(self, cfg, cascade):
        self.cfg, self.cas = cfg, cascade
        self.scale = 0.01
        self.plan = None
        return

    def step(self, year):
        ec, cas = self.cfg['ec'], self.cas
        out = []
        for d in ec['Distributions']:
            if not (d['Start_Year'] <= year < d['End_Year']):
                continue
            key = id(d)
            if self.plan is None or self.plan[0] != key:
                steps = int(round((d['End_Year'] - d['Start_Year'])*12))
                targets = []
                for k, ar in enumerate(d['Age_Ranges_Years']):
                    for sex, col in ((0, 'Num_Targeted_Males'), (1, 'Num_Targeted_Females')):
                        n = int(d[col][k]*self.scale)
                        per = np.full(steps, n//steps) + (np.arange(steps) < n % steps)
                        targets.append((ar, sex, per))
                self.plan = (key, targets, 0)
            key, targets, k = self.plan
            for ar, sex, per in targets:
                n = per[k] if k < len(per) else 0
                if n:
                    out += self.choose(d, ar, sex, n)
            self.plan = (key, targets, k + 1)
        return out

    def choose(self, d, ar, sex, n):
        cas, ppl, hiv = self.cas, self.cas.sim.people, self.cas.sim.diseases.hiv
        u = ppl.auids
        u = u[np.isin(ppl.node.raw[u] + 1, self.cfg['nodes']) & (ppl.female.raw[u] == sex) & ppl.alive.raw[u] & ~cas.dying[u]]
        age = ppl.emod_age.raw[u]/365
        u = u[(age >= ar['Min']) & (age < ar['Max'])]
        state = d['Target_Disease_State'][0]
        if 'HIV_Negative' in state:
            u = u[~hiv.infected.raw[u]]
        if 'Not_Have_Intervention' in state:
            name = d.get('Target_Disease_State_Has_Intervention_Name')
            has = {'Any_MC': cas.vmmc_given, 'PrEP': np.zeros(len(ppl.alive.raw), bool)}.get(name)
            if has is not None:
                u = u[~has[u]]
        u = u[cas.has_props(u, d.get('Property_Restrictions_Within_Node') or None)]
        pick = u[cas.rng.permutation(len(u))[:n]] # Those refused (circumcised traditionally) count towards the target
        return cas.give(self.cfg['tid'], pick)
