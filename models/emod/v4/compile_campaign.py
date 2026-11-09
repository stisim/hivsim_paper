"""
EMOD′ care cascade compiled into Starsim deliveries: the events of Daniel's campaign file (the parts cascade.py runs)
become `ss.Delivery` objects, with EMOD's semantics expressed natively where Starsim has them (needs Starsim's
deliveries with cycles through delays, 4.0.0rc3 or later).

How each part of the campaign maps:
- A listener (NodeLevelHealthTriggeredIV) → deliveries triggered by the event (`when=` the deliveries that broadcast it,
  joined with `|`), with its node, sex and age as conditions on people, `& ~` its Disqualifying_Properties as
  conditions on `people.istatus`, and not dying this step.
- HIVSigmoidByYearAndSexDiagnostic, HIVPiecewiseByYearAndSexDiagnostic → `coverage=` a function of the year (the yes
  branch, which broadcasts the positive event), and the complement as a delivery on `trigger & ~yes` (the negative event).
- HIVRandomChoice → a chain of deliveries, each with the conditional probability as `coverage=`, on the trigger and not
  the ones before.
- New_Property_Value and PropertyValueChanger → `SetIP`, a product that sets `istatus`, `artstatus` or `onprep`.
- BroadcastEvent → a delivery with no product but its status, whose event is the broadcast.
- HIVMuxer → two deliveries: on the trigger, `MuxEntry`, which sets the status for those not already waiting (EMOD
  refuses a second entry, Max_Entries 1), then the delay as `when=ss.on(..., delay=)`, whose event is the muxer's broadcast. Its
  Disqualifying_Properties are checked when due, as a delay does, and also on each step while waiting, which a delay
  does not: `CompiledCampaign.step()` ends the wait of those disqualified, as EMOD expires the muxer.
- Decisions on the person's state (HIVRapidHIVDiagnostic, STIIsPostDebut, HIVDrawBlood, the two ART staging
  diagnostics) → `Decide`, a product with outcomes `positive` and `negative`, using cascade.py's decision code.
- ARTMortalityTable, ARTDropout, PMTCT, MaleCircumcision, ControlledVaccine → `Act`, a product calling the HIV module;
  PMTCT's reset after 280 days is a delivery with a delay.
- A one-off event (Start_Year) → a delivery `when=` that year. NChooserEventCoordinatorHIV → a delivery at every step
  whose `who=` is cascade.py's NChooser selection.
- Events from coded modules (STIDebut, NewlySymptomatic, SixWeeksOld, TwelveWeeksPregnant) → deliveries at every step
  with `who=` (the last two a step late, as EMOD first updates what they give on the next step).

Deliveries are given in the order EMOD would update them for one person (breadth first along the chains, in campaign
order), which Starsim keeps where its triggers allow, so that, e.g., ART staging's status disqualifies the
combination-PMTCT decision as in EMOD.
"""
from collections import defaultdict
import numpy as np
import starsim as ss
from starsim.internals.conditions import Cond

from cascade import EmodCascade, sigmoid, piecewise, DELAYS

PEOPLE = ss.person


def ref(name, outcome=None):
    """ The event of a delivery (or one of its outcomes), by name, so a delivery can name one made later """
    return Cond('event', name, outcome)


def any_of(conds):
    out = None
    for c in conds:
        out = c if out is None else out | c
    return out


def all_of(conds):
    out = None
    for c in conds:
        if c is not None:
            out = c if out is None else out & c
    return out


class CampaignProduct(ss.Product):
    """ A product of the compiled campaign, acting through the state holder (`CompiledCampaign`) """
    def __init__(self, tid, name):
        super().__init__(name=name)
        self.tid = tid

    @property
    def cas(self):
        return self.sim.interventions.emodcascade


class SetIP(CampaignProduct):
    """ New_Property_Value, or a PropertyValueChanger: sets InterventionStatus (or ARTStatus, OnPrEP) """
    def __init__(self, tid, name, key=None, value=None):
        super().__init__(tid, name)
        self.key, self.value = key, value

    def administer(self, uids):
        cas = self.cas
        if self.key is None:
            cas.set_status(self.tid, uids)
        elif self.key == 'InterventionStatus':
            cas.istatus[uids] = cas.statuses.index(self.value)
        elif self.key == 'ARTStatus':
            cas.artstatus[uids] = cas.artstatuses.index(self.value)
        elif self.key == 'OnPrEP':
            cas.onprep[uids] = self.value == 'True'
        return


class Decide(CampaignProduct):
    """ A decision on the person's state: its status, then positive or negative (cascade.py's EmodCascade.decide) """
    def outcome_names(self):
        return ['positive', 'negative']

    def administer(self, uids):
        cas = self.cas
        cas.set_status(self.tid, uids)
        pos = cas.decide(cas.templates[self.tid], uids, cas.t.now('year'))
        return dict(positive=uids[pos], negative=uids[~pos])


class Act(CampaignProduct):
    """ What acts on the person (ART, dropout, PMTCT, circumcision, PrEP): cascade.py's EmodCascade.give, and the status """
    def administer(self, uids):
        cas = self.cas
        cas.give(self.tid, uids)
        cas.set_status(self.tid, uids)
        return


class MuxEntry(CampaignProduct):
    """ Entering a muxer: the status, for those not already waiting in it """
    def __init__(self, tid, name, wait):
        super().__init__(tid, name)
        self.wait = wait

    def administer(self, uids):
        waiting = ~np.isnan(self.sim.interventions[self.wait].ti_due.raw[uids])
        self.cas.set_status(self.tid, uids[~waiting])
        return


class ResetPMTCT(CampaignProduct):
    def administer(self, uids):
        self.sim.diseases.hiv.pmtct[uids] = 0.0
        return


class CompiledCampaign(EmodCascade):
    """
    The state of the compiled campaign (InterventionStatus, ARTStatus, OnPrEP, the programme's circumcisions) and what
    Starsim's deliveries don't do: who is dying this step, muxer waits ended by disqualification, the events of the last
    step that are given a step late, and the NChooser selections. Its deliveries are `compile_campaign()`'s.
    """
    def __init__(self, inp, **kwargs):
        super().__init__(inp, **kwargs)
        self.muxer_waits = {} # Name of the delayed delivery -> template, for the expiry of waits
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.dying = sim.people.add_state(ss.BoolArr('dying_hiv', default=False))
        self.prev_twelve = ss.uids()
        return

    def step(self):
        ppl = self.sim.people
        self.dt_days = self.t.dt_year*365
        self.dying[ppl.auids] = False
        self.dying[self.sim.diseases.hiv.dead_uids] = True
        for dname, tid in self.muxer_waits.items(): # EMOD expires a muxer whose person becomes disqualified while waiting
            d = self.sim.interventions[dname]
            u = ppl.auids[~np.isnan(d.ti_due.raw[ppl.auids])]
            gone = u[np.isin(self.istatus.raw[u], self.templates[tid]['disq']) | self.dying.raw[u]]
            d.ti_due[gone] = np.nan
        age = ppl.emod_age.raw[ppl.auids]
        self.six = ppl.auids[(age - 2*self.dt_days < 42) & (age - self.dt_days >= 42)] # Six weeks old on the step before
        self.twelve, self.prev_twelve = self.prev_twelve, self.sim.demographics.emoddemog.twelve_weeks
        self.year = self.t.now('year')
        return


def conditions(cas, target, tid, nodes=None):
    """ A listener's or coordinator's targeting, the template's Disqualifying_Properties, and not dying, as conditions on people """
    out = [~PEOPLE.dying_hiv]
    if nodes is not None and sorted(nodes) != [1, 2, 3]:
        out.append(any_of([PEOPLE.node == n - 1 for n in nodes]))
    td = target.get('Target_Demographic', 'Everyone')
    if td in ('ExplicitGender', 'ExplicitAgeRangesAndGender') and target.get('Target_Gender', 'All') != 'All':
        out.append(PEOPLE.female if target['Target_Gender'] == 'Female' else ~PEOPLE.female)
    if td == 'ExplicitAgeRangesAndGender':
        out += [PEOPLE.emod_age >= target.get('Target_Age_Min', 0)*365, PEOPLE.emod_age < target.get('Target_Age_Max', 125)*365]
    for d in target.get('Property_Restrictions_Within_Node') or []:
        assert list(d) == ['Accessibility'] and d['Accessibility'] == 'Yes', d # Everyone is Accessibility:Yes
    assert target.get('Demographic_Coverage', 1) == 1
    disq = cas.templates[tid]['disq']
    if disq:
        out.append(~any_of([PEOPLE.istatus == k for k in disq]))
    return all_of(out)


def compile_campaign(inp):
    """ The state holder and the deliveries of the cascade; returns (holder, [deliveries]) """
    cas = CompiledCampaign(inp)
    T = cas.templates
    specs = [] # Each: name, trigger (EMOD events, or another delivery), extra condition, kwargs of the Delivery, and what it broadcasts (outcome -> event)
    def spec(name, events, cond, emits=None, after=None, **kw):
        specs.append(dict(name=name, events=events, after=after, cond=cond, emits=emits or {}, kw=kw))
        return name

    # Events from coded modules, and the one-offs and NChoosers
    adapters = dict(STIDebut=lambda sim: sim.interventions.emodrisk.new_debuts, NewlySymptomatic=lambda sim: sim.diseases.hiv.newly_symptomatic,
                    SixWeeksOld=lambda sim: sim.interventions.emodcascade.six, TwelveWeeksPregnant=lambda sim: sim.interventions.emodcascade.twelve)
    for ev, fn in adapters.items():
        spec(f'ev_{ev}', None, None, emits={None: ev}, who=fn, when=each_step())
    for k, o in enumerate(cas.oneoffs):
        c = T[o['tid']]
        assert c['class'] == 'BroadcastEvent'
        spec(f'once{k}', None, conditions(cas, o['target'], o['tid'], o['nodes']), emits={None: c['Broadcast_Event']}, when=[o['start']], what=SetIP(o['tid'], f'once{k}_p'))
    for k, ch in enumerate(cas.choosers_state):
        spec(f'nchooser{k}', None, None, when=each_step(), who=lambda sim, k=k: sim.interventions.emodcascade.choosers_state[k].select(sim.interventions.emodcascade.year), what=Act(ch.cfg['tid'], f'nchooser{k}_p'))

    # Listeners
    for k, li in enumerate(cas.listeners):
        tid, c = li['tid'], T[li['tid']]
        cls, base = c['class'], f'L{k}'
        cond = conditions(cas, li['target'], tid, li['nodes'])
        start = dict(start=li['start']) if li['start'] > inp.start else {}
        ev = li['triggers']
        if cls in ('HIVSigmoidByYearAndSexDiagnostic', 'HIVPiecewiseByYearAndSexDiagnostic'):
            assert c.get('Female_Multiplier', 1) == 1
            if cls == 'HIVSigmoidByYearAndSexDiagnostic':
                p = lambda sim, c=c: max(0.0, float(sigmoid(sim.t.now('year'), c)))
            else:
                p = lambda sim, c=c: float(np.clip(piecewise(c['Time_Value_Map'], sim.t.now('year'), c.get('Default_Value', 0), c.get('Interpolation_Order', 0)), 0, 1))
            yes = spec(f'{base}_yes', ev, cond, emits={None: c.get('Positive_Diagnosis_Event')}, coverage=p, what=SetIP(tid, f'{base}_yes_p'), **start)
            spec(f'{base}_no', ev, cond & ~ref(yes), emits={None: c.get('Negative_Diagnosis_Event')}, what=SetIP(tid, f'{base}_no_p'), **start)
        elif cls == 'HIVRandomChoice':
            probs, names, before, left = c['Choice_Probabilities'], c['Choice_Names'], [], 1.0
            for j, (name, pr) in enumerate(zip(names, probs)):
                d = spec(f'{base}_c{j}', ev, all_of([cond] + [~ref(b) for b in before]), emits={None: name}, coverage=min(1.0, pr/left) if left > 0 else 0.0, what=SetIP(tid, f'{base}_c{j}_p'), **start)
                before.append(d)
                left -= pr
            if left > 1e-9:
                spec(f'{base}_rest', ev, all_of([cond] + [~ref(b) for b in before]), what=SetIP(tid, f'{base}_rest_p'), **start)
        elif cls in ('HIVRapidHIVDiagnostic', 'STIIsPostDebut', 'HIVDrawBlood', 'HIVARTStagingCD4AgnosticDiagnostic', 'HIVARTStagingByCD4Diagnostic'):
            spec(base, ev, cond, emits={'positive': c.get('Positive_Diagnosis_Event'), 'negative': c.get('Negative_Diagnosis_Event')}, what=Decide(tid, f'{base}_p'), **start)
        elif cls == 'PropertyValueChanger':
            spec(base, ev, cond, what=SetIP(tid, f'{base}_p', c['Target_Property_Key'], c['Target_Property_Value']), **start)
        elif cls == 'BroadcastEvent':
            spec(base, ev, cond, emits={None: c['Broadcast_Event']}, what=SetIP(tid, f'{base}_p'), **start)
        elif cls in ('ARTMortalityTable', 'ARTDropout', 'MaleCircumcision', 'ControlledVaccine', 'PMTCT'):
            d = spec(base, ev, cond, what=Act(tid, f'{base}_p'), **start)
            if cls == 'PMTCT': # Reset after 280 days: PMTCT::Update ends when its time left is within the step
                steps = int(np.ceil(280/inp.dt)) - 1
                spec(f'{base}_end', None, None, after=d, delay=ss.days(steps*inp.dt), what=ResetPMTCT(tid, f'{base}_end_p'))
        elif cls == 'HIVMuxer':
            wait = f'{base}_wait'
            assert c.get('Max_Entries', 1) == 1
            enter = spec(f'{base}_in', ev, cond, what=MuxEntry(tid, f'{base}_in_p', wait), **start)
            spec(wait, None, conditions(cas, {}, tid), after=enter, emits={None: c.get('Broadcast_Event')}, delay=delay_dist(c))
            cas.muxer_waits[wait] = tid
        else:
            raise NotImplementedError(cls)

    # Drop dead branches (events nothing broadcasts, and what only they trigger), then make the deliveries
    while True:
        alive = {x['name'] for x in specs}
        emitters = defaultdict(list)
        for x in specs:
            for outcome, event in x['emits'].items():
                if event:
                    emitters[event].append(ref(x['name'], outcome))
        keep = [x for x in specs if (x['events'] is None or any(e in emitters for e in x['events'])) and (x['after'] is None or x['after'] in alive)]
        if len(keep) == len(specs):
            break
        specs = keep
    trig = {}
    for x in specs:
        refs = [ref(x['after'])] if x['after'] else [r for e in (x['events'] or []) for r in emitters.get(e, [])]
        trig[x['name']] = refs
    seq = order(specs, trig)
    extra = after_previous(seq, specs, trig, cas)
    out = []
    for name in seq:
        x = next(y for y in specs if y['name'] == name)
        t = any_of(trig[name])
        if extra.get(name) and t is not None: # Runs after the deliveries before it in EMOD's order: event terms that select nobody
            t = t | any_of([ref(p) & ~ref(p) for p in extra[name]])
        kw = dict(x['kw'])
        trigger, cond = ss.conditions.split_events(all_of([t, x['cond']])) if t is not None or x['cond'] is not None else (None, None) # Its events trigger it, the rest says who is eligible
        delay = kw.pop('delay', None)
        when = kw.pop('when', None) or (ss.on(trigger, delay=delay) if delay is not None else trigger)
        out.append(ss.Delivery(name=name, when=when, who=kw.pop('who', None) or cond, disease='hiv', **kw))
    return cas, out


def order(specs, trig):
    """ Names of the deliveries breadth first along the chains of events, from those triggered by nothing, in campaign order """
    children = defaultdict(list)
    roots = []
    for x in specs:
        parents = {r.args[0] for r in trig[x['name']]} | {leaf.args[0] for leaf in ss.conditions.event_leaves(x['cond'])}
        if not parents:
            roots.append(x['name'])
        for p in parents:
            children[p].append(x['name'])
    seen, out, queue = set(), [], roots
    while queue:
        nxt = []
        for name in queue:
            if name not in seen:
                seen.add(name)
                out.append(name)
                nxt += children[name]
        queue = nxt
    return out + [x['name'] for x in specs if x['name'] not in seen]


def statuses(x, cas):
    """ The InterventionStatus a delivery sets (or None), and those it is disqualified by """
    p = x['kw'].get('what')
    sets = None
    if isinstance(p, SetIP) and p.key == 'InterventionStatus':
        sets = cas.statuses.index(p.value)
    elif isinstance(p, CampaignProduct) and not isinstance(p, ResetPMTCT):
        sets = cas.templates[p.tid]['npv']
    disq = set(cas.templates[p.tid]['disq']) if isinstance(p, CampaignProduct) else set()
    return sets, disq


def after_previous(seq, specs, trig, cas, depth=4):
    """
    Starsim runs a delivery after every delivery that can trigger it, and otherwise in the order given; EMOD updates one
    person's interventions breadth first along the chain they are in, and what one sets can disqualify the next. So
    along each chain (the deliveries a delivery's events lead to within the step, breadth first, in campaign order), a
    delivery waits for an earlier one if they interact through InterventionStatus (one sets a status the other is
    disqualified by, or both set it), where that makes no cycle. E.g. after a positive test at ANC, ART staging's
    status, set by a grandchild through one branch, then disqualifies the decision on combination PMTCT, a grandchild
    through the other, as in EMOD. Returns, for each delivery, the ones it must wait for.
    """
    delayed = {x['name'] for x in specs if 'delay' in x['kw']}
    by = {x['name']: x for x in specs}
    st = {n: statuses(by[n], cas) for n in seq}
    def interact(a, b):
        (sa, da), (sb, db) = st[a], st[b]
        return (sa is not None and sa in db) or (sb is not None and sb in da) or (sa is not None and sb is not None and sa != sb)
    leaves = lambda n: {r.args[0] for r in trig[n]} | {leaf.args[0] for leaf in ss.conditions.event_leaves(by[n]['cond'])}
    deps = {n: (leaves(n) if n not in delayed else set()) for n in seq}
    children = defaultdict(list)
    for n in seq:
        for p in deps[n]:
            children[p].append(n)
    def ancestors(n):
        stack, seen = list(deps[n]), set()
        while stack:
            m = stack.pop()
            if m not in seen:
                seen.add(m)
                stack += list(deps[m])
        return seen
    extra = defaultdict(list)
    for root in seq:
        level, chain = [root], []
        for _ in range(depth):
            level = [c for n in level for c in children[n] if c not in chain]
            chain += level
        for j, name in enumerate(chain):
            for prev in chain[:j]:
                if interact(prev, name) and prev not in ancestors(name) and name not in ancestors(prev):
                    deps[name].add(prev)
                    extra[name].append(prev)
    return dict(extra)


def each_step():
    """ Rounds at every step (EMOD′'s steps are months), for a delivery to those a function of the sim chooses """
    return ss.schedule(every=ss.months(1))


def delay_dist(c):
    d = c.get('Delay_Period_Distribution', 'CONSTANT_DISTRIBUTION')
    if d == 'CONSTANT_DISTRIBUTION':
        return ss.days(c['Delay_Period_Constant'])
    if d == 'UNIFORM_DISTRIBUTION':
        return ss.uniform(ss.days(c['Delay_Period_Min']), ss.days(c['Delay_Period_Max']))
    if d == 'EXPONENTIAL_DISTRIBUTION':
        return ss.expon(ss.days(c['Delay_Period_Exponential']))
    raise NotImplementedError(d)
