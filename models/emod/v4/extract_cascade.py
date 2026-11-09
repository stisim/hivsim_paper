"""
The condensed cascade's parameter tables (common/cascade/params.py) for Malawi, read from Daniel's campaign for one
parameter set: each row names the campaign event it comes from.

    python extract_cascade.py --pid 0   # writes common/cascade/tables/malawi_emod_set000_{params,flows}.csv
"""
import argparse
from pathlib import Path
import pandas as pd

import inputs
from emod_cascade import EmodCascade

OUT = Path(__file__).resolve().parents[3] / 'common' / 'cascade' / 'tables'


def tvm(m):
    return ', '.join(f'{t:g}: {v:g}' for t, v in zip(m['Times'], m['Values']))


def extract(inp):
    """ (params, flows) DataFrames for one parameter set """
    cas = EmodCascade(inp)
    T = cas.templates
    def find(trigger, cls=None, **target):
        for k, li in enumerate(cas.listeners):
            c, t = T[li['tid']], li['target']
            if trigger in li['triggers'] and (cls is None or c['class'] == cls) and all(t.get(a) == b for a, b in target.items()):
                return c, f'listener on {trigger} ({c["class"]})'
        raise KeyError(trigger)
    rows = []
    def row(pathway, parameter, form, value, source, sex='all', default=None):
        rows.append(dict(pathway=pathway, parameter=parameter, sex=sex, form=form, value=value, default=default, source=source))
    def sig(c):
        return f'{c.get("Ramp_Min", 0):g}, {c.get("Ramp_Max", 1):g}, {c.get("Ramp_MidYear", 2000):g}, {c.get("Ramp_Rate", 1):g}'
    def pw(c):
        return ('linear' if c.get('Interpolation_Order', 0) == 1 else 'step'), tvm(c['Time_Value_Map'])
    def delay(c):
        d = c['Delay_Period_Distribution']
        if d == 'EXPONENTIAL_DISTRIBUTION':
            return 'exp_days', f'{c["Delay_Period_Exponential"]}'
        if d == 'UNIFORM_DISTRIBUTION':
            return 'uniform_days', f'{c["Delay_Period_Min"]}, {c["Delay_Period_Max"]}'
        return 'days', f'{c["Delay_Period_Constant"]}'
    one = next(o for o in cas.oneoffs if T[o['tid']]['class'] == 'BroadcastEvent')

    # Routine testing
    row('hct', 'start', 'const', one['start'], 'one-off BroadcastEvent of HCTUptakePostDebut0 (everyone past debut into uptake)')
    c, s = find('STIDebut', 'HIVSigmoidByYearAndSexDiagnostic'); row('hct', 'p_debut', 'sigmoid', sig(c), s)
    c, s = find('HCTUptakePostDebut2'); row('hct', 'p_uptake', 'sigmoid', sig(c), s)
    c, s = find('HCTUptakePostDebut1'); row('hct', 'uptake_interval', *delay(c), s)
    c, s = find('HCTTestingLoop0'); row('hct', 'test_interval', *delay(c), s)
    for sex, g in (('f', 'Female'), ('m', 'Male')):
        c, s = find('ARTStaging9', Target_Gender=g); row('hct', 'p_immediate_art', *pw(c), s + ', positive: OnART1 (no staging)', sex=sex, default=c.get('Default_Value', 0))

    # Testing on symptoms, at ANC, of infants
    c, s = find('NewlySymptomatic'); row('symptomatic', 'p_present', *pw(c), s, default=0)
    row('symptomatic', 'start', 'const', 1990, 'Start_Year of its listener')
    c, s = find('ARTStaging8'); row('symptomatic', 'p_present_again', *pw(c), s + ' (those who did not present go straight to staging)', default=0)
    fixes = getattr(inp, 'fixes', ())
    c, s = find('TwelveWeeksPregnant')
    if 'anc' in fixes:
        row('anc', 'p_test', *pw(c), s + ': its Time_Value_Map (fix anc; EMOD ignores it)', default=0)
    else:
        row('anc', 'p_test', 'sigmoid', sig(c), s + ': its Time_Value_Map is ignored by the class, so the default sigmoid applies (as in EMOD)')
    row('anc', 'start', 'const', next(li['start'] for li in cas.listeners if 'TwelveWeeksPregnant' in li['triggers']), 'Start_Year of its listener')
    c, s = find('HIV_Positive_at_ANC', 'HIVPiecewiseByYearAndSexDiagnostic'); row('anc', 'p_stage', *pw(c), s, default=0)
    c, s = find('SixWeeksOld'); row('infant', 'p_test', *pw(c), s, default=0)
    row('infant', 'start', 'const', next(li['start'] for li in cas.listeners if 'SixWeeksOld' in li['triggers']), 'Start_Year of its listener')

    # PMTCT
    c, s = find('HIV_Positive_at_ANC', 'HIVSigmoidByYearAndSexDiagnostic'); row('pmtct', 'p_combination', 'sigmoid', sig(c), s)
    c, s = find('Needs_Combination_PMTCT'); row('pmtct', 'p_option_b', *pw(c), s, default=c.get('Default_Value', 0))
    for name, trig in (('eff_sdnvp', 'Needs_sdNVP_PMTCT'), ('eff_option_a', 'Needs_Option_A'), ('eff_option_b', 'Needs_Option_B')):
        c, s = find(trig); row('pmtct', name, 'const', c['Efficacy'], s)
    row('pmtct', 'days', 'const', 280, 'PMTCT::Update: the efficacy ends after 40 weeks')
    if 'pmtct' in fixes:
        row('pmtct', 'blocked_if_staged', 'const', 0, 'fix pmtct: the combination decision is not disqualified by ARTStaging')
    else:
        row('pmtct', 'blocked_if_staged', 'const', 1, 'the combination decision is disqualified by ARTStaging, set first by staging (EMOD update order)')

    # Staging, eligibility, linkage, pre-ART
    c, s = find('ARTStaging1'); row('staging', 'p_lost', 'const', c['Choice_Probabilities'][c['Choice_Names'].index('HCTUptakePostDebut3')], s)
    c, s = find('ARTStaging2'); row('staging', 'delay', *delay(c), s + ' (fires the step after)')
    ag, s1 = find('ARTStaging4'); cd, s2 = find('ARTStaging6')
    row('eligibility', 'adult_age', 'const', ag['Adult_Treatment_Age'], s1)
    row('eligibility', 'who_adult', 'step', tvm(ag['Adult_By_WHO_Stage']), s1, default=6)
    row('eligibility', 'pregnant', 'step', tvm(ag['Adult_By_Pregnant']), s1, default=0)
    row('eligibility', 'child_age', 'step', tvm(ag['Child_Treat_Under_Age_In_Years_Threshold']), s1, default=-1)
    row('eligibility', 'who_child', 'step', tvm(ag['Child_By_WHO_Stage']), s1, default=6)
    row('eligibility', 'cd4', 'step', tvm(cd['Threshold']), s2, default=-1e9)
    row('eligibility', 'cd4_pregnant', 'step', tvm(cd['If_Pregnant']), s2, default=-1e9)
    c, s = find('LinkingToART0'); row('staging', 'p_link_art', 'sigmoid', sig(c), s)
    c, s = find('LinkingToPreART0'); row('staging', 'p_link_preart', 'sigmoid', sig(c), s)
    c, s = find('OnPreART0'); row('preart', 'interval', *delay(c), s)
    c, s = find('OnPreART1'); row('preart', 'p_lost', 'const', c['Choice_Probabilities'][c['Choice_Names'].index('HCTUptakePostDebut9')], s + ' (else eligibility again: the pre-ART recheck, whose rule differs a little from staging\'s; the staging rule is used)')

    # ART, dropout and return
    c, s = find('OnART0'); row('art', 'start_delay', *delay(c), s)
    row('art', 'age_split', 'const', 35, 'Target_Age_Min/Max of the listeners on OnART1')
    for sex, g in (('f', 'Female'), ('m', 'Male')):
        for name, lo in (('p_suppressive_under', 0), ('p_suppressive_over', 35)):
            c, s = find('OnART1', 'HIVRandomChoice', Target_Gender=g, Target_Age_Min=lo)
            row('art', name, 'const', c['Choice_Probabilities'][c['Choice_Names'].index('EffectiveART')], s, sex=sex)
    c, s = find('OnART1', 'HIVMuxer'); row('art', 'dropout_delay', *delay(c), s)
    c, s = find('OnART3', 'HIVRandomChoice'); lost = c['Choice_Probabilities'][c['Choice_Names'].index('LostForever9')]
    c2, s2 = find('LostForever9'); row('dropout', 'p_lost', 'const', lost*(1 - piecewise0(c2)), s + ' × (1 - ' + s2 + ')')
    c, s = find('HCTUptakePostDebut8'); row('dropout', 'p_out', *pw(c), s + ': positive is OnART1, whose start of ART is disqualified by InterventionStatus OnEffective/NonSuppressiveART (the dropout trap); negative, back to uptake', default=0)
    if any(T[li['tid']]['class'] == 'PropertyValueChanger' and T[li['tid']]['Target_Property_Key'] == 'InterventionStatus' for li in cas.listeners if 'OnART3' in li['triggers']):
        row('dropout', 'p_restart', 'const', 1, 'listener on OnART3 (PropertyValueChanger, fix reenrol): InterventionStatus None at dropout, so p_out start ART again at OnART1')
    else:
        row('dropout', 'p_restart', 'const', 0, 'the dropout trap: those who decide to restart (p_out) get no ART, and are offered dropout and return again')

    # Prevention
    pc = next(ch for ch in cas.choosers if T[ch['tid']]['class'] == 'ControlledVaccine')
    w = T[pc['tid']]['Waning_Config']
    row('prep', 'efficacy', 'const', w['Initial_Effect'], 'ControlledVaccine (oral PrEP)')
    row('prep', 'days', 'const', w['Durability_Map']['Times'][-2], 'its Durability_Map')
    row('flows', 'scale', 'const', inp.scale, 'x_Base_Population: targets are scaled and truncated')
    params = pd.DataFrame(rows)

    flows = []
    for ch in cas.choosers:
        kind = 'vmmc' if T[ch['tid']]['class'] == 'MaleCircumcision' else 'prep'
        for d in ch['ec']['Distributions']:
            restrict = ['not:vmmc'] if d.get('Target_Disease_State_Has_Intervention_Name') == 'Any_MC' else []
            within = d.get('Property_Restrictions_Within_Node') or []
            if within:
                restrict.append('risk:' + '|'.join(r['Risk'] for r in within))
            for k, ar in enumerate(d['Age_Ranges_Years']):
                for sex, col in (('m', 'Num_Targeted_Males'), ('f', 'Num_Targeted_Females')):
                    if d[col][k]:
                        flows.append(dict(intervention=kind, node=ch['nodes'][0] if len(ch['nodes']) == 1 else 'all', sex=sex, age_lo=ar['Min'], age_hi=ar['Max'],
                                          start=d['Start_Year'], end=d['End_Year'], n=d[col][k], restrict=';'.join(restrict)))
    return params, pd.DataFrame(flows)


def piecewise0(c):
    """ The value of a piecewise map at its first time (e.g. a probability that is 0 throughout) """
    return c['Time_Value_Map']['Values'][0]


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pid', type=int, default=0)
    a = ap.parse_args()
    p, f = extract(inputs.load(a.pid))
    OUT.mkdir(parents=True, exist_ok=True)
    p.to_csv(OUT / f'malawi_emod_set{a.pid:03d}_params.csv', index=False)
    f.to_csv(OUT / f'malawi_emod_set{a.pid:03d}_flows.csv', index=False)
    print(p.to_string())
