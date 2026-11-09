"""
EMOD-HIV Malawi inputs for one of Daniel's parameter sets, as the parameters of the EMOD′ network, risk and
demography modules. The inputs are built exactly as for the EMOD run (reference/run_emod.py), so his quirks
(e.g. the node order of the Risk overlay) come with them.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import sciris as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference'))
import run_emod # noqa: E402

TYPES = ['TRANSITORY', 'INFORMAL', 'MARITAL', 'COMMERCIAL'] # EMOD's RelationshipType order
RISKS = ['LOW', 'MEDIUM', 'HIGH']
SEXES = ['Male', 'Female'] # EMOD's Gender: 0 male, 1 female
DT = 365/12 # EMOD's Simulation_Timestep, days


def build(pid):
    """ The EMOD input files (as dicts) for parameter set `pid` """
    sets = pd.read_csv(run_emod.orig / 'resampled_parameter_sets_20230424.csv')
    row = sets[sets.parameterization_id == pid].iloc[0]
    meta = ['iteration_number', 'likelihood', 'parameterization_id', 'run_number', 'sim_id']
    table = run_emod.map_sample(row.drop(meta).to_dict(), run_emod.load_param_map())
    return run_emod.build_inputs(table, row.run_number), int(row.run_number)


def society(pfa, assort, node_id):
    """ The Society block for one node: the defaults with the node's values on top, then the assortativity overlay """
    soc = sc.dcp(pfa['Defaults']['Society'])
    node = next(n for n in pfa['Nodes'] if n['NodeID'] == node_id)
    soc = sc.mergenested(soc, node.get('Society', {}))
    return sc.mergenested(soc, assort['Defaults']['Society'])


def edges(pfp, sex):
    """ Upper bin edges of a PFA's age bins (years): age below edges[0] is bin 0, at or above edges[-1] the last bin """
    first, inc, n = pfp[f'Age_of_First_Bin_Edge_{sex}'], pfp[f'Years_Between_Bin_Edges_{sex}'], pfp[f'Number_Age_Bins_{sex}']
    return first + inc*np.arange(n)


def table2d(dist):
    """ An EMOD 2-axis demographics distribution (age, year): the groups of each axis and the values (age × year) """
    g = dist['PopulationGroups']
    return sc.objdict(age=np.array(g[0], float), year=np.array(g[1], float), vals=np.array(dist['ResultValues'], float)*dist['ResultScaleFactor']*365) # Per year


def campaign_sexwork(camp):
    """ The sex-work events of the campaign: who will ever sell or buy sex, by node and sex, and the delays """
    out = sc.objdict(p_ever=np.zeros((3, 2)), start_year=None, delay_max=np.zeros(2), dropout={})
    for ev in camp['Events']:
        ic = ev['Event_Coordinator_Config'].get('Intervention_Config', {})
        act = ic.get('Actual_IndividualIntervention_Config', {})
        trig = ic.get('Trigger_Condition_List', [])
        sex = {'Male': 0, 'Female': 1}.get(ic.get('Target_Gender'))
        if act.get('Positive_Diagnosis_Event') == 'Commercial_DelayFromDebutToUptake':
            assert trig == ['STIDebut'] and len(act['Time_Value_Map']['Values']) == 1
            for node in ev['Nodeset_Config']['Node_List']:
                out.p_ever[node-1, sex] = act['Time_Value_Map']['Values'][0]
        elif trig == ['Commercial_DelayFromDebutToUptake']: # Listeners that start the delay to uptake; they start in 1961
            assert act['Delay_Period_Distribution'] == 'UNIFORM_DISTRIBUTION' and act['Delay_Period_Min'] == 0
            out.delay_max[sex] = act['Delay_Period_Max']
            out.start_year = ev['Start_Year']
        elif trig == ['Commercial_Uptake'] and act['class'] == 'HIVDelayedIntervention':
            out.dropout[sex] = {k.replace('Delay_Period_', ''): v for k, v in act.items() if k.startswith('Delay_Period_')}
    assert out.dropout[0]['Distribution'] == 'UNIFORM_DISTRIBUTION' and out.dropout[1]['Distribution'] == 'WEIBULL_DISTRIBUTION'
    return out


def campaign_hiv(camp):
    """ The campaign's STI co-infection (at debut, by Risk), traditional circumcision (by node) and seeding """
    out = sc.objdict(p_sti=np.zeros(3), p_circ=np.zeros(3), circ_year=None, seed_year=None, seed_coverage=None)
    for ev in camp['Events']:
        ec = ev['Event_Coordinator_Config']
        ic = ec.get('Intervention_Config', {})
        act = ic.get('Actual_IndividualIntervention_Config', {})
        if ic.get('Trigger_Condition_List') == ['STIDebut'] and act.get('class') == 'ModifyStiCoInfectionStatus':
            risks = [r['Risk'] for r in ic.get('Property_Restrictions_Within_Node', [])] or [r.split(':')[1] for r in ic['Property_Restrictions']]
            for r in risks:
                out.p_sti[RISKS.index(r)] = ic['Demographic_Coverage']
        elif ic.get('Trigger_Condition_List') == ['Births'] and act.get('Choice_Names', [''])[0] == 'WillReceiveTraditionalMaleCircumcision':
            for node in ev['Nodeset_Config']['Node_List']:
                out.p_circ[node-1] = act['Choice_Probabilities'][0]
            out.circ_year = ev['Start_Year']
        elif ic.get('class') == 'OutbreakIndividual':
            assert ec['Property_Restrictions'] == ['Risk:HIGH'] and ic['Incubation_Period_Override'] == 0
            out.seed_year, out.seed_coverage = ev['Start_Year'], ec['Demographic_Coverage']
        elif ic.get('Trigger_Condition_List') == ['WillReceiveTraditionalMaleCircumcision']:
            out.circ_efficacy = act['Circumcision_Reduced_Acquire']
    assert out.p_sti.all() and out.p_circ.all() and out.seed_year
    return out


FIXES = dict( # Fixes of the published model's input files, by name (see fix()); none by default
    mortality='Demographics_Malawi.json: the mortality tables\' ResultScaleFactor 1/365 (deaths per person per year, as labelled), not 1/365,000',
    reenrol='the campaign: on OnART3 (ART dropout), InterventionStatus to None from 1990, so those who drop out can start ART again',
    art_male='the campaign: the male non-suppressive ART mortality table\'s 182-365-day rows from CD4 25 (about 7 times the female ones) as the female table has them',
    risk_nodes='the Risk overlay: each region\'s LOW-risk fraction on its own node (the overlay lists them as Southern, Northern, Central)',
    # Fixes of how EMOD reads the campaign, in the condensed cascade only (extract_cascade.py; EmodCascade refuses them)
    anc='ANC testing uptake from its Time_Value_Map (EMOD\'s sigmoid class ignores the table and uses its default curve)',
    pmtct='combination PMTCT (Option B) not blocked for women sent to ART staging by the same positive ANC test',
)
REENROL = { # The `reenrol` fix's campaign event, as in Daniel's fix (hivsim_paper PR #4)
    'Event_Name': 'OnART: switch InterventionStatus off so people can reenroll', 'class': 'CampaignEventByYear',
    'Nodeset_Config': {'class': 'NodeSetAll'}, 'Start_Year': 1990,
    'Event_Coordinator_Config': {'class': 'StandardInterventionDistributionEventCoordinator', 'Intervention_Config': {
        'class': 'NodeLevelHealthTriggeredIV', 'Trigger_Condition_List': ['OnART3'],
        'Actual_IndividualIntervention_Config': {'class': 'PropertyValueChanger', 'Target_Property_Key': 'InterventionStatus', 'Target_Property_Value': 'None'}}},
}


def fix(files, fixes):
    """ The EMOD input files with the named FIXES applied (in place) """
    assert set(fixes) <= set(FIXES), f'unknown fixes {set(fixes) - set(FIXES)}'
    if 'mortality' in fixes:
        ia = files['Demographics_Malawi.json']['Defaults']['IndividualAttributes']
        for k in ('MortalityDistributionMale', 'MortalityDistributionFemale'):
            assert ia[k]['ResultScaleFactor'] == 2.7397260273972604e-06
            ia[k]['ResultScaleFactor'] = 2.7397260273972604e-03
    events = files[run_emod.CAMPAIGN]['Events']
    if 'reenrol' in fixes: # Where Daniel put it: before "OnART: switch ARTStatus off"
        events.insert(next(i for i, ev in enumerate(events) if ev.get('Event_Name') == 'OnART: switch ARTStatus off'), sc.dcp(REENROL))
    if 'art_male' in fixes:
        tables = {}
        for ev in events:
            ic = ev['Event_Coordinator_Config'].get('Intervention_Config', {})
            act = ic.get('Actual_IndividualIntervention_Config', {})
            if act.get('class') == 'ARTMortalityTable' and 'NonSuppressiveART' in ic['Trigger_Condition_List']:
                tables[ic['Target_Gender']] = act['MortalityTable']
        assert tables['Male'][1][0][2] == 0.56
        for m, f in zip(tables['Male'][1], tables['Female'][1]): # By age; CD4 from 25 (the row below 25 matches)
            m[2:] = f[2:]
    if 'risk_nodes' in fixes:
        nodes = {n['NodeID']: n for n in files['Accessibility_and_Risk_IP_Overlay_Malawi.json']['Nodes']}
        risk = {i: next(ip for ip in n['IndividualProperties'] if ip['Property'] == 'Risk') for i, n in nodes.items()}
        dist = {i: r['Initial_Distribution'] for i, r in risk.items()}
        for node, listed in ((1, 2), (2, 3), (3, 1)): # NodeID 1 is Northern in every other file; the overlay's 2 is its Northern, 3 Central, 1 Southern
            risk[node]['Initial_Distribution'] = dist[listed]
    return files


def load(pid, fixes=()):
    """ Parameters of the EMOD′ network, risk and demography for parameter set `pid`, with the named FIXES of its input files """
    files, run_number = build(pid)
    if fixes:
        files = fix(files, fixes)
    cp = files['config.json']['parameters']
    pfa, assort = files['PFA_Overlay_Malawi.json'], files['Risk_Assortivity_Overlay_Malawi.json']
    socs = [society(pfa, assort, n) for n in (1, 2, 3)]
    p = sc.objdict(pid=pid, run_number=run_number, dt=cp['Simulation_Timestep'], start=cp['Base_Year'], scale=cp['x_Base_Population'])
    if fixes:
        p.fixes = tuple(fixes)
    assert abs(p.dt - DT) < 1e-6

    # Network: one entry per relationship type, in EMOD's order
    s0 = socs[0]
    p.types = []
    for t in TYPES:
        pfp, rp = s0[t]['Pair_Formation_Parameters'], s0[t]['Relationship_Parameters']
        same = lambda s: {k: v for k, v in s[t]['Relationship_Parameters'].items() if k != 'Condom_Usage_Probability'}
        assert all(s[t]['Pair_Formation_Parameters'] == pfp and same(s) == same(s0) for s in socs) # Everything but condom use is the same in every node
        conc = s0[t]['Concurrency_Parameters']
        joint = np.array(pfp['Joint_Probabilities'], float)
        rows = joint.sum(axis=1, keepdims=True)
        assert pfp['Formation_Rate_Type'] == 'CONSTANT'
        p.types.append(sc.objdict(
            name = t,
            rate = pfp['Formation_Rate_Constant'], # Per day
            edges = [edges(pfp, 'Male'), edges(pfp, 'Female')],
            marginal = [joint.sum(axis=1)/joint.sum(), joint.sum(axis=0)/joint.sum()],
            joint = np.divide(joint, rows, out=np.zeros_like(joint), where=rows > 0), # Rows normalized, as EMOD uses them for matching
            rate_ratio = [pfp['Extra_Relational_Rate_Ratio_Male'], pfp['Extra_Relational_Rate_Ratio_Female']],
            update_period = pfp['Update_Period'],
            assort = np.array(pfp['Assortivity']['Weighting_Matrix_RowMale_ColumnFemale'], float),
            dur_scale = rp['Duration_Weibull_Scale'], dur_het = rp['Duration_Weibull_Heterogeneity'],
            coital_rate = rp['Coital_Act_Rate'],
            condoms = [s[t]['Relationship_Parameters']['Condom_Usage_Probability'] for s in socs], # By node
            p_extra = np.array([[conc[r][f'Prob_Extra_Relationship_{x}'] for x in SEXES] for r in RISKS]), # Risk × sex
            max_rels = np.array([[conc[r][f'Max_Simultaneous_Relationships_{x}'] for x in SEXES] for r in RISKS]),
        ))
        assert pfp['Assortivity']['Group'] == 'INDIVIDUAL_PROPERTY' and pfp['Assortivity']['Axes'] == RISKS
    cc = s0['Concurrency_Configuration']
    assert cc['Individual_Property_Name'] == 'Risk' and cc['Probability_Person_Is_Behavioral_Super_Spreader'] == 0
    p.flag_order = [[TYPES.index(t) for t in cc[r]['Correlated_Relationship_Type_Order']] for r in RISKS]
    assert all(cc[r]['Extra_Relational_Flag_Type'] == 'Correlated' for r in RISKS)
    p.burnin = cp['PFA_Burnin_Duration_In_Days']
    p.threshold = cp['PFA_Cum_Prob_Selection_Threshold']
    p.dilution = [1, 1, cp['Coital_Dilution_Factor_2_Partners'], cp['Coital_Dilution_Factor_3_Partners'], cp['Coital_Dilution_Factor_4_Plus_Partners']] if cp['Enable_Coital_Dilution'] else [1]*5
    p.debut = sc.objdict(min=cp['Sexual_Debut_Age_Min'], scale=[cp[f'Sexual_Debut_Age_{x}_Weibull_Scale'] for x in SEXES],
                         het=[cp[f'Sexual_Debut_Age_{x}_Weibull_Heterogeneity'] for x in SEXES])
    assert cp['Sexual_Debut_Age_Setting_Type'] == 'WEIBULL' and cp['Min_Days_Between_Adding_Relationships'] == 0

    # Risk: the initial distribution by node, as the overlay has it (by NodeID, whatever its NodeName says)
    ips = files['Accessibility_and_Risk_IP_Overlay_Malawi.json']['Nodes']
    p.p_low = np.zeros(3)
    for n in ips:
        risk = next(ip for ip in n['IndividualProperties'] if ip['Property'] == 'Risk')
        assert risk['Values'] == RISKS and risk['Initial_Distribution'][2] == 0
        p.p_low[n['NodeID'] - 1] = risk['Initial_Distribution'][0]
    p.sexwork = campaign_sexwork(files[run_emod.CAMPAIGN])
    p.campaign = campaign_hiv(files[run_emod.CAMPAIGN])
    p.events = files[run_emod.CAMPAIGN]['Events'] # The whole campaign, for the cascade (emod_cascade.py)
    ip = next(n for n in ips if n['NodeID'] == 1)['IndividualProperties']
    p.ip_values = {d['Property']: d['Values'] for d in ip}

    # HIV (InfectionHIV, SusceptibilityHIV, IndividualHumanSTI)
    keys = ['Base_Infectivity', 'Acute_Duration_In_Months', 'Acute_Stage_Infectivity_Multiplier', 'AIDS_Duration_In_Months', 'AIDS_Stage_Infectivity_Multiplier',
            'Condom_Transmission_Blocking_Probability', 'STI_Coinfection_Acquisition_Multiplier', 'STI_Coinfection_Transmission_Multiplier',
            'Male_To_Female_Relative_Infectivity_Ages', 'Male_To_Female_Relative_Infectivity_Multipliers', 'Maternal_Infection_Transmission_Probability',
            'CD4_Post_Infection_Weibull_Scale', 'CD4_At_Death_LogLogistic_Scale', 'HIV_Adult_Survival_Scale_Parameter_Intercept', 'HIV_Adult_Survival_Scale_Parameter_Slope',
            'HIV_Adult_Survival_Shape_Parameter', 'HIV_Age_Max_for_Adult_Age_Dependent_Survival', 'HIV_Age_Max_for_Child_Survival_Function',
            'HIV_Child_Survival_Rapid_Progressor_Fraction', 'HIV_Child_Survival_Rapid_Progressor_Rate', 'HIV_Child_Survival_Slow_Progressor_Scale', 'HIV_Child_Survival_Slow_Progressor_Shape',
            'Days_Between_Symptomatic_And_Death_Weibull_Scale', 'Days_Between_Symptomatic_And_Death_Weibull_Heterogeneity', 'Maternal_Transmission_ART_Multiplier']
    p.hiv = sc.objdict({k: cp[k] for k in keys})
    assert cp['Max_Individual_Infections'] == 1 and cp['Heterogeneous_Infectiousness_LogNormal_Scale'] == 0 and cp['Enable_Maternal_Infection_Transmission']
    assert cp['CD4_Post_Infection_Weibull_Heterogeneity'] == 0 and cp['CD4_At_Death_LogLogistic_Heterogeneity'] == 0 and cp['Base_Mortality'] == 0

    # Demography
    dem = files['Demographics_Malawi.json']
    ia = dem['Defaults']['IndividualAttributes']
    p.init_pop = np.array([n['NodeAttributes']['InitialPopulation'] for n in dem['Nodes']])*cp['x_Base_Population']
    p.age_cdf = sc.objdict(age=np.array(ia['AgeDistribution']['ResultValues'], float), cdf=np.array(ia['AgeDistribution']['DistributionValues'], float))
    p.fertility = table2d(ia['FertilityDistribution']) # Pregnancies per woman per year
    p.mortality = [table2d(ia['MortalityDistributionMale']), table2d(ia['MortalityDistributionFemale'])] # Deaths per person per year, as EMOD reads them
    assert cp['Birth_Rate_Dependence'] == 'INDIVIDUAL_PREGNANCIES_BY_AGE_AND_YEAR' and cp['Death_Rate_Dependence'] == 'NONDISEASE_MORTALITY_BY_YEAR_AND_AGE_FOR_EACH_GENDER'
    return p
