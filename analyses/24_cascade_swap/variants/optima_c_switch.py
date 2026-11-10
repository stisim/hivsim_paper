"""
Optima's cascade (models/optima/v4/components.py's make_C) with treatment failure and regain acting on the host's ART
arm: `switch_arm` as the host's `start_art` in the new arm. On EMOD′ that is EmodHIV.go_on_art for someone already on
ART (EMOD's GoOnART: the new arm's table and transmission multiplier, a new prognosis, time on ART kept). With
PREVENTION (optima_c_switch_prev.py), EMOD's own VMMC and PrEP flows (its campaign's, through the condensed cascade's
holder) run as well. With PMTCT (optima_c_pmtct.py, optima_c_switch_prev_pmtct.py), the table's PMTCT number is drawn
per person (`TableCascade(pmtct=True)`); with SWITCH False, the arm doesn't switch (as models/optima/v4/components.py).
Sensitivity runs of analysis 24 only.
"""
import sys
import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('optima_components', REPO / 'models' / 'optima' / 'v4' / 'components.py')
oc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oc)
PREVENTION = False
SWITCH = True
PMTCT = True # As models/optima/v4/components.py (False: skipped, as the stored runs before 10 Nov)


def pars(pid=None):
    return oc.pars(pid)


def switching(host):
    class Switching(host):
        def switch_arm(self, uids, arm):
            self.start_art(uids, arm)
    Switching.__name__ = f'{host.__name__}Switching'
    return Switching


def make_C(params, ctx):
    from common.cascade.pathways import TableCascade
    vm = ctx.models[ctx.choice['V']]
    scale = getattr(ctx.models[ctx.choice['D']], 'SCALE', None)
    intvs = [TableCascade(params, switching(vm.HOST) if SWITCH else vm.HOST, scale=scale, pmtct=PMTCT)]
    if PREVENTION:
        intvs = prevention(ctx, vm.HOST) + intvs
    return dict(interventions=intvs)


def prevention(ctx, host):
    """ EMOD′'s VMMC and PrEP flows for the cell's parameter set, as condensed_cascade() runs them """
    import starsim as ss
    from common.cascade.params import Params
    from common.cascade.cascade import Holder, Circumcise, PrEP, each_step
    emod = REPO / 'models' / 'emod' / 'v4'
    sys.path.insert(0, str(emod))
    spec = importlib.util.spec_from_file_location('emod_extract_cascade', emod / 'extract_cascade.py')
    ex = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ex)
    params, flows = ex.extract(ctx.pars[ctx.choice['V']])
    D = [Holder(Params(params, flows), host)]
    for name, prod in (('cx_vmmc', Circumcise('cx_vmmc_p')), ('cx_prep', PrEP('cx_prep_p'))):
        kind = name.split('_')[1]
        D.append(ss.Delivery(prod, name=name, when=each_step(), who=lambda sim, k=kind: sim.interventions.cascade.flow(k), disease='hiv'))
    return D
