"""
EMOD′'s slot factories (docs/swap_design.md §1.1): `make_<slot>(inp, context)` gives that slot's part of the sim, as
keywords of `ss.Sim` (lists of modules, and the sim's own parameters), for one parameter set's inputs (`inputs.load(pid)`);
`make_sim()` puts them together in the canonical order. `harness.make_sim()` is built on these.

    context: seed, stop, record (the network's and HIV's logs); network (N and R), hiv (T, H, S), cascade (C: None,
             'interpreted' (cascade.py), 'compiled' (compile_campaign.py) or 'condensed' (common/cascade));
             analyzers, sim_kwargs (passed to ss.Sim)

Where one module covers several slots, the first slot it belongs to builds it, and the others give nothing:
T (transmission) and S (seeding) are inside H's `EmodHIV`.
"""
import sciris as sc
import starsim as ss

from demography import EmodDemography
from network import EmodNetwork
from risk import EmodRisk
from hiv import EmodHIV
from cascade import EmodCascade # At import, so that `cascade` is this folder's (common/cascade is a package of the same name)

SLOTS = ['D', 'R', 'N', 'T', 'H', 'S', 'C', 'M']
KINDS = ['demographics', 'networks', 'diseases', 'interventions', 'analyzers'] # The keywords that hold lists of modules


def make_D(inp, context):
    """ Demography: EMOD's population at 1% scale, from 1960.5 in monthly steps, ageing as EMOD does (not Starsim's) """
    return dict(n_agents=int(inp.init_pop.sum()), start=inp.start, dt=1/12, aging=False, demographics=EmodDemography(inp))


def make_R(inp, context):
    """ Risk and sex work, and the STI flag: EmodRisk (an intervention, first of them) """
    return dict(interventions=EmodRisk(inp)) if context.network or context.hiv or context.cascade else {}


def make_N(inp, context):
    """ Sexual contact: EmodNetwork """
    return dict(networks=EmodNetwork(inp, record=context.record)) if context.network or context.hiv or context.cascade else {}


def make_T(inp, context):
    """ Transmission: inside EmodHIV (make_H) """
    return {}


def make_H(inp, context):
    """ Natural history (and, as EmodHIV, transmission and the 1983 seeding) """
    return dict(diseases=EmodHIV(inp, record=context.record)) if context.hiv or context.cascade else {}


def make_S(inp, context):
    """ Seeding: inside EmodHIV (make_H) """
    return {}


def make_C(inp, context):
    """ Care cascade and prevention: the campaign interpreted (EmodCascade), compiled into deliveries, or as the condensed cascade """
    if context.cascade == 'compiled':
        from compile_campaign import compile_campaign
        holder, deliveries = compile_campaign(inp)
        return dict(interventions=[holder] + deliveries)
    elif context.cascade == 'condensed':
        from cascade_host import EmodHost
        from extract_cascade import extract
        from common.cascade.cascade import condensed_cascade
        from common.cascade.params import Params
        return dict(interventions=condensed_cascade(Params(*extract(inp)), EmodHost))
    elif context.cascade == 'interpreted':
        return dict(interventions=EmodCascade(inp))
    elif context.cascade:
        raise ValueError(f'Unknown cascade {context.cascade!r}')
    return {}


def make_M(inp, context):
    """ Measurement: Snapshots, and with a cascade, EmodReport and then any analyzers given; without, the analyzers given or else Snapshots """
    from harness import Snapshots, EmodReport
    if context.cascade:
        return dict(analyzers=[Snapshots(), EmodReport()] + sc.tolist(context.analyzers))
    return dict(analyzers=context.analyzers if context.analyzers is not None else Snapshots())


FACTORIES = dict(D=make_D, R=make_R, N=make_N, T=make_T, H=make_H, S=make_S, C=make_C, M=make_M)


def assemble(parts):
    """ The slots' parts as ss.Sim's keywords: module lists joined in slot order, other keywords given once """
    kw = {}
    for part in parts:
        for k, v in part.items():
            if k in KINDS and not isinstance(v, str): # A string (e.g. STIsim's location for demographics) is given as it is
                kw[k] = kw.get(k, []) + sc.tolist(v)
            elif k in kw:
                raise ValueError(f'Two slots set {k}')
            else:
                kw[k] = v
    return kw


def make_sim(inp, context=None, **kwargs):
    """ EMOD′ from its slot factories, for the inputs `inp` (inputs.load(pid)); `context` keys as in the module docstring, or as keywords """
    context = sc.objdict(dict(seed=1, stop=2025, record=True, network=True, hiv=False, cascade=None, analyzers=None, sim_kwargs={}) | dict(context or {}) | kwargs)
    kw = assemble(FACTORIES[s](inp, context) for s in SLOTS)
    return ss.Sim(**kw, stop=context.stop, rand_seed=context.seed, verbose=0, **context.sim_kwargs)
