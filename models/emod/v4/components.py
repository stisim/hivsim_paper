"""
EMOD′'s slot factories (docs/swap_design.md §1.1, in the runner's interface, common/swaps/README.md): `make_<slot>(inp, ctx)`
gives that slot's part of the sim as `ss.Sim` keywords (lists of modules, and the sim's own parameters), for one
parameter set's inputs (`pars(pid)`, i.e. `inputs.load(pid)`). `make_sim()` puts them together in the canonical order;
`harness.make_sim()` is built on it, and the runner's `build()` calls the same factories.

    ctx: from the runner, `choice`, `models`, `pars`, `seed`, `stop`, `cell`; from make_sim(), also: record (the
         network's and HIV's logs; default False), network (R and N), hiv (T, H, S), cascade (C: None, 'interpreted'
         (emod_cascade.py), 'compiled' (compile_campaign.py) or 'condensed' (common/cascade), the default), analyzers

Fine slots, and the runner's coarse ones made of them: RN = R + N, V = T + H + S. Where one module covers several
slots, the first slot it belongs to builds it and the others give nothing: T and S are inside H's `EmodHIV`.
"""
from types import SimpleNamespace
import sciris as sc
import starsim as ss

import inputs
from demography import EmodDemography
from network import EmodNetwork
from risk import EmodRisk
from hiv import EmodHIV
from emod_cascade import EmodCascade
from cascade_host import EmodHost

HOST = EmodHost # The CascadeHost, when EMOD′'s V is in the hybrid
SCALE = 100 # People per agent (1% scale)
SLOTS = ['D', 'RN', 'V', 'C', 'M'] # The order make_sim() calls them in (the runner's canonical order)
LISTS = ['demographics', 'networks', 'diseases', 'interventions', 'connectors', 'analyzers'] # The keywords that hold lists of modules


def pars(pid):
    return inputs.load(pid)


def opt(ctx, key, default):
    """ A make_sim() option; the runner gives none, so the defaults are the full model with the condensed cascade """
    return getattr(ctx, key, default)


def cascade(ctx):
    return opt(ctx, 'cascade', 'condensed')


def with_hiv(ctx):
    return opt(ctx, 'hiv', True) or bool(cascade(ctx))


def with_network(ctx):
    return opt(ctx, 'network', True) or with_hiv(ctx)


def make_D(inp, ctx):
    """ Demography: EMOD's population at 1% scale, from 1960.5 in monthly steps, ageing as EMOD does (not Starsim's) """
    return dict(n_agents=int(inp.init_pop.sum()), start=inp.start, dt=1/12, aging=False, demographics=[EmodDemography(inp)])


def make_R(inp, ctx):
    """ Risk and sex work, and the STI flag: EmodRisk (an intervention, first of them) """
    return dict(interventions=[EmodRisk(inp)]) if with_network(ctx) else {}


def make_N(inp, ctx):
    """ Sexual contact: EmodNetwork """
    return dict(networks=[EmodNetwork(inp, record=opt(ctx, 'record', False))]) if with_network(ctx) else {}


def make_T(inp, ctx):
    """ Transmission: inside EmodHIV (make_H) """
    return {}


def make_H(inp, ctx):
    """ Natural history (and, as EmodHIV, transmission and the 1983 seeding) """
    return dict(diseases=[EmodHIV(inp, record=opt(ctx, 'record', False))]) if with_hiv(ctx) else {}


def make_S(inp, ctx):
    """ Seeding: inside EmodHIV (make_H) """
    return {}


def make_C(inp, ctx):
    """ Care cascade and prevention: the campaign interpreted (EmodCascade), compiled into deliveries, or as the condensed
    cascade (on the host of V's model, in a hybrid); with EMOD's reports (Snapshots, EmodReport), as harness.make_sim has """
    cx = cascade(ctx)
    if not cx:
        return {}
    from harness import Snapshots, EmodReport
    if cx == 'compiled':
        from compile_campaign import compile_campaign
        holder, deliveries = compile_campaign(inp)
        intvs = [holder] + deliveries
    elif cx == 'condensed':
        from extract_cascade import extract
        from common.cascade.cascade import condensed_cascade
        from common.cascade.params import Params
        choice = opt(ctx, 'choice', {})
        host = ctx.models[choice['V']].HOST if 'V' in choice else HOST
        intvs = condensed_cascade(Params(*extract(inp)), host)
    elif cx == 'interpreted':
        intvs = [EmodCascade(inp)]
    else:
        raise ValueError(f'Unknown cascade {cx!r}')
    return dict(interventions=intvs, analyzers=[Snapshots(), EmodReport()])


def make_M(inp, ctx):
    """ Measurement: any analyzers given (after C's reports); without a cascade, the analyzers given or else Snapshots """
    from harness import Snapshots
    analyzers = opt(ctx, 'analyzers', None)
    if analyzers is None and not cascade(ctx):
        analyzers = Snapshots()
    return dict(analyzers=sc.tolist(analyzers))


def join(*parts):
    """ Slots' parts as one: module lists joined in order, other keywords given once """
    kw = {}
    for part in parts:
        for k, v in part.items():
            if k in LISTS:
                kw[k] = kw.get(k, []) + sc.tolist(v)
            elif k in kw:
                raise ValueError(f'Two slots set {k}')
            else:
                kw[k] = v
    return kw


def make_RN(inp, ctx):
    """ Risk and sexual contact (the runner's coarse slot) """
    return join(make_R(inp, ctx), make_N(inp, ctx))


def make_V(inp, ctx):
    """ The virus: transmission, natural history and seeding (the runner's coarse slot) """
    return join(make_T(inp, ctx), make_H(inp, ctx), make_S(inp, ctx))


def make_sim(inp, ctx=None, sim_kwargs=None, **kwargs):
    """ EMOD′ from its slot factories, for the inputs `inp` (pars(pid)); `ctx` keys as in the module docstring, or as keywords """
    ctx = SimpleNamespace(**(dict(seed=1, stop=2025) | dict(ctx or {}) | kwargs))
    kw = join(*[globals()[f'make_{s}'](inp, ctx) for s in SLOTS])
    return ss.Sim(**kw, stop=ctx.stop, rand_seed=ctx.seed, verbose=0, **(sim_kwargs or {}))
