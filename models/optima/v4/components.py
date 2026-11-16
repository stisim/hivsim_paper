"""
Optima′'s slot factories for the swap runner (common/swaps/README.md): Optima′ as agent modules (modules.py,
analysis 26), `make_<slot>(params, ctx)` giving each slot's part of the sim as `ss.Sim` keywords, and `make_sim()`
putting them together in the runner's canonical order, so `build(all Optima′)` is `make_sim()`.

- D: `OptimaDR`, Optima's demography and its 39 populations, with R (risk transitions between them): on Optima paths
  D and R are one step (docs/swap_design.md §1.4), so `make_R` gives nothing.
- N: none. Optima's transmission has no network: N and T are one step (D12), so `make_N` gives nothing and T is
  `OptimaT`, the force of infection by population pair. `make_RN` gives nothing either.
- T: `OptimaT`; H: `OptimaH`; S: `OptimaS` (Optima's 1990 state); V: T, H and S (the runner's coarse slot).
- C: Optima's cascade as a table (common/cascade/tables/malawi_optima_params.csv, analysis 23), run as agents by
  `common.cascade.pathways.TableCascade` on the host of H's model (or V's), with the people per agent of D's model
  (its `SCALE`) for the table's numbers (number on ART, viral-load monitoring).

On another model's people (analysis 23's provisional decisions): Optima's groups by population come from the host's
`population()` (poplabel); care and lost live only in the cascade's record where the host keeps none (EMOD′, HIVsim),
and `switch_arm` changes nothing there; the PrEP and PEP coverages are not applied (no per-person form). Optima's PMTCT
number is drawn per person (`TableCascade(pmtct=True)`: diagnosed pregnant women not on ART, Optima's share), on by
default (`PMTCT`; False skips it, as the runs before 10 Nov 2026 did). It adds nothing on Optima′'s own agents, whose
births have no mothers (no one is pregnant; Optima′ applies its PMTCT inside births, here in D).

Optima′ declared (one disease) runs from `optima_prime.make_sim` (ODE) or `agents.make_sim` (agents); the modules
here are checked against the latter in analysis 26.

The runner puts this folder on sys.path: models/emod/v4 has a module of the same name as this folder's
`extract_cascade`, so code of another model in the same hybrid should not import it flat after this file is loaded.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TABLE = REPO / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv'
PMTCT = True # Optima's PMTCT number per person on the host's pregnant women
N_AGENTS = 100_000 # Agents in 1990 (about 1% of Malawi), as analyses 20-24
DT = 0.02 # Years, as Optima′ as agents (agents.make_sim)
START = 1990
SLOTS = ['D', 'S', 'T', 'H', 'C'] # The order make_sim() calls them in (the runner's canonical order)
LISTS = ['demographics', 'networks', 'diseases', 'interventions', 'connectors', 'analyzers']


def __getattr__(name):
    if name == 'HOST': # The CascadeHost, when Optima′'s H (or V) is in the hybrid (imported only when asked for)
        from modules import ModulesHost
        return ModulesHost
    raise AttributeError(name)


def pars(pid=None):
    """ Optima's Malawi project has one parameter set: its cascade table (the other slots read the project's inputs, models/optima/reference) """
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from common.cascade.params import Params
    return Params(TABLE)


def host_of(ctx, slot='H'):
    """ The CascadeHost of the model filling H (or V) in the hybrid; Optima′'s without a choice """
    choice = getattr(ctx, 'choice', None) or {}
    m = choice.get(slot, choice.get('V'))
    if m is None:
        from modules import ModulesHost
        return ModulesHost
    return ctx.models[m].HOST


def make_D(params, ctx):
    """ Optima's demography and populations (D and R) """
    import modules as om
    total = om.ref().raw['people'][:, :, 0].sum()
    return dict(n_agents=getattr(ctx, 'n_agents', N_AGENTS), start=START, dt=DT, total_pop=total, demographics=[om.OptimaDR(host_of(ctx))])


def make_R(params, ctx):
    """ Nothing: Optima's risk transitions are part of its populations (D) """
    return {}


def make_N(params, ctx):
    """ Nothing: Optima's transmission has no network (T) """
    return {}


def make_RN(params, ctx):
    return {}


def make_T(params, ctx):
    """ Optima's force of infection by population pair """
    import modules as om
    return dict(diseases=[om.OptimaT(host=host_of(ctx))])


def make_H(params, ctx):
    """ Optima's natural history and ART biology """
    import modules as om
    return dict(diseases=[om.OptimaH(host=host_of(ctx))])


def make_S(params, ctx):
    """ Optima's state in 1990 """
    import modules as om
    return dict(interventions=[om.OptimaS(host=host_of(ctx))])


def make_V(params, ctx):
    """ T, H and S (the runner's coarse slot) """
    return join(make_T(params, ctx), make_H(params, ctx), make_S(params, ctx))


def make_C(params, ctx):
    """ Optima's cascade as a table, on the host of H's (or V's) model, after what that model's `for_cascade()` gives """
    from common.cascade.pathways import TableCascade
    choice = getattr(ctx, 'choice', None) or {}
    key = 'V' if 'V' in choice else 'H'
    vm = ctx.models[choice[key]] if key in choice else None
    pre = vm.for_cascade(ctx.pars[choice[key]], ctx) if hasattr(vm, 'for_cascade') else {}
    scale = getattr(ctx.models[choice['D']], 'SCALE', None) if 'D' in choice else None
    host = vm.HOST if vm is not None else host_of(ctx)
    intvs = list(pre.get('interventions', [])) + [TableCascade(params, host, scale=scale, pmtct=PMTCT)]
    return dict(pre, interventions=intvs)


def join(*parts):
    """ Slots' parts as one: module lists joined in order, other keywords given once """
    kw = {}
    for part in parts:
        for k, v in part.items():
            if k in LISTS:
                kw[k] = kw.get(k, []) + list(v if isinstance(v, (list, tuple)) else [v])
            elif k in kw:
                raise ValueError(f'Two slots set {k}')
            else:
                kw[k] = v
    return kw


def make_sim(seed=1, stop=2021, n_agents=N_AGENTS, outputs=True, analyzers=(), **kwargs):
    """ Optima′ as modules, from its slot factories, as the runner builds it (with `common.outputs.Outputs` last if `outputs`, after `analyzers`) """
    import starsim as ss
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    ctx = SimpleNamespace(seed=seed, stop=stop, n_agents=n_agents)
    params = pars()
    kw = join(*[globals()[f'make_{s}'](params, ctx) for s in SLOTS])
    kw['analyzers'] = kw.get('analyzers', []) + list(analyzers)
    if outputs:
        from common.outputs import Outputs
        from modules import ModulesHost
        kw['analyzers'].append(Outputs(ModulesHost))
    return ss.Sim(**kw, rand_seed=seed, stop=stop, verbose=0, **kwargs)
