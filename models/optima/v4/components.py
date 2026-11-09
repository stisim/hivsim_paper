"""
Optima′'s slot factories for the swap runner (common/swaps/README.md). Only the cascade (C) so far: Optima's own
cascade as a table (common/cascade/tables/malawi_optima_params.csv, analysis 23), run as agents by
`common.cascade.pathways.TableCascade` on the host of V's model, with the people per agent of D's model (its `SCALE`)
for the table's numbers (number on ART, viral-load monitoring).

On another model's people (analysis 23's provisional decisions): Optima's groups by population come from the host's
`population()` (poplabel); care and lost live only in the cascade's record where the host keeps none (EMOD′, HIVsim),
and `switch_arm` changes nothing there; the PMTCT, PrEP and PEP coverages are not applied (no per-person form).
Optima′'s other slots (D, RN, V) aren't built for the runner: Optima′ runs from `optima_prime.make_sim` (ODE) or
`agents.make_sim` (agents).

The runner puts this folder on sys.path: models/emod/v4 has a module of the same name as this folder's
`extract_cascade`, so code of another model in the same hybrid should not import it flat after this file is loaded.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TABLE = REPO / 'common' / 'cascade' / 'tables' / 'malawi_optima_params.csv'


def __getattr__(name):
    if name == 'HOST': # OptimaHost, when Optima′'s V is in the hybrid (imported only when asked for)
        from host import OptimaHost
        return OptimaHost
    raise AttributeError(name)


def pars(pid=None):
    """ Optima's Malawi project has one parameter set: its cascade table """
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from common.cascade.params import Params
    return Params(TABLE)


def make_C(params, ctx):
    """ Optima's cascade as a table, on the host of V's model, after what that model's `for_cascade()` gives """
    from common.cascade.pathways import TableCascade
    vm = ctx.models[ctx.choice['V']]
    pre = vm.for_cascade(ctx.pars[ctx.choice['V']], ctx) if hasattr(vm, 'for_cascade') else {}
    scale = getattr(ctx.models[ctx.choice['D']], 'SCALE', None) if 'D' in ctx.choice else None
    intvs = list(pre.get('interventions', [])) + [TableCascade(params, vm.HOST, scale=scale)]
    return dict(pre, interventions=intvs)
