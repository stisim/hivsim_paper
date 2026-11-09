"""
EMOD′'s slot factories as harness.make_sim(pid, seed, cascade=True, condensed=True, record=False) builds it, for
testing the runner until models/emod/v4/components.py exists. C's cascade runs on the host of V's model.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'models' / 'emod' / 'v4'))
import harness # noqa: E402,F401 (first, so `cascade` is EMOD′'s)
import inputs # noqa: E402
from demography import EmodDemography # noqa: E402
from network import EmodNetwork # noqa: E402
from risk import EmodRisk # noqa: E402
from hiv import EmodHIV # noqa: E402
from cascade_host import EmodHost # noqa: E402
from extract_cascade import extract # noqa: E402
from common.cascade.cascade import condensed_cascade # noqa: E402
from common.cascade.params import Params # noqa: E402

HOST = EmodHost
SCALE = 100 # People per agent (1% scale)


def pars(pid):
    return inputs.load(pid)


def make_D(inp, ctx):
    return dict(n_agents=int(inp.init_pop.sum()), start=inp.start, dt=1/12, aging=False, demographics=[EmodDemography(inp)])


def make_RN(inp, ctx):
    return dict(networks=[EmodNetwork(inp, record=False)], interventions=[EmodRisk(inp)])


def make_V(inp, ctx):
    return dict(diseases=[EmodHIV(inp, record=False)])


def make_C(inp, ctx):
    host = ctx.models[ctx.choice['V']].HOST
    return dict(interventions=condensed_cascade(Params(*extract(inp)), host), analyzers=[harness.Snapshots(), harness.EmodReport()])
