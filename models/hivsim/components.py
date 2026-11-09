"""
HIVsim's slot factories (docs/swap_design.md §1.1, in the runner's interface, common/swaps/README.md): `make_<slot>(pars, ctx)`
gives that slot's part of the sim as `sti.Sim` keywords (lists of modules, and the sim's own parameters). `make_sim()`
puts them together in the canonical order; `malawi.make_sim()` and `swap.make_sim()` are built on it, and the runner's
`build()` calls the same factories (with `SIM`, sti.Sim, as the sim, since D gives a location for STIsim to build).

    pars: the run's parameters (`pars(pid)`), as in swap.BEST: `hiv_*` for sti.HIV, `nw_*` for sti.StructuredSexual,
          and test_scale, test_rel_m, test_start (testing), condom_scale (condoms)
    ctx:  from the runner, `choice`, `models`, `pars`, `seed`, `stop`, `cell`; from make_sim(), also: n_agents (default
          10,000); `cascade`: 'own' (HIVsim's testing, ART and VMMC; the default) or 'condensed' (the condensed cascade,
          with `params`, by default EMOD's Malawi table)

Fine slots, and the runner's coarse ones made of them: RN = R + N, V = T + H + S. Where one module covers several
slots, the first slot it belongs to builds it, and the others give nothing:
- R (risk groups, sex work) is inside N's `sti.StructuredSexual` (`risk_group`, `fsw`, `client`).
- T (transmission) and S (seeding) are inside H's `sti.HIV` (`beta_m2f`, `rel_trans`/`rel_sus`; `init_prev_data`).
- ART is in C (`make_art()`), as the cascade's treatment, though its biology is H's (swap_design §1.1).
- M (measurement) is sti.HIV's own results, and V's `ss.infection_log()` (for common/outputs.py).
"""
import sys
from pathlib import Path
from types import SimpleNamespace
import sciris as sc
import starsim as ss
import stisim as sti

here = Path(__file__).resolve().parent
for p in (here / 'malawi', here, here.parents[1]):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import malawi # noqa: E402
from stisim_host import StisimHost, art_for_cascade # noqa: E402

HOST = StisimHost # The CascadeHost, when HIVsim's V is in the hybrid
SIM = sti.Sim # The sim the runner makes when HIVsim's D is in the hybrid
SLOTS = ['D', 'RN', 'V', 'C', 'M'] # The order make_sim() calls them in (the runner's canonical order)
LISTS = ['demographics', 'networks', 'diseases', 'interventions', 'connectors', 'analyzers'] # The keywords that hold lists of modules


def pars(pid):
    """ Parameter set `pid`: for now every pid is hivsim_malawi experiment 03's best run (swap.BEST), until models/hivsim/reference/ has the sets """
    import swap
    return dict(swap.BEST)


def opt(ctx, key, default):
    """ A make_sim() option; the runner gives none, so the defaults are HIVsim Malawi with its own care """
    return getattr(ctx, key, default)


def sub(pars, prefix):
    return {k[len(prefix):]: v for k, v in pars.items() if k.startswith(prefix)}


def make_D(pars, ctx):
    """ Demography: STIsim's Malawi demographics (pregnancy, deaths, the 1985 age structure), from 1985 """
    malawi.check_data()
    return dict(n_agents=opt(ctx, 'n_agents', 10e3), start=malawi.start, datafolder=malawi.datadir, location=malawi.location)


def make_R(pars, ctx):
    """ Risk and sex work: inside sti.StructuredSexual (make_N) """
    return {}


def make_N(pars, ctx):
    """ Sexual contact: sti.StructuredSexual with Malawi's condom use (times condom_scale, at most 1), and the prenatal network (mother to fetus) """
    condoms = malawi.read('condom_use.csv')
    if condoms is not None:
        condoms = condoms.set_index('partnership').mul(pars.get('condom_scale', 1)).clip(upper=1).reset_index()
    return dict(networks=[sti.StructuredSexual(condom_data=condoms, **sub(pars, 'nw_')), ss.PrenatalNet()])


def make_T(pars, ctx):
    """ Transmission: inside sti.HIV (make_H) """
    return {}


def make_H(pars, ctx):
    """ Natural history (and, as sti.HIV, transmission and seeding): beta_m2f, the 1985 prevalence by group, the surveys' age groups """
    hiv = sti.HIV(**{
        'beta_m2f': 0.01, # To calibrate
        'init_prev_data': malawi.read('init_prev_hiv.csv'),
        'age_bins': [0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100], # The surveys' age groups
    } | sub(pars, 'hiv_'))
    return dict(diseases=[hiv])


def make_S(pars, ctx):
    """ Seeding: inside sti.HIV (make_H, `init_prev_data`) """
    return {}


def make_art(kind):
    """ HIVsim's ART, the one place it is made: `own`, to the numbers on ART and the share suppressed; `cascade`, holding ART's states for the condensed cascade """
    if kind == 'own':
        return sti.ART(coverage=malawi.read('n_art.csv', index_col='year'), vls_coverage=malawi.read('vls.csv', index_col='year'))
    return art_for_cascade()


def make_C(pars, ctx):
    """ Care cascade and prevention: HIVsim's testing, ART and VMMC (`own`), or the condensed cascade (`condensed`) """
    if opt(ctx, 'cascade', 'own') == 'own':
        intvs = malawi.make_testing(pars.get('test_scale', 1), pars.get('test_rel_m', 1), pars.get('test_start', 2004)) + [make_art('own')]
        vmmc = malawi.read('vmmc.csv') # Share of men circumcised, by age
        if vmmc is not None:
            intvs.append(sti.VMMC(coverage=vmmc))
        return dict(interventions=intvs)
    import swap
    from common.cascade.cascade import condensed_cascade
    params = opt(ctx, 'params', None) or swap.emod_table()
    return dict(interventions=[make_art('cascade'), swap.traditional_vmmc()] + condensed_cascade(params, StisimHost))


def make_M(pars, ctx):
    """ Measurement: sti.HIV's results (the infection log is V's); any analyzers given """
    return dict(analyzers=sc.tolist(opt(ctx, 'analyzers', None)))


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


def make_RN(pars, ctx):
    """ Risk and sexual contact (the runner's coarse slot) """
    return join(make_R(pars, ctx), make_N(pars, ctx))


def make_V(pars, ctx):
    """ The virus: transmission, natural history and seeding (the runner's coarse slot), with Starsim's infection log (sources of transmission, for common/outputs.py; it doesn't change the run) """
    return join(make_T(pars, ctx), make_H(pars, ctx), make_S(pars, ctx), dict(analyzers=[ss.infection_log()]))


def make_sim(pars=None, ctx=None, sim_kwargs=None, **kwargs):
    """ HIVsim Malawi from its slot factories; `ctx` keys as in the module docstring, or as keywords; `sim_kwargs` go to sti.Sim (analyzers among them after the slots') """
    ctx = SimpleNamespace(**(dict(seed=1, stop=2021, verbose=1/12) | dict(ctx or {}) | kwargs))
    extra = dict(sim_kwargs or {})
    ctx.analyzers = extra.pop('analyzers', None)
    kw = join(*[globals()[f'make_{s}'](dict(pars or {}), ctx) for s in SLOTS])
    for k in LISTS: # Other modules given by the caller come after the slots'
        if k in extra:
            kw[k] = kw.get(k, []) + sc.tolist(extra.pop(k))
    return SIM(**kw, stop=ctx.stop, rand_seed=ctx.seed, verbose=ctx.verbose, **extra)
