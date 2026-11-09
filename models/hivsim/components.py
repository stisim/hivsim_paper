"""
HIVsim's slot factories (docs/swap_design.md §1.1): `make_<slot>(pars, context)` gives that slot's part of the sim, as
keywords of `sti.Sim` (lists of modules, and the sim's own parameters); `make_sim()` puts them together in the canonical
order. `malawi.make_sim()` and `swap.make_sim()` are built on these.

    pars:    the run's parameters, as in swap.BEST: `hiv_*` for sti.HIV, `nw_*` for sti.StructuredSexual, and
             test_scale, test_rel_m, test_start (testing), condom_scale (condoms)
    context: seed, n_agents, stop, verbose; `cascade`: 'own' (HIVsim's testing, ART and VMMC) or 'condensed' (the
             condensed cascade, with `params`, by default EMOD's Malawi table); `sim_kwargs`, passed to sti.Sim

Where one module covers several slots, the first slot it belongs to builds it, and the others give nothing:
- R (risk groups, sex work) is inside N's `sti.StructuredSexual` (`risk_group`, `fsw`, `client`).
- T (transmission) and S (seeding) are inside H's `sti.HIV` (`beta_m2f`, `rel_trans`/`rel_sus`; `init_prev_data`).
- ART is in C (`make_art()`), as the cascade's treatment, though its biology is H's (swap_design §1.1).
- M (measurement) is sti.HIV's own results; no analyzers.
"""
import sys
from pathlib import Path
import sciris as sc
import starsim as ss
import stisim as sti

here = Path(__file__).resolve().parent
for p in (here / 'malawi', here, here.parents[1]):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import malawi # noqa: E402

SLOTS = ['D', 'R', 'N', 'T', 'H', 'S', 'C', 'M']
KINDS = ['demographics', 'networks', 'diseases', 'interventions', 'analyzers'] # The keywords that hold lists of modules


def sub(pars, prefix):
    return {k[len(prefix):]: v for k, v in pars.items() if k.startswith(prefix)}


def make_D(pars, context):
    """ Demography: STIsim's Malawi demographics (pregnancy, deaths, the 1985 age structure), from 1985 """
    malawi.check_data()
    return dict(n_agents=context.n_agents, start=malawi.start, datafolder=malawi.datadir, demographics=malawi.location)


def make_R(pars, context):
    """ Risk and sex work: inside sti.StructuredSexual (make_N) """
    return {}


def make_N(pars, context):
    """ Sexual contact: sti.StructuredSexual with Malawi's condom use (times condom_scale, at most 1), and the prenatal network (mother to fetus) """
    condoms = malawi.read('condom_use.csv')
    if condoms is not None:
        condoms = condoms.set_index('partnership').mul(pars.get('condom_scale', 1)).clip(upper=1).reset_index()
    return dict(networks=[sti.StructuredSexual(condom_data=condoms, **sub(pars, 'nw_')), ss.PrenatalNet()])


def make_T(pars, context):
    """ Transmission: inside sti.HIV (make_H) """
    return {}


def make_H(pars, context):
    """ Natural history (and, as sti.HIV, transmission and seeding): beta_m2f, the 1985 prevalence by group, the surveys' age groups """
    hiv = sti.HIV(**{
        'beta_m2f': 0.01, # To calibrate
        'init_prev_data': malawi.read('init_prev_hiv.csv'),
        'age_bins': [0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100], # The surveys' age groups
    } | sub(pars, 'hiv_'))
    return dict(diseases=hiv)


def make_S(pars, context):
    """ Seeding: inside sti.HIV (make_H, `init_prev_data`) """
    return {}


def make_art(kind):
    """ HIVsim's ART, the one place it is made: `own`, to the numbers on ART and the share suppressed; `cascade`, holding ART's states for the condensed cascade """
    if kind == 'own':
        return sti.ART(coverage=malawi.read('n_art.csv', index_col='year'), vls_coverage=malawi.read('vls.csv', index_col='year'))
    from cascade_host import art_for_cascade
    return art_for_cascade()


def make_C(pars, context):
    """ Care cascade and prevention: HIVsim's testing, ART and VMMC (`own`), or the condensed cascade (`condensed`) """
    if context.get('cascade', 'own') == 'own':
        intvs = malawi.make_testing(pars.get('test_scale', 1), pars.get('test_rel_m', 1), pars.get('test_start', 2004)) + [make_art('own')]
        vmmc = malawi.read('vmmc.csv') # Share of men circumcised, by age
        if vmmc is not None:
            intvs.append(sti.VMMC(coverage=vmmc))
        return dict(interventions=intvs)
    import swap
    from cascade_host import StisimHost
    from common.cascade.cascade import condensed_cascade
    params = context.get('params') or swap.emod_table()
    return dict(interventions=[make_art('cascade'), swap.traditional_vmmc()] + condensed_cascade(params, StisimHost))


def make_M(pars, context):
    """ Measurement: sti.HIV's results; no analyzers """
    return {}


FACTORIES = dict(D=make_D, R=make_R, N=make_N, T=make_T, H=make_H, S=make_S, C=make_C, M=make_M)


def assemble(parts):
    """ The slots' parts as sti.Sim's keywords: module lists joined in slot order, other keywords given once """
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


def make_sim(pars=None, context=None, **kwargs):
    """ HIVsim Malawi from its slot factories; `context` keys as in the module docstring, or as keywords """
    pars = dict(pars or {})
    context = sc.objdict(dict(seed=1, n_agents=10e3, stop=2021, verbose=1/12, cascade='own', params=None, sim_kwargs={}) | dict(context or {}) | kwargs)
    kw = assemble(FACTORIES[s](pars, context) for s in SLOTS)
    extra = dict(context.sim_kwargs)
    for k in KINDS: # Modules given by the caller (e.g. analyzers) come after the slots'
        if k in extra:
            kw[k] = kw.get(k, []) + sc.tolist(extra.pop(k))
    return sti.Sim(**kw, stop=context.stop, rand_seed=context.seed, verbose=context.verbose, **extra)
