""" A toy model for the runner's own tests: SIS on a random network, in slots D, N and V """
import starsim as ss

N_CONTACTS, BETA = 4, 0.05


def pars(pid):
    return dict(pid=pid, n_contacts=N_CONTACTS, beta=BETA*(1 + 0.1*(pid or 0)))


def make_D(p, ctx):
    return dict(n_agents=500, start=2000, dt=1/12, demographics=[ss.Births(birth_rate=20), ss.Deaths(death_rate=15)])


def make_N(p, ctx):
    return dict(networks=[ss.RandomNet(n_contacts=p['n_contacts'])])


def make_V(p, ctx):
    return dict(diseases=[ss.SIS(beta=ss.perday(p['beta']), init_prev=0.1)])
