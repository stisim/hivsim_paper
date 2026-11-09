""" The other toy model: more contacts, lower transmission """
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from toy_a import make_D, make_N, make_V, BETA # noqa: E402,F401

N_CONTACTS = 8


def pars(pid):
    return dict(pid=pid, n_contacts=N_CONTACTS, beta=BETA/2)
