"""
What the condensed cascade needs from the model it runs in: the HIV disease and the people's life events.

A model implements `CascadeHost` (one object, bound to the sim when the cascade starts); the cascade uses nothing else.
EMOD′ implements it in models/emod/v4/cascade_host.py.

HIVsim (STIsim's `sti.HIV` with `sti.ART`) implements it in models/hivsim/cascade_host.py, with WHO stage and first
symptoms from the CD4 count.

Optima′ is compartmental (counts by CD4 and cascade state, moved by rates), so this per-agent interface can't drive
it. Instead `rates.py` turns the same parameter table into rates between its compartments (undx → dx → care →
treated, and lost); see README.md and models/optima/v4/optima_prime.py (`make_sim(cascade=...)`).
"""
from abc import ABC, abstractmethod
import numpy as np


class CascadeHost(ABC):
    """
    The model's side of the cascade. All `uids` are arrays of agent IDs; times are in days.

    Each step, the cascade reads the events of the step (`debuts()`, `newly_symptomatic()`, `twelve_weeks_pregnant()`,
    `six_weeks_old()`), decides who is tested, staged and treated, and acts through the methods below.
    """
    def __init__(self, sim):
        self.sim = sim

    # People
    @abstractmethod
    def age_years(self, uids): ...
    @abstractmethod
    def female(self, uids): ...
    @abstractmethod
    def pregnant(self, uids): ...
    @abstractmethod
    def post_debut(self, uids): ...
    @abstractmethod
    def dying(self, uids):
        """ Whether each is dying in this step (so is given nothing) """
    def node(self, uids):
        """ Each one's node (region), from 1, for flows given by node; by default one node """
        return np.ones(len(uids), int)
    def risk(self, uids):
        """ Each one's risk group, 0-2 (low, medium, high), for flows restricted by risk """
        raise NotImplementedError(f'{type(self).__name__} has no risk groups, which a flow is restricted to')

    # Events of the step (UIDs)
    @abstractmethod
    def debuts(self): ...
    @abstractmethod
    def newly_symptomatic(self): ...
    @abstractmethod
    def twelve_weeks_pregnant(self):
        """ Pregnancies at 12 weeks (the first ANC visit) """
    @abstractmethod
    def six_weeks_old(self):
        """ Infants at 6 weeks (early infant diagnosis) """

    # HIV
    @abstractmethod
    def infected(self, uids): ...
    @abstractmethod
    def test(self, uids):
        """ An HIV test: records it (and a diagnosis if positive); returns whether each is positive """
    @abstractmethod
    def cd4(self, uids): ...
    @abstractmethod
    def record_cd4(self, uids):
        """ A CD4 count taken now; returns the lowest each has had, which staging uses """
    @abstractmethod
    def who_stage(self, uids): ...
    @abstractmethod
    def on_art(self, uids): ...
    @abstractmethod
    def start_art(self, uids, arm):
        """ Start ART in arm `arm` (e.g. 'suppressive' or 'nonsuppressive') """
    @abstractmethod
    def stop_art(self, uids): ...
    @abstractmethod
    def start_pmtct(self, uids, efficacy, days):
        """ Reduce transmission at birth by `efficacy` for `days` """
    @abstractmethod
    def circumcise(self, uids):
        """ Circumcise men (others are ignored); returns those circumcised now """
    @abstractmethod
    def start_prep(self, uids, efficacy, days): ...
