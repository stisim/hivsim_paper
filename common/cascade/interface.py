"""
What the condensed cascade needs from the model it runs in: the HIV disease and the people's life events.

A model implements `CascadeHost` (one object, bound to the sim when the cascade starts); the cascade uses nothing else.
EMOD′ implements it in models/emod/v4/cascade_host.py.

Fit with other models (not implemented):
- STIsim's `sti.HIV` with `sti.ART`: `diagnosed` and `cd4` are states of `sti.HIV`; `start_art(uids, arm)` maps to
  `ART.start_art(uids, p_effective_art=1 or 0)` (its arms are suppressive or not, by a probability), `stop_art` to
  `ART.stop_art(uids)`, `start_prep(uids, efficacy, days)` to `HIV.start_prep(uids, eff, dur, ...)`, `circumcise` to
  `HIV.circumcise(uids)`. STIsim has no WHO stage, symptom onset event, or PMTCT efficacy per mother: those would need
  adding (or WHO stage derived from CD4, symptoms from its AIDS stage).
- Optima′ is compartmental: people are counts by CD4 and cascade state, moved by rates, so a cascade that acts on
  agents (`start_art(uids)`) cannot drive it; it would need the same parameters turned into rates between its
  compartments (undx → dx → care → treated, and lost), as Optima does.
"""
from abc import ABC, abstractmethod


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
