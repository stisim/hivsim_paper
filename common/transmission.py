"""
Transmissions by the source's cascade stage (docs/swap_design.md 2.2), after Daniel Citron's
models/emod/reference/source/Analysis/mwi_transmission.Rmd, as one partition of the transmissions.

Any model gives a log with one row per transmission: `year`, `src`, `dst`, and the source's state then:
`diagnosed`, `on_art`, `suppressive` (on the suppressive ART arm), `n_starts` (ART starts so far), `since_start`
(years since the last ART start) and `since_stop` (years since the last stop; NaN if never). Events in the same time
step as the transmission count as before it (EMOD's TransmissionReport is mixed on this: for 0.5% of transmissions the
source's ARTStatus there is from before an ART start or stop in that step, and as often from after). `classify()` turns
the state into categories, `shares()` into shares by year.

- `from_emod(folder)`: the log from EMOD's TransmissionReport and ReportEventRecorder (run_emod.py --transmission).
- `TransmissionLog`: a Starsim analyzer that writes the log for an agent-based model with a `CascadeHost` that has
  `cascade_stage()`; it keeps the ART history itself by watching who is on ART each step.
"""
import numpy as np
import pandas as pd

HALF = 0.5 - 1e-6 # Years: the 6-month windows (six monthly steps, 182.5 days, count as 6 months in both EMOD and EMOD′)
CATEGORIES = ['undiagnosed', 'diagnosed_naive', 'art_first_lt6m', 'art_reinit_lt6m', 'art_supp', 'art_nonsupp',
              'off_first', 'off_later_lt6m', 'off_later_ge6m']
COARSE = dict(undiagnosed='undiagnosed', diagnosed_naive='diagnosed_naive', art_first_lt6m='art_lt6m', art_reinit_lt6m='art_lt6m',
              art_supp='art_supp', art_nonsupp='art_nonsupp', off_first='interrupted', off_later_lt6m='interrupted', off_later_ge6m='interrupted')
COARSE_ORDER = ['undiagnosed', 'diagnosed_naive', 'art_lt6m', 'art_supp', 'art_nonsupp', 'interrupted']


def classify(log):
    """ The category of each transmission (a Series of CATEGORIES), from the source's state """
    d = log
    on, n = d.on_art.to_numpy(bool), d.n_starts.to_numpy()
    assert not (on & (n == 0)).any(), 'on ART with no ART start'
    lt6 = d.since_start.to_numpy() < HALF
    cat = np.select([(n == 0) & ~d.diagnosed.to_numpy(bool), n == 0,
                     on & lt6 & (n == 1), on & lt6, on & d.suppressive.to_numpy(bool), on,
                     n == 1, d.since_stop.to_numpy() < HALF],
                    CATEGORIES[:-1], CATEGORIES[-1])
    return pd.Series(pd.Categorical(cat, CATEGORIES), index=d.index, name='category')


def shares(log, coarse=False, by=()):
    """ Share of each year's transmissions by category (columns), with `n` the number of transmissions """
    cat = classify(log)
    if coarse:
        cat = pd.Series(pd.Categorical(cat.map(COARSE), COARSE_ORDER), index=cat.index)
    keys = [np.floor(log.year).astype(int).rename('year')] + [log[k] for k in by]
    counts = pd.crosstab(keys, cat).reindex(columns=cat.cat.categories, fill_value=0)
    n = counts.sum(1)
    return counts.div(n, axis=0).assign(n=n)


def from_emod(folder):
    """ The log from one EMOD run's output folder (TransmissionReport.csv, ReportEventRecorder.csv with HIVNewlyDiagnosed,
    StartedART and StoppedART from the start). On ART is from the start and stop events; the arm is the source's
    ARTStatus at the transmission (`art_status`, for checking) """
    tr = pd.read_csv(f'{folder}/TransmissionReport.csv', usecols=['SIM_TIME', 'YEAR', 'SRC_ID', 'DEST_ID', 'SRC_GENDER', 'SRC_AGE', 'SRC_IP'])
    art = tr.SRC_IP.str.extract(r'ARTStatus:(\w+)')[0]
    log = pd.DataFrame(dict(time=tr.SIM_TIME, year=tr.YEAR, src=tr.SRC_ID, dst=tr.DEST_ID, src_female=tr.SRC_GENDER == 1, src_age=tr.SRC_AGE/365,
                            art_status=art, suppressive=art == 'OnEffectiveART'))
    ev = pd.read_csv(f'{folder}/ReportEventRecorder.csv', usecols=['Time', 'Event_Name', 'Individual_ID'])
    log = log.sort_values('time').reset_index(drop=True)
    diag = ev[ev.Event_Name == 'HIVNewlyDiagnosed'].groupby('Individual_ID').Time.min()
    for name, what in (('StartedART', 'start'), ('StoppedART', 'stop')):
        e = ev[ev.Event_Name == name].sort_values('Time')
        e = pd.DataFrame({'time': e.Time.to_numpy(float), 'src': e.Individual_ID.to_numpy(), f'n_{what}s': e.groupby('Individual_ID').cumcount().to_numpy() + 1,
                          f't_{what}': e.Time.to_numpy(float)})
        log = pd.merge_asof(log, e, on='time', by='src', direction='backward', allow_exact_matches=True) # The last one at or before
    log[['n_starts', 'n_stops']] = log[['n_starts', 'n_stops']].fillna(0).astype(int)
    log['diagnosed'] = (log.src.map(diag) <= log.time) | (log.n_starts > 0)
    log['since_start'] = (log.time - log.t_start)/365
    log['since_stop'] = (log.time - log.t_stop)/365
    log['on_art'] = log.n_starts > log.n_stops
    return log.drop(columns=['t_start', 't_stop'])


try:
    import starsim as ss
except ImportError: # from_emod() and classify() don't need it
    ss = None

if ss is not None:
    class TransmissionLog(ss.Analyzer):
        """
        One row per transmission (`df`, after the run) with the source's cascade stage at that moment.

        Args:
            host (class): the model's `CascadeHost`, which needs `cascade_stage()`
            pairs (function): `pairs(sim)` gives this step's transmissions as (recipients, sources), sources -1 if none
        """
        def __init__(self, host, pairs, name='transmissionlog', **kwargs):
            super().__init__(name=name, **kwargs)
            self.host_class, self.pairs = host, pairs
            self.define_states(ss.BoolArr('on', default=False), ss.IntArr('n_starts', default=0), ss.FloatArr('t_start'), ss.FloatArr('t_stop'))
            self.rows = []
            return

        def init_post(self):
            super().init_post()
            self.host = self.host_class(self.sim)
            return

        def step(self):
            ppl, now = self.sim.people, self.t.now('year')
            au = ppl.auids
            on = self.host.cascade_stage(au)['on_art']
            on = np.where(ppl.ti_dead.raw[au] == self.ti, self.on.raw[au], on) # Those dying this step have left ART already: as before
            was = self.on.raw[au]
            self.n_starts[au[on & ~was]] += 1
            self.t_start[au[on & ~was]] = now
            self.t_stop[au[~on & was]] = now
            self.on[au] = on
            dst, src = self.pairs(self.sim)
            keep = np.asarray(src) >= 0
            dst, src = ss.uids(np.asarray(dst)[keep]), ss.uids(np.asarray(src)[keep])
            if len(src):
                st = self.host.cascade_stage(src)
                self.rows.append(pd.DataFrame(dict(year=now, src=src, dst=dst, src_female=ppl.female.raw[src], src_age=self.host.age_years(src),
                                                   diagnosed=st['diagnosed'], on_art=self.on.raw[src], suppressive=st['suppressive'],
                                                   n_starts=self.n_starts.raw[src], since_start=now - self.t_start.raw[src],
                                                   since_stop=now - self.t_stop.raw[src])))
            return

        def finalize(self):
            super().finalize()
            self.df = pd.concat(self.rows, ignore_index=True) if self.rows else pd.DataFrame()
            return
