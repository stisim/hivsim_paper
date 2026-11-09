"""
Optima HIV's model (the step-explicit copy models/optima/reference/model_dt.py, the identity at dt=0.2 with FIX empty)
with a tracer of people who have ever been on ART, to split Optima's in-care and lost compartments into never and
ever treated. The tracer `ever` (states × pops × time) holds, in the care and lost states, the ever-treated people in
each compartment, and in the ART states everyone in them. It is moved by the same operations as Optima's people;
nothing it does feeds back into the model. Infections by source (`raw_incibypop`) are kept, and its sexual part
separately (before mother-to-child transmission is added).

`model()` stores what's needed in TRACE; use with `op.project.model = optima_attr.model`.
"""
from pathlib import Path
import types

SRC = Path(__file__).resolve().parents[2] / 'models' / 'optima' / 'reference' / 'model_dt.py'

# (anchor: a unique piece of a line, where: 'after'/'before', code to insert at the anchor line's indentation)
PATCHES = [
    ("raw_incibypop       = zeros((nstates, npops, npts))", 'after',
     ["ever = zeros((nstates, npops, npts)) # Tracer: ever-treated people (care, lost); everyone in the ART states",
      "raw_incibypop_sex = zeros((nstates, npops, npts))"]),
    ("for t in range(startind, npts): # Loop over time", 'before',
     ["ever[alltx,:,startind] = people[alltx,:,startind]",
      "assert compareversions(version, '2.12.0') < 0 and 'agebug' not in FIX and 'linear' not in FIX, 'tracer covers the old ageing code at the calibrated step only'"]),
    ("raw_incibypop[:,:,t]        = einsum('ij,ijkl->kl',people[sus,:,t], forceinffull)/dt", 'after',
     ["raw_incibypop_sex[:,:,t] = raw_incibypop[:,:,t]"]),
    ("people[:,:,t+1] += einsum('ij,ikj->kj',people[:,:,t],thistransit[:,:,:])", 'after',
     ["ever[:,:,t+1] += einsum('ij,ikj->kj',ever[:,:,t],thistransit[:,:,:])"]),
    ("people[:, p1, t+1] -= peopleleaving # Take away from pop1...", 'before',
     ["everleaving = ever[:, p1, t+1] * thisagerate", "ever[:, p1, t+1] -= everleaving", "ever[:, p2, t+1] += everleaving"]),
    ("people[:, p1, t+1] += peoplemoving2 - peoplemoving1", 'before',
     ["evermoving1 = ever[:, p1, t+1] * thisrisktransprob",
      "evermoving2 = ever[:, p2, t+1] * thisrisktransprob * (sum(people[:, p1, t+1])/sum(people[:, p2, t+1]))",
      "ever[:, p1, t+1] += evermoving2 - evermoving1", "ever[:, p2, t+1] += evermoving1 - evermoving2"]),
    ("people[susnotonart,p,t+1] *= ratio # It's OK, so scale to match", 'after',
     ["ever[notonart,p,t+1] *= ratio"]),
    # Treatment: people taken from care onto ART, proportionally within each compartment
    ("people[care,:,t+1] -= newmovers # Shift people out of care", 'before',
     ["ever[care,:,t+1] -= newmovers*ever[care,:,t+1]/maximum(people[care,:,t+1], 1e-300)"]),
    # Linkage: people taken from dx and lost into care; the lost ones carry their ever-treated share
    ("people[lowerstate,:,t+1] -= newmovers # Shift people out of the less progressed state...", 'before',
     ["if name == 'propcare':",
      "    evermovers = newmovers[ncd4:,:]*ever[lost,:,t+1]/maximum(people[lost,:,t+1], 1e-300)",
      "    ever[lost,:,t+1] -= evermovers",
      "    ever[care,:,t+1] += evermovers"]),
    # Fewer on ART than wanted: from ART into care, all ever treated
    ("people[care,:,t+1] += newmoversusvl+newmoverssvl # Add both groups of movers into care", 'after',
     ["ever[care,:,t+1] += newmoversusvl+newmoverssvl"]),
    # Fewer in care than wanted: from care into lost
    ("people[tostate,:,t+1]    -= newmovers # Shift people out of the more progressed state...", 'before',
     ["if name == 'propcare':",
      "    evermovers = newmovers*ever[care,:,t+1]/maximum(people[care,:,t+1], 1e-300)",
      "    ever[care,:,t+1] -= evermovers",
      "    ever[lost,:,t+1] += evermovers"]),
    ("if debug: checkfornegativepeople(people, tind=t+1)", 'before',
     ["ever[alltx,:,t+1] = people[alltx,:,t+1]"]),
    ("return raw # Return raw results", 'before',
     ["TRACE.update(tvec=tvec, people=people, ever=ever, incibypop=raw_incibypop, incibypop_sex=raw_incibypop_sex, inci=raw_inci, mtct=raw_mtct, popkeys=popkeys)"]),
]


def build():
    lines = SRC.read_text().splitlines()
    for anchor, where, code in PATCHES:
        hits = [i for i, l in enumerate(lines) if anchor in l]
        assert len(hits) == 1, (anchor, hits)
        i = hits[0]
        indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
        new = [indent + c for c in code]
        lines[i:i] = new if where == 'before' else []
        if where == 'after':
            lines[i+1:i+1] = new
    mod = types.ModuleType('optima_attr_model')
    mod.__file__ = str(SRC)
    mod.TRACE = {}
    exec(compile('\n'.join(lines), str(SRC), 'exec'), mod.__dict__)
    return mod


_mod = build()
model, TRACE, FIX = _mod.model, _mod.TRACE, _mod.FIX
