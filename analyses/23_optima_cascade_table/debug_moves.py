""" Moves by cascade state, table (TableCascade's moves, agents) or own (transition results, people), to 1994: python debug_moves.py table|own SEED """
import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import run as R, numpy as np, sciris as sc, host as H
from common.cascade.pathways import TableCascade
from common.cascade.params import Params
mode, seed = sys.argv[1], int(sys.argv[2])
class TC(TableCascade):
    def act(self, uids, src, dst):
        self.log.setdefault((src,dst), np.zeros(2000))[self.ti] += len(uids)
        super().act(uids, src, dst)
P = Params(R.TABLE)
kw = dict(n_agents=100_000, rand_seed=seed, end=1994)
if mode == 'table':
    tc = TC(P, H.OptimaHost); tc.log = {}
    sim = R.ag.make_sim(cascade=P, interventions=[tc, H.OptimaTracker()], **kw).run()
    flows = {k: v[:len(sim.t)] for k, v in sim.interventions.table_cascade.log.items()}
else:
    sim = R.ag.make_sim(**kw).run()
    r = sim.results.optima
    flows = {}
    for k in r.keys():
        if '_to_' in k:
            a, b = k.split('_to_'); sa, sb = a.split('_')[0], b.split('_')[0]
            if sa != sb: flows[(sa, sb)] = flows.get((sa, sb), 0) + np.asarray(r[k])
names=[s.replace('-','_') for s in R.opp.load().settings['statelabels']]
cnt = {c: sum(np.asarray(sim.results.optima[n]) for n in names if n.startswith(c+'_')) for c in ['undx','dx','care','lost','usvl','svl']}
sc.save(f'' + str(R.OUT) + '/debug_{mode}_{seed}.obj', dict(flows=flows, cnt=cnt))
