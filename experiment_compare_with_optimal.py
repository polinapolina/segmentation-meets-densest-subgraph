import sys
sys.path.append("..")
import numpy as np
import networkx as nx
from approxDP import *
from optimalDP import *
from charikar import *
from utils import *
import time
 
if __name__ == "__main__":

    alg_pars = {}
    alg_pars['densest_eps'] = 0.1
    alg_pars['dp_eps'] = 0.1
    
    #range_k = range(2,31,1)
    range_k = range(2,10,2)
    
    G = nx.fast_gnp_random_graph(50, 0.06)
    TS = generateTS_dict(G, len = 50)
    
    res = np.empty((len(range_k), 2))
    
    for i in xrange(len(range_k)):
        k = range_k[i]
        st = time.time()
        ADP = optimalDP(range_k[i], TS)
        ADP.run_DP()
        el_time = time.time() - st
        graphs, densities, intervals = ADP.get_sol_graphs()
        
        print 'Optimal:'
        print 'running time', el_time, 'total density', sum(densities)
        res[i][:] = [el_time, sum(densities)]
        
        st = time.time()
        ADP = ApproxDP(alg_pars['dp_eps'], k, TS)
        ADP.run_DP(alg_pars['densest_eps'])
        el_time = time.time() - st
        graphs, densities, intervals = ADP.get_sol_graphs()
        
        print 'Approximate (kGapprox):'
        print 'running time', el_time, 'total density', sum(densities)
    