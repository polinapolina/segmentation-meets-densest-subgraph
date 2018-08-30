import numpy as np
from utils import *
import networkx as nx
from incremental_densest import *
from approxDP import *
from pylab import *
from charikar import *
import generate
import time


def get_node_quality(graphs, generated_C):
    Q = np.empty((len(graphs), 3))
    for i in xrange(len(graphs)):
        n1, n2 = graphs[i].nodes(), generated_C[i][0].nodes()
        TP = float(len(set(n1).intersection(set(n2))))
        P = TP/len(n1)
        R = TP/len(n2)
        F = 2.0*(P * R)/(P + R) if P + R > 0 else 0
        Q[i,:] = [P,R,F]
    return Q.mean(0)

def get_interval_quality(intervals, generated_C):
    Q = np.empty((len(graphs), 3))
    for i in xrange(len(intervals)):
        s, t = intervals[i]
        s_true, t_true = generated_C[i][1]
        
        tmp = float(min(t, t_true) - max(s, s_true))
        TP = tmp if tmp >= 0 else 0.
        P = TP/(t - s) if s!=t else 0.
        R = TP/(t_true - s_true)
        F = 2.0*(P * R)/(P + R) if P + R > 0 else 0
        Q[i,:] = [P,R,F]
    return Q.mean(0)
    
def get_true_density(generated_C):
    true_density = 0
    for i in xrange(len(generated_C)):
        true_density += 2.*generated_C[i][0].number_of_edges()/generated_C[i][0].number_of_nodes()
    return true_density
     
def plot_PRF(noise_range, innerdegree_range, res):
    plt.figure()
    plt.plot(noise_range, res[5, :,0], 'r:')
    plt.plot(noise_range, res[5, :,1], 'b--')    
    plt.plot(noise_range, res[5, :,2], 'k-')
    plt.legend(['P', 'R', 'F'])
    
    plt.xlabel('avg. degree of planted graph')
    plt.title('avg. degree of background: ' + str(noise_range[5]))
    plt.ylim(ymin = 0.0, ymax=1.01)
    plt.show()
    
    plt.figure()
    plt.plot(innerdegree_range, res[:, 5,0], 'r:')
    plt.plot(innerdegree_range, res[:, 5,1], 'b--')    
    plt.plot(innerdegree_range, res[:, 5,2], 'k-')
    plt.legend(['P', 'R', 'F'])
    
    plt.xlabel('avg. degree of background')
    plt.title('avg. degree of planted graph: ' + str(noise_range[5]))
    plt.ylim(ymin = 0.0, ymax=1.01)
    plt.show()


if __name__ == "__main__":
    generator_pars = {}
    
    generator_pars['k'] = 3
    generator_pars['B'] = int(8*(7)/2.0)    
    generator_pars['innerdegree'] = 7.
    generator_pars['nodesInCom'] = 8
    generator_pars['backgoundN'] = 100
    generator_pars['wholeSpan'] = 1000
    #noise_range = np.linspace(0.5, 4.0, 50)
    noise_range = np.linspace(0.5, 4.0, 10)
    #innerdegree_range = [5]
    innerdegree_range = np.linspace(2., 7.0, 10)
    
    alg_pars = {}
    alg_pars['densest_eps'] = 0.1
    alg_pars['dp_eps'] = 0.1  
    
    res = np.empty((len(innerdegree_range), len(noise_range), 9))
    
    for j in xrange(len(innerdegree_range)):
        generator_pars['innerdegree'] = innerdegree_range[j]
        print 'inner noise (avg. degree of the planted graph):', innerdegree_range[j]
        for i in xrange(len(noise_range)):
            print 'outer noise (avg. degree of the background network):', noise_range[i]
            noise = noise_range[i]
            generator_pars['noise'] =  noise_range[i]     
            TS, backNoise, innerNoise, generated_C = generate.generate(generator_pars)
            
            st = time.time()
            ADP = ApproxDP(alg_pars['dp_eps'], generator_pars['k'], TS)
            ADP.run_DP(alg_pars['densest_eps'])
            el_time = time.time() - st
            graphs, densities, intervals = ADP.get_sol_graphs()
            nq = get_node_quality(graphs, generated_C)
            iq = get_interval_quality(intervals, generated_C)
            res[j, i, 0:3] = nq
            res[j, i, 3:6] = iq
            res[j, i,6] = sum(densities)
            print 'total density', res[j, i,6]
            print
            res[j, i,7] = get_true_density(generated_C)
            res[j, i,8] = el_time
    plot_PRF(noise_range, innerdegree_range, res)
    