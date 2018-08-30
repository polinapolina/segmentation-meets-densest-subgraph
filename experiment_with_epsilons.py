import numpy as np
from utils import *
import copy
import networkx as nx
import pickle
from incremental_densest import *
from approxDP import *
#from pylab import *
from charikar import *
#import generate
import time
import uuid



def get_node_quality(graphs, generated_C):
    Q = np.empty((len(graphs), 3))
    for i in xrange(len(graphs)):
        n1, n2 = graphs[i].nodes(), generated_C[i][0].nodes()
        #print n1, n2
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
        
    
        
    
def plot_F_measure(noise_range, innerdegree_range, res):
    
    plt.figure()   
    plt.plot(innerdegree_range, res[:, 0,2], 'k:')   
    plt.plot(innerdegree_range, res[:, 1,2], 'k--')    
    plt.plot(innerdegree_range, res[:, 2,2], 'k-')
    plt.legend(['background avg. deg.: '+str(noise_range[0]),
    'background avg. deg.: '+str(noise_range[1]),
    'background avg. deg.: '+str(noise_range[2])])
    plt.xlabel('avg. degree of planted graph')
    plt.ylim(ymax=1.01)
    plt.ylabel('F-measure')
    plt.title('planted 8-subgraph, k=3')
    plt.show()
    
def plot_F_measure_intervals(noise_range, innerdegree_range, res):    
    plt.figure() 
    plt.plot(innerdegree_range, res[:,0,5], 'k:')   
    plt.plot(innerdegree_range, res[:,1,5], 'k--')    
    plt.plot(innerdegree_range, res[:,2,5], 'k-')
    plt.legend(['background avg. deg.: '+str(noise_range[0]),
    'background avg. deg.: '+str(noise_range[1]),
    'background avg. deg.: '+str(noise_range[2])])
    plt.xlabel('avg. degree of planted graph')
    plt.ylim(ymax=1.01)
    plt.ylabel('F-measure')
    plt.title('planted 8-cliques, k=3, e_dens = 0.1, e_DP = 0.1')
    plt.show()
    
def plot_PRF(noise_range, innerdegree_range, res):
    
    plt.figure()
    plt.plot(innerdegree_range, res[:, 1,0], 'r:')
    plt.plot(innerdegree_range, res[:, 1,1], 'b--')    
    plt.plot(innerdegree_range, res[:, 1,2], 'k-')
    plt.legend(['P', 'R', 'F'])
    
    plt.xlabel('avg. degree of planted subgraphs')
    plt.ylim(ymax=1.01)
    #plt.ylabel('F-measure')
    plt.title('planted 8-subgraphs, '+'avg. background degree = ' + str(noise_range[1]) +', k=3, e_dens = 0.1, e_DP = 0.1')
    plt.show()
    
def plot_PRF_intrvals(noise_range, innerdegree_range, res):
    
    plt.figure()
    plt.plot(innerdegree_range, res[:, 1,3], 'r:')
    plt.plot(innerdegree_range, res[:, 1,4], 'b--')    
    plt.plot(innerdegree_range, res[:, 1,5], 'k-')
    plt.legend(['P', 'R', 'F'])
    
    plt.xlabel('avg. degree of planted subgraphs')
    plt.ylim(ymax=1.01)
    plt.ylabel('quality of intervals')
    plt.title('planted 8-subgraphs, '+'avg. background degree = ' + str(noise_range[1]) +', k=3, e_dens = 0.1, e_DP = 0.1')
    plt.show()
    
def plot_density(noise_range, innerdegree_range, res):
    
    plt.figure()
   
    plt.plot(innerdegree_range, res[:, 1,7], 'b--')    
    plt.plot(innerdegree_range, res[:, 1,6], 'k-')
    
    plt.legend(['true total cost', 'found total cost'])
    
    plt.xlabel('avg. degree of planted subgraphs')
    #plt.ylim(ymax=1.01)
    plt.ylabel('cost')
    plt.title('planted 8-subgraphs, '+'avg. background degree = ' + str(noise_range[1]) +', k=3, e_dens = 0.1, e_DP = 0.1')
    plt.show()

    
def plot_density_compare(densest_eps_range, dp_eps_range,  res):
    
    plt.figure()
    
    
    plt.plot(densest_eps_range, res[:,9,7].astype(float)/res[:,9,6].astype(float), 'r:')
    plt.plot(densest_eps_range, res[:,2,7].astype(float)/res[:,2,6].astype(float), 'b--') 
    plt.plot(densest_eps_range, res[:,0,7].astype(float)/res[:,0,6].astype(float), 'k-')
    
    plt.legend(['dp_eps: '+str(dp_eps_range[9]),
    'dp_eps: ' + str(dp_eps_range[2]),
    'dp_eps: ' + str(dp_eps_range[0])])
    plt.xlabel('densest eps')
    #plt.ylim(ymax=1.01)
    plt.ylabel('true cost/ found cost')
    plt.title('planted 8-cliques, avg. inner density 5, avg. background density 2, , k=3')
    
def plot_density_compare_vary_dp(densest_eps_range, dp_eps_range,  res):
    
    plt.figure()
    
    
    #plt.plot(dp_eps_range, res[19,:,7]/res[19,:,6], 'g-.')
    #plt.plot(dp_eps_range, res[19,:,6])

    plt.plot(dp_eps_range, res[9,:,7].astype(float)/res[9,:,6].astype(float), 'r:')
    plt.plot(dp_eps_range, res[2,:,7].astype(float)/res[2,:,6].astype(float), 'b--') 
    plt.plot(dp_eps_range, res[0,:,7].astype(float)/res[0,:,6].astype(float), 'k-')
    
    plt.legend(['densest_eps: '+str(densest_eps_range[19]),
    'densest_eps: '+str(densest_eps_range[9]),
    'densest_eps: ' + str(densest_eps_range[2]),
    'densest_eps: ' + str(densest_eps_range[0])])
    plt.xlabel('DP epsilon')
    #plt.ylim(ymax=1.01)
    plt.ylabel('true cost/ found cost')
    plt.title('planted 8-cliques, avg. inner density = 5, avg. background density = 2, k = 3')
  
def compare_cost_heatmap():
    plt.imshow(res[:,:,7]/res[:,:,6], extent=[min(dp_eps_range), max(dp_eps_range), max(densest_eps_range), min(densest_eps_range)])
    plt.xlabel('dp_eps_range')
    plt.ylabel('densest_eps_range')
    
    plt.colorbar()
    plt.show()
    exit()  
    
def compare_F_measure_heatmap():
    plt.imshow(res[:,:,2], extent=[min(dp_eps_range), max(dp_eps_range), max(densest_eps_range), min(densest_eps_range)])
    plt.xlabel('dp_eps_range')
    plt.ylabel('densest_eps_range')
    plt.title('node F-measure')
    plt.colorbar()
    plt.show()
    exit()  
    
def compare_P_heatmap():
    plt.imshow(res[:,:,0], extent=[min(dp_eps_range), max(dp_eps_range), max(densest_eps_range), min(densest_eps_range)])
    plt.xlabel('dp_eps_range')
    plt.ylabel('densest_eps_range')
    plt.title('node Precision')
    plt.colorbar()
    plt.show()
    exit()  
    
def compare_R_heatmap():
    plt.imshow(res[:,:,1], extent=[min(dp_eps_range), max(dp_eps_range), max(densest_eps_range), min(densest_eps_range)])
    plt.xlabel('dp_eps_range')
    plt.ylabel('densest_eps_range')
    plt.title('node Recall')
    plt.colorbar()
    plt.show()
    exit()  
    
def plot_running_time():   
    
    plt.figure()
   
    plt.plot(noise_range, res[0, :, 8], 'r:')
    plt.plot(noise_range, res[1, :, 8], 'b--') 
    plt.plot(noise_range, res[2, :, 8], 'k-')
    plt.legend(['planted avg. degree: '+str(innerdegree_range[0]),
    'planted avg. degree: '+str(innerdegree_range[1]),
    'planted avg. degree: '+str(innerdegree_range[2])])
    plt.xlabel('avg. degree of backgroud')
    #plt.ylim(ymax=1.01)
    plt.ylabel('running time (sec)')
    plt.title('planted 8-cliques, k=3, e_dens = 0.1, e_DP = 0.1')
    
    




if __name__ == "__main__":
    
    alg_pars = {}
    
    densest_eps_range = np.linspace(0.01, 1., 20)
    dp_eps_range = np.linspace(0.01, 1., 20)
    res = np.empty((len(densest_eps_range), len(dp_eps_range), 1))
    
    dataset = 'students1000'
    filename = os.path.join('data', dataset + '.txt')
    TS = readdata_dict(filename, unix = False)
    k = 20
    for j in xrange(len(densest_eps_range)):        
        alg_pars['densest_eps'] = densest_eps_range[j]
        print 'epsilon for densest subgraph:', densest_eps_range[j]
        for i in xrange(len(dp_eps_range)):
            print 'epsilon for DP:', dp_eps_range[i]  
            #exp = Experiment(generator_pars, alg_pars)
            #exp.noise = noise
            alg_pars['dp_eps'] =  dp_eps_range[i]            
            st = time.time()
            ADP = ApproxDP(alg_pars['dp_eps'], k, TS)
            ADP.run_DP(alg_pars['densest_eps'])
            el_time = time.time() - st
            print 'running time:', el_time
            graphs, densities, intervals = ADP.get_sol_graphs()
            res[j, i, 0] = sum(densities)
            print 'total density:', res[j, i, 0]