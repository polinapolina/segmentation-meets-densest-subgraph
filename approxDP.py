import numpy as np
from copy import deepcopy
import networkx as nx
import pickle
from incremental_densest import *
from charikar import charikarHeap
import os.path
from optimal_max_density import *
import generate


class ApproxDP:
    def __init__(self, eps, k, timestamps):
        self.k = k
        self.eps = eps
        self.timestamps = timestamps
        self.sorted_TS = sorted(timestamps.keys())
        self.m = len(self.sorted_TS)
        self.DP = np.empty((self.k, self.m))
        self.C = np.empty((self.k, self.m))
        self.AD = {}
    
    def run_DP(self, densest_eps):
        ds = IncDensest(densest_eps)
        for i in xrange(0, self.m):
            for (n1,n2) in self.timestamps[self.sorted_TS[i]]:
                best_density =  ds.add_edge(n1, n2)
            self.DP[0, i] = best_density
            self.C[0, i] = 0
            
        for l in xrange(1, self.k):
            self.AD = {}
            for i in xrange(0, self.m):
                self.AD[i] = IncDensest(densest_eps)
                best_candidate = -1
                best_index = -1
                for idx in self.AD: 
                    for (n1,n2) in self.timestamps[self.sorted_TS[i]]:
                        best_density =  self.AD[idx].add_edge(n1, n2)
                    candidate = self.DP[l-1, idx-1] + best_density if idx > 0 else best_density
                    if candidate >= best_candidate:
                        best_candidate = candidate
                        best_index = idx                
                if i > 0 and self.DP[l, i-1] > self.DP[l-1, i] and self.DP[l, i-1] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    best_index = self.C[l, i-1]
                elif ((i > 0 and self.DP[l-1, i] > self.DP[l, i-1]) or i == 0) and self.DP[l-1, i] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    best_index = -1
                    
                self.DP[l, i] = best_candidate
                self.C[l, i] = best_index
               
                self.__SPRS(best_candidate, l)        

    def __SPRS(self, best_candidate, l):
        delta = best_candidate * self.eps / (self.k + l * self.eps)
        sorted_A = sorted(self.AD.keys()) 
        j = 0
        while j <= len(sorted_A) - 3:
            
            if self.DP[l-1, sorted_A[j+2]] - self.DP[l-1, sorted_A[j]] <= delta:
                del self.AD[sorted_A[j+1]]
                del sorted_A[j+1]
            else:
                j += 1               
        return
        
        
    def get_sol_intervals(self):
        starts = []
        end = self.m - 1
        for l in xrange(self.k-1,-1, -1):
            choice = int(self.C[l, end])
            if choice > -1:
                starts.append(choice)
                end = choice -1
            if end < 0:
                return starts[::-1]
        return starts[::-1]
        
    def get_sol_graphs(self):
        starts = self.get_sol_intervals()
        graphs = [-1]*len(starts)
        densities = [-1]*len(starts)
        intervals = [-1]*len(starts)
        for s in xrange(len(starts)):
            G = nx.Graph()
            end = len(self.sorted_TS)-1 if s == len(starts)-1 else starts[s+1]-1
            intervals[s] = (starts[s], end)
            for j in xrange(starts[s], end+1):
                for (n1, n2) in self.timestamps[self.sorted_TS[j]]:
                    G.add_edge(n1, n2)
            G = nx.Graph(G)
            G.remove_edges_from(G.selfloop_edges())
            S, best_avg = charikarHeap(G)
            graphs[s] = S
            densities[s] = best_avg
        return graphs, densities, intervals

def get_overlap(graphs):
    matches = []
    for i in xrange(len(graphs)):
        best_match = 0
        set1 = set(graphs[i].nodes())
        for j in xrange(len(graphs)):
            if j != i:
                set2 = set(graphs[j].nodes())
                H = 1.0*len(set1.intersection(set2))/len(set1.union(set2))
                if H > best_match:
                    best_match = H
        matches.append(best_match)
    return matches
        
if __name__ == "__main__":

    generator_pars = {}
    generator_pars['k'] = 2 #10
    generator_pars['B'] = 100 #100
    generator_pars['noise'] = 0.5
    generator_pars['innerdegree'] = 9 # 5.
    generator_pars['nodesInCom'] = 10 #10
    generator_pars['backgoundN'] = 100 #1000
    generator_pars['wholeSpan'] = 1000 #10000
        
    alg_pars = {}
    alg_pars['densest_eps'] = 0.1
    alg_pars['dp_eps'] = 0.1
    
    TS, backNoise, innerNoise, generated_C = generate.generate(generator_pars)
    
    ADP = ApproxDP(alg_pars['dp_eps'], generator_pars['k'], TS)
    ADP.run_DP(alg_pars['densest_eps'])
    
    graphs, densities, intervals = ADP.get_sol_graphs()
    print 'intervals: ', intervals
    print 'nodes in graphs: ', [i.nodes() for i in graphs]
    print 'densities: ', densities
    
    