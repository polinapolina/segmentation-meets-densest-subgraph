import numpy as np
import utils
from copy import deepcopy
import networkx as nx
import pickle
from incremental_densest_overlap import *
from kcharikar import kcharikarHeap
import matplotlib.pyplot as plt
import os.path
import madoka

class ApproxDPoverlap:
    def __init__(self, eps, k, gamma, timestamps):
        self.k = k
        self.gamma = gamma
        self.eps = eps
        self.timestamps = timestamps
        self.sorted_TS = sorted(timestamps.keys())
        self.m = len(self.sorted_TS)
        self.DP = np.empty((self.k, self.m))
        tmp = [madoka.Sketch()] * self.m
        self.sketch = [tmp[:] for i in range(self.k)]
        
        self.C = np.empty((self.k, self.m))
        self.AD = {}
        
    def __add_to_sketch(self, best_set, l, i):
        for e in best_set:            
            self.sketch[l][i].inc(str(e))
        
    
    def run_DP(self, densest_eps):
        ds = IncDensestOverlap(densest_eps, madoka.Sketch(), self.gamma)
        for i in xrange(0, self.m):
            for (n1,n2) in self.timestamps[self.sorted_TS[i]]:
                best_density, best_set = ds.add_edge(n1, n2)
            self.DP[0, i] = best_density
            self.C[0, i] = 0
            self.__add_to_sketch(best_set, 0, i)
            
            
        for l in xrange(1, self.k):
            self.AD = {}
            for i in xrange(0, self.m):
                self.AD[i] = IncDensestOverlap(densest_eps, self.sketch[l-1][i-1] if i > 0 else madoka.Sketch(), self.gamma)
                best_candidate = -1
                best_index = -1
                for idx in self.AD: 
                    for (n1,n2) in self.timestamps[self.sorted_TS[i]]:
                        best_density, best_set =  self.AD[idx].add_edge(n1, n2)
                    candidate = self.DP[l-1, idx-1] + best_density if idx > 0 else best_density
                    if candidate > best_candidate:
                        best_candidate = candidate
                        self.sketch[l][i].clear()
                        self.__add_to_sketch(best_set, l, i)
                        self.C[l, i] = idx
                self.DP[l, i] = best_candidate
                self.sketch[l][i].merge(self.sketch[l-1][idx-1])
                
                if i > 0 and self.DP[l, i-1] > self.DP[l-1, i] and self.DP[l, i-1] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    self.sketch[l][i] = self.sketch[l][i-1] #ok to have shallow copy, we never alter filled entries of DP table
                    best_index = self.C[l, i-1]
                elif ((i > 0 and self.DP[l-1, i] > self.DP[l, i-1]) or i == 0) and self.DP[l-1, i] > best_candidate:
                    best_candidate = self.DP[l, i-1]
                    best_index = -1
                    self.sketch[l][i] = self.sketch[l][i-1]
               
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
        selected_nodes = set()
        for s in xrange(len(starts)):
            G = nx.Graph()
            end = len(self.sorted_TS)-1 if s == len(starts)-1 else starts[s+1]-1
            intervals[s] = (starts[s], end)
            for j in xrange(starts[s], end+1):
                for (n1, n2) in self.timestamps[self.sorted_TS[j]]:
                    G.add_edge(n1, n2)
            G = nx.Graph(G)
            G.remove_edges_from(G.selfloop_edges())
            node_weights = {}
            for n in G.nodes_iter():
                if n not in selected_nodes:
                    node_weights[n] = 1.0
                else:
                    node_weights[n] = 0.0
            S, best_avg = kcharikarHeap(G, node_weights, self.gamma)
            graphs[s] = S
            selected_nodes.update(set(S.nodes()))
            densities[s] = 2.* S.number_of_edges()/S.number_of_nodes()
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
    densest_eps = 0.1
    dp_eps = 0.1
    k = 10
    filename = os.path.join('data', 'enron.txt')
    timestamps = utils.readdata_dict(filename, unix = False)
    
    ADP = ApproxDP(dp_eps, k, timestamps)
    ADP.run_DP(densest_eps)
    graphs, densities = ADP.get_sol_graphs()
    matches = get_overlap(graphs)
    print 'matches', np.mean(matches), np.min(matches), np.max(matches)
    print 'densities', np.mean(densities), np.min(densities), np.max(densities)
    