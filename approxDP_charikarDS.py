import numpy as np
from copy import deepcopy
import networkx as nx
import pickle
import incremental_densest_charikar as chid
from charikar import charikarHeap
import os.path
import generate
from utils import *
import time

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
    
    def run_DP(self):
        ds = chid.IncDensest()
        for i in xrange(0, self.m):
            for (n1,n2) in self.timestamps[self.sorted_TS[i]]:
                best_density =  ds.add_edge(n1, n2)
            self.DP[0, i] = best_density
            self.C[0, i] = 0
            
        for l in xrange(1, self.k):
            self.AD = {}
            for i in xrange(0, self.m):
                self.AD[i] = chid.IncDensest()
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

    dataset = 'enron1000'
    k = 5
    dp_eps = 0.1
    
    alg_pars = {}
    alg_pars['dp_eps'] = dp_eps
    
    filename = os.path.join('data', dataset+'.txt')
    TS = readdata_dict(filename, unix = False)
    
    
    t1 = time.time()
    ADP = ApproxDP(alg_pars['dp_eps'], k, TS)
    ADP.run_DP()
    t2 = time.time()
    print 'running time:', t2 - t1
    
    graphs, densities, intervals = ADP.get_sol_graphs()
    print 'total density:', sum(densities)
    print 'intervals:', intervals
    
    