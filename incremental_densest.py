import numpy as np
import math
import fibonacci_heap_mod as fb
import copy
import networkx as nx
import pickle

class IncDensest:
    def __init__(self, eps):
        self.G = nx.Graph()
        self.node2set = {}
        self.delta = {}
        self.beta = 0.5
        self.eps = eps
        self.best_density = 0.5
        self.bestiter = -1
        
    def __add_new_node(self, n1):
        self.G.add_node(n1)    
        self.node2set[n1] = 0
        self.delta[n1] = {}
        return 
     
    def add_edge(self, n1, n2):
        if not self.G.has_node(n1):
            self.__add_new_node(n1)
        if not self.G.has_node(n2):
            self.__add_new_node(n2)           
        if not self.G.has_edge(n1, n2):
            self.G.add_edge(n1,n2)
            rebuild = self.__add(n1, n2)
            if rebuild:
                self.beta, density = self.__find_densest(self.best_density)  
                if density > self.best_density:
                    self.best_density = density
                self.__get_degrees()
        return self.best_density
    
    def __add(self, u, v):
        rebuild = False
        if self.node2set[u] >= self.node2set[v]:
            if self.node2set[u] not in self.delta[v]:
                self.delta[v][self.node2set[u]] = 0.
            self.delta[v][self.node2set[u]] += 1.
        if self.node2set[u] <= self.node2set[v]:
            if self.node2set[v] not in self.delta[u]:
                self.delta[u][self.node2set[v]] = 0.
            self.delta[u][self.node2set[v]] += 1.

        S = [u, v]
        maxt = int(np.ceil(math.log(float(self.G.number_of_nodes()), self.eps + 1.)))
        while S:
            w = S.pop()
            idx = self.node2set[w]
            s = sum([self.delta[w][k] for k in self.delta[w].keys()])
            if s < 2.*(1. + self.eps)*self.beta:
                pass
            else:
                for key in sorted(self.delta[w].keys()):
                    s -= self.delta[w][key]
                    if s < 2.*(1. + self.eps)*self.beta:
                        new_idx = key + 1
                        break          
                if new_idx == maxt:
                    return True
                else:
                    for u in self.G.neighbors_iter(w):
                        S.append(u)
                        if self.node2set[w] >= self.node2set[u]:
                            self.delta[u][idx] -= 1.
                            if new_idx not in self.delta[u]:
                                self.delta[u][new_idx] = 0.
                            self.delta[u][new_idx] += 1.
                            
                        elif self.node2set[w] < self.node2set[u]:
                            u_idx = self.node2set[u]
                            if u_idx <= new_idx:                            
                                if new_idx not in self.delta[u]:
                                    self.delta[u][new_idx] = 0.
                                self.delta[u][new_idx] += 1.
                            if u_idx < new_idx:
                                self.delta[w][u_idx] -= 1.                    
                    self.delta[w] = {k: self.delta[w][k] for k in self.delta[w] if k >= new_idx}
                                
                    self.node2set[w] = new_idx
        return False
                

    def __find_densest(self, rho):
        local_beta = max(0.25/(1 + self.eps), rho*(1. + self.eps))
        
        out_density = rho
        while True:
            best_avg = self.__find(local_beta)            
            if best_avg >= local_beta:
                local_beta = (1. + self.eps) * best_avg
                out_density = best_avg
            else:
                self.__find(local_beta)
                return local_beta, out_density
        return 
        
    def __find(self, beta):             
        E = self.G.number_of_edges()
        N = self.G.number_of_nodes()
        number_of_nodes = self.G.number_of_nodes()
        fib_heap = fb.Fibonacci_heap()
        entries = {}
        self.node2set, set2nodes = {}, {}
        S = copy.deepcopy(self.G)
        
        for node, deg in S.degree_iter():
            entries[node] = fib_heap.enqueue(node, deg)
            
        best_avg = 0.0    
        iter = 0
        avg_degree = 1.0*S.number_of_edges()/S.number_of_nodes()    
        if best_avg <= avg_degree:
            best_avg = avg_degree
            best_iter = iter
        
        maxiter = np.ceil(math.log(float(number_of_nodes), self.eps+1))
        while fib_heap and iter < maxiter and fib_heap.min().get_priority() < 2. * (1 + self.eps) * beta:
            set2nodes[iter] = []
            while fib_heap and fib_heap.min().get_priority() < 2. * (1 + self.eps) * beta:
                min_deg_obj = fib_heap.dequeue_min()
                min_deg_node = min_deg_obj.get_value()
                
                set2nodes[iter].append(min_deg_node)
                self.node2set[min_deg_node] = iter
            for n in set2nodes[iter]:            
                for neigh in S.neighbors_iter(n):
                    if neigh not in set2nodes[iter]:
                        fib_heap.decrease_key(entries[neigh], 1)                    
            S.remove_nodes_from(set2nodes[iter])            
            iter += 1
            
            avg_degree = 1.0*S.number_of_edges()/S.number_of_nodes() if S else 0.0 
            
            if best_avg <= avg_degree:
                best_avg = avg_degree
                best_iter = iter
                
        i = min(iter, maxiter)
        while fib_heap:            
            min_deg_node = fib_heap.dequeue_min().get_value()
            self.node2set[min_deg_node] = i
        self.bestiter = best_iter
        return best_avg
    
    def __get_degrees(self):
        self.delta = {}
        for node in self.G.nodes_iter():
            self.delta[node] = {}
            for neigh in self.G.neighbors_iter(node):
                if self.node2set[neigh] >= self.node2set[node]:
                    if self.node2set[neigh] not in self.delta[node]:
                        self.delta[node][self.node2set[neigh]] = 0.
                    self.delta[node][self.node2set[neigh]] += 1.
        return
    
def get_best_graph(G, set2nodes, best_iter):
    S = copy.deepcopy(G)       
    for i in xrange(best_iter):
        S.remove_nodes_from(set2nodes[i])
    return S 

def check_invariant(G, node2set, eps, beta):
    S = nx.Graph()
    correct = True
    maxiter = np.ceil(math.log(float(G.number_of_nodes()), eps+1))
    max_set = int(max(node2set.values())) if node2set else 0
    if not max_set <= maxiter-1:
        print max_set, maxiter
        print 'incorrect range of iterations'
        exit()
    for i in xrange(0, max_set+1):
        nodes = [k for k,v in node2set.iteritems() if v >= i]
        S = G.subgraph(nodes)
        tmp = set([n for n, d in S.degree_iter() if d < 2.*beta*(1.+eps)])
        if tmp != set([k for k,v in node2set.iteritems() if v == i]):
            print 'wrong iteration:', i
            print 'should be:', tmp
            print 'mine:',  nodes, node2set
            correct = False
            
    return correct
    
def check_delta(G, delta):
    correct = True
    true_delta = get_degrees(G)
    for u in true_delta:
        for i in true_delta[u]:
            if (i not in delta[u] and true_delta[u][i]!=0) or true_delta[u][i] != delta[u][i]:
                print u, i, delta, true_delta
                correct = False
                print 'delta is incorrect'
                exit()
    for u in delta:
        for i in delta[u]:
            if (i in true_delta[u] or delta[u][i]!=0) and true_delta[u][i] != delta[u][i]:
                print u, i
                correct = False
                print 'delta is incorrect'
                exit()
    return correct


if __name__ == "__main__":
    eps = 0.4
    nodes = 100
    beta = 0.5

    H = nx.erdos_renyi_graph(nodes, 0.7)
    iterations = 10000
    timestamps = utils.generateTS(H, iterations)
    
    ds = IncDensest(eps)
    
    for (tst, n1, n2) in timestamps:
        best_density =  ds.add_edge(n1, n2)
        print best_density
 