from networkx import *
import fibonacci_heap_mod
from copy import deepcopy

def get_solution(G, g = None):
    if not g:
         g = find_densest(G)

    m, n = 1.*G.number_of_edges(), 1.*G.number_of_nodes()
    degrees = G.degree()
    nodes = G.nodes()
    
    H = G.to_directed()
    for u,v in H.edges():
        H[u][v]['cap'] = 1.
    for v in nodes:
        H.add_edge('s', v, cap = m)
        H.add_edge( v, 't', cap = m + 2.*g - 1.*degrees[v])
    cut_value, (S, T) = nx.minimum_cut(H, 's', 't', capacity='cap')
    S.remove('s')
    D = G.subgraph(S)
    return D, 2.*D.number_of_edges()/D.number_of_nodes()
    
    
def find_densest(G, g = None):
    
    ub = 0.5*(G.number_of_nodes()-1.)
    lb = 0.
    if not g:     
        g = (lb + ub)/2.
    m, n = 1.*G.number_of_edges(), 1.*G.number_of_nodes()
    degrees = G.degree()   
    nodes = G.nodes()
    
    H = G.to_directed()    
    for u,v in H.edges():
        H[u][v]['cap'] = 1.
    for v in nodes:
        H.add_edge('s', v, cap = m)
        H.add_edge( v, 't', cap = m + 2.*g - 1.*degrees[v])
    eps = 1e-8
    
    while ub - lb > eps:
        cut_value, (S, T) = nx.minimum_cut(H, 's', 't', capacity='cap')
        if len(S) - 1 > 0:
            lb = g
        else:
            ub = g
        g = (lb + ub)/2.
        for v in nodes:            
            H[v]['t']['cap'] = m + 2.*g - 1.*degrees[v]            
    return lb
    
if __name__ == "__main__":
    
    G = nx.erdos_renyi_graph(100, 0.7)
    print 1.*G.number_of_edges()/G.number_of_nodes()
    densest_density = find_densest(G)
    D = get_solution(G, densest_density)
    print 1.*D.number_of_edges()/D.number_of_nodes()
    
      