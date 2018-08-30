from networkx import *
from networkx.generators.random_graphs import *
import numpy as np
from datetime import datetime, timedelta
import os
import random
from copy import deepcopy
    
def plant_uniform(interval, C):

    edges = C.edges()
    random.shuffle(edges)
    TS = []
    t = interval[0]
    for i in xrange(interval[1] + 1 - interval[0]):
        if i == len(edges):
            return TS
        TS.append((t, edges[i]))
        t += 1   
    return TS
    
def plantBackground_uniform(span, G):
    edges = G.edges()
    
    TS = []
    for i in xrange(span):
        if i%len(edges) == 0:
            random.shuffle(edges)
        TS.append((i, edges[i%len(edges)]))    
    return TS
    
        
def generateRG(n, edges_p, nodes):
    G = fast_gnp_random_graph(n, edges_p, seed=None, directed=False)
    gnodes = G.nodes()
    
    relabel = {gnodes[i]: nodes[i] for i in xrange(n)}    
    G = relabel_nodes(G, relabel) 
    return G
    
def generate(generator_pars):
    
    k = generator_pars['k']
    B = generator_pars['B']
    noise = generator_pars['noise']
    innerdegree = generator_pars['innerdegree']
    nodesInCom = generator_pars['nodesInCom']
    backgoundN = generator_pars['backgoundN']
    wholeSpan = generator_pars['wholeSpan']
    
    com_number = k
    
    nCom = nodesInCom
    desired_avgCom = innerdegree

    n = backgoundN
    desired_avg_degree = noise   
        
    start = 0
    span = wholeSpan
    end = start + span

    TS = []
    edges_p = float(desired_avg_degree)/(n-1)
    G = generateRG(n, edges_p, range(n))
    TS = plantBackground_uniform(span, G)

    innerNoise = []
    allocated_interval_length = wholeSpan/k
    generated_C = []
    for i in xrange(com_number):
        interval = (allocated_interval_length*(i), allocated_interval_length*(i)+B-1)
        edges_pCom = float(desired_avgCom)/(nCom - 1)
        C = generateRG(nCom, edges_pCom, range(i*nCom,(i+1)*nCom))
        
        generated_C.append((copy.deepcopy(C), interval))
        TS += plant_uniform(interval, C)  
        t = 2.0*C.number_of_edges()/C.number_of_nodes()
        innerNoise.append(t)
            
    TS.sort()
    TS_out = {}
    for i in TS:
        if i[0] not in TS_out:
            TS_out[i[0]] = []
        TS_out[i[0]].append(i[1])
        
    backNoise = 2.0*G.number_of_edges()/G.number_of_nodes()
    innerNoiseout = np.mean(innerNoise)
    
    return TS_out, backNoise, innerNoiseout, generated_C
