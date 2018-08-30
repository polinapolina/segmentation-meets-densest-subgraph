import numpy as np
import networkx as nx
from networkx.utils import powerlaw_sequence
import time
import sys
import pickle
import os.path

reload(sys)  
sys.setdefaultencoding('utf8')
from datetime import datetime, timedelta

    
def generateGraph(n = 100, seed = 1.0):
    generated = False
    while not generated:
        try:
            z = nx.utils.create_degree_sequence(n, powerlaw_sequence)
            G = nx.configuration_model(z, seed = seed)
            generated = True
        except:
            pass
    G = nx.Graph(G)
    G.remove_edges_from(G.selfloop_edges())
    return G
    
def readdata(filename, unix = True):
    timestamps = []
    with open(filename, 'r') as f:    
        for line in f:                
            line = line.strip().split(' ')
            if not unix:
                tstr =  line[0][1:] + ' ' + line[1][0:-1]
                t = datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S')
                timestamp = time.mktime(datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S').timetuple())
                
                tst, n1, n2 = int(timestamp), line[2], line[3]
            else:
                timestamp = int(line[0])
                tst, n1, n2 = int(line[0]), int(line[1]), int(line[2])
            
            if n1 == n2:
                continue
            
            if n2 < n1:
                n1, n2 = n2, n1            
            timestamps.append((tst, n1, n2))            
        timestamps.sort()
    return timestamps
    
def readdata_dict(filename, unix = True):
    timestamps = {}
    with open(filename, 'r') as f:    
        for line in f:                
            line = line.strip().split(' ')
            if not unix:
                tstr =  line[0][1:] + ' ' +line[1][0:-1]
                t = datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S')
                timestamp = time.mktime(datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S').timetuple())
                
                tst, n1, n2 = int(timestamp), line[2], line[3]
            else:
                timestamp = int(line[0])
                tst, n1, n2 = int(line[0]), int(line[1]), int(line[2])
            
            if n1 == n2:
                continue
            
            if n2 < n1:
                n1, n2 = n2, n1
            if tst not in timestamps:
                timestamps[tst] = []
            timestamps[tst].append((n1, n2))
    return timestamps
    
def readdata_dict_limit(filename, limit = 1000, unix = True):
    timestamps = {}
    with open(filename, 'r') as f:
        c = 0
        for line in f:                
            line = line.strip().split(' ')
            if not unix:
                tstr =  line[0][1:] + ' ' +line[1][0:-1]
                t = datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S')
                timestamp = time.mktime(datetime.strptime(tstr, '%Y-%m-%d %H:%M:%S').timetuple())
                
                tst, n1, n2 = int(timestamp), line[2], line[3]
            else:
                timestamp = int(line[0])
                tst, n1, n2 = int(line[0]), int(line[1]), int(line[2])
            
            if n1 == n2:
                continue
            
            if n2 < n1:
                n1, n2 = n2, n1
            if tst not in timestamps:
                timestamps[tst] = []
            timestamps[tst].append((n1, n2))
            c += 1
            if c >= limit:
                return timestamps
    return timestamps

    
def readgraph(filename, unix = True):
    G = nx.Graph()
    with open(filename, 'r') as f:    
        for line in f:                
            line = line.strip().split(' ')
            n1, n2 = int(line[2]), int(line[3])            
            if n1 != n2:
                if n2 < n1:
                    n1, n2 = n2, n1  
                G.add_edge(n1, n2)
    return G        

def generateTS(G, len = 10000):
    timestamps = []
    edges = G.edges()
    idx = np.random.choice(G.number_of_edges(), len)
    for tst in xrange(len):        
        n1, n2 = edges[idx[tst]]
        if n2 < n1:
            n1, n2 = n2, n1            
        timestamps.append((tst, n1, n2))
    return timestamps
    
def generateTS_dict(G, len = 10000):
    timestamps = {}
    edges = G.edges()
    idx = np.random.choice(G.number_of_edges(), len)
    for tst in xrange(len):
        n1, n2 = edges[idx[tst]]
        if n2 < n1:
            n1, n2 = n2, n1
        if tst not in timestamps:
            timestamps[tst] = []
        timestamps[tst].append((n1, n2))
    return timestamps
        
    
def generateTS_planted(G, len = 100, k = 10, clique_size = 5, interval = 10):
    timestamps = []
    edges = G.edges()
    idx = np.random.choice(G.number_of_nodes(), len)
    for tst in xrange(len):        
        n1, n2 = edges[idx[tst]]
        if n2 < n1:
            n1, n2 = n2, n1            
        timestamps.append((tst, n1, n2))
    
    cl = nx.complete_graph(clique_size)    
    start = interval
    max_ID = max(G.nodes())  
    tst = start
    for i in xrange(k):
        cl_tmp = nx.relabel_nodes(cl, {n: i+max_ID+1 for (i, n) in enumerate(cl.nodes())})
        for (n1, n2) in cl_tmp.edges():
            if n2 < n1:
                n1, n2 = n2, n1
            timestamps.append((tst, n1, n2))
        tst += interval
        max_ID = max(cl_tmp.nodes())
        
    timestamps.sort()
    return timestamps
