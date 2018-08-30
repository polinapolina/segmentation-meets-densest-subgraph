import numpy as np
import math
import copy
import networkx as nx
import pickle
import charikar

class IncDensest:
    def __init__(self):
        self.G = nx.Graph()
        self.best_density = 0.5
        
    def __add_new_node(self, n1):
        self.G.add_node(n1)
        return 
     
    def add_edge(self, n1, n2):
        if not self.G.has_node(n1):
            self.__add_new_node(n1)
        if not self.G.has_node(n2):
            self.__add_new_node(n2)           
        if not self.G.has_edge(n1, n2):
            self.G.add_edge(n1,n2)
            self.best_density = charikar.charikar_densest(self.G)
        return self.best_density
 