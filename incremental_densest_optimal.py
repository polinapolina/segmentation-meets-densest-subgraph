import numpy as np
import math
import fibonacci_heap_mod as fb
import copy
import networkx as nx
import pickle
import optimal_max_density

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
            self.best_density = optimal_max_density.find_densest(self.G, self.best_density)
        return self.best_density
 