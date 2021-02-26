import numpy as np
import networkx as nx
import math


class Node:
    def __init__(self, label, generation, state, recover_rate):
        self.generation = generation
        self.label = label
        self.state = state
        self.event_rate = recover_rate

    def infect(self):
        self.state = 1

    def recover(self):
        self.state = 2

    def set_generation(self, g):
        self.generation = g

    def set_recover_rate(self, recover_rate):
        self.event_rate = event_rate

    def display_info(self):
        print('Node index: ', self.label, ' state: ', self.state, ' event_rate: ', self.event_rate, ' gen: ',
              self.generation)

    def equals(self, node):
        if self.label == node.label:
            return True
        else:
            return False


class Edge:
    def __init__(self, left_node, right_node, infect_rate):
        self.left_node = left_node  # not just an index, this is a whole Node object
        self.right_node = right_node
        self.event_rate = infect_rate

    def infect(self):
        self.right_node.infect()
        self.right_node.set_generation(self.left_node.generation + 1)

    def set_event_rate(self, event_rate):
        self.event_rate = event_rate

    def display_info(self):
        print('Edge with event rate: ', self.event_rate, ' nodes:')
        self.left_node.display_info()
        self.right_node.display_info()

    def equals(self, other_edge):
        # Imperative to use Node class equality here
        if self.left_node.equals(other_edge.left_node) and self.right_node.equals(other_edge.right_node):
            return True
        else:
            return False


class NetworkBuilder:
    def __init__(self, N):
        self.N = N

    @staticmethod
    def from_degree_distribution(N, degree_dist, return_pos=False):
        number_of_nodes = N * np.array(degree_dist)
        degree_sequence = []
        for i in range(int(math.floor(len(number_of_nodes)))):
            number_with_that_degree = number_of_nodes[i]
            for k in range(int(math.floor(number_with_that_degree))):
                degree_sequence.append(i)
        graphical = nx.is_graphical(degree_sequence)
        if not graphical:
            degree_sequence.append(1)
        G = nx.configuration_model(degree_sequence)
        try:
            G.remove_edges_from(nx.selfloop_edges(G))
        except RuntimeError:
            print('No self loops to remove')
        if return_pos:
            pos = nx.spring_layout(G)
            return G, pos
        return G, None

    @staticmethod
    def from_adjacency_matrix(A, return_pos=False):
        G = nx.from_numpy_matrix(A)
        if return_pos:
            pos = nx.spring_layout(G)
            return G, pos
        return G, None

