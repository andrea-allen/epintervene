import numpy as np
import networkx as nx
import math
from epintervene.simobjects import nodestate


class Node:
    def __init__(self, label, generation, state, event_rate):
        self._generation = generation
        self._label = label
        self._state = state
        self._event_rate = event_rate
        self._membership = None

    def infect(self):
        self._state = nodestate.NodeState.INFECTED

    def recover(self):
        self._state = nodestate.NodeState.RECOVERED

    def expose(self):
        self._state = nodestate.NodeState.EXPOSED

    def get_label(self):
        return self._label

    def get_generation(self):
        return self._generation

    def get_state(self):
        return self._state

    def get_event_rate(self):
        return self._event_rate

    def get_membership(self):
        return self._membership

    def set_generation(self, g):
        self._generation = g

    def set_event_rate(self, event_rate):
        self._event_rate = event_rate

    def display_info(self):
        print('Node index: ', self._label, ' state: ', self._state, ' event_rate: ', self._event_rate, ' gen: ',
              self._generation, 'membership: ', self._membership)

    # Give Nodes an option for network class membership, for example in multilayer network or SBM
    # Because Python isn't strongly typed, this can be an int, a float, a string or an Enum
    def set_membership(self, membership):
        self._membership = membership

    def equals(self, node):
        if self._label == node.get_label():
            return True
        else:
            return False


class Edge:
    def __init__(self, left_node, right_node, event_rate):
        self._left_node = left_node  # not just an index, this is a whole Node object
        self._right_node = right_node
        self._event_rate = event_rate

    def infect(self):
        self._right_node.infect()
        self._right_node.set_generation(self._left_node.get_generation() + 1)

    def expose(self):
        self._right_node.expose()
        self._right_node.set_generation(self._left_node.get_generation() + 1)

    def set_event_rate(self, event_rate):
        self._event_rate = event_rate

    def get_event_rate(self):
        return self._event_rate

    def get_right_node(self):
        return self._right_node

    def get_left_node(self):
        return self._left_node

    def display_info(self):
        print('Edge with event rate: ', self._event_rate, ' nodes:')
        self._left_node.display_info()
        self._right_node.display_info()

    def equals(self, other_edge):
        # Imperative to use Node class equality here
        if self._left_node.equals(other_edge._left_node) and self._right_node.equals(other_edge._right_node):
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

