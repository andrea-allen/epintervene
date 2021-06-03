import numpy as np
import networkx as nx
import math
from epintervene.simobjects import nodestate
import matplotlib.pyplot as plt


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

    def vaccinate(self):
        self._state = nodestate.NodeState.VACCINATED

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

    @staticmethod
    def erdos_renyi(N, p):
        return nx.erdos_renyi_graph(N, p)

    @staticmethod
    def create_adjacency_list(G):
        adjlist = list(list(map(int, minilist)) for minilist in list(map(str.split, nx.generate_adjlist(G))))

        corrected_adjlist = list(list() for i in range(len(adjlist)))
        for entry in adjlist:
            try:
                corrected_adjlist[entry[0]] = entry
            except IndexError:
                for i in range(entry[0]-len(adjlist) + 1):
                    corrected_adjlist.append(list())
                corrected_adjlist[entry[0]] = entry
                print(entry)

        adjlist = corrected_adjlist
        # make adjlist symmetric:
        len_adj = len(adjlist)
        current_source = 0
        for row in adjlist:
            try:
                source = row[0]
                current_source += 1
            except IndexError:
                row.append(current_source)
                source = row[0]
                current_source += 1
                print(row)
            if len(row) > 1:
                for target in row[1:]:
                    if sum([1 for entry in row if entry==target]) > 1:
                        row.remove(target)
                    try:
                        if source not in adjlist[target]:
                            adjlist[target].append(source)
                    except IndexError:
                        print(source, target)
                        adjlist.append([source])
                # unique_neighbors = list(np.unique(row[1:]))
                # adjlist[-1][1:] = unique_neighbors
        total_entry_count = 0
        for i in range(len_adj):
            if len(adjlist[i]) > 1:
                total_entry_count += len(adjlist[i][1:])
        return adjlist

def visualize(N, graph, pos, gen_collection):
    G = graph
    val_map = {}
    gen_labels_map = {}
    vmax=10
    max_gen = max(gen_collection.keys()) + 1
    for gen in gen_collection.keys():
        nodes = gen_collection[gen]
        for node in nodes:
            val_map[node] = vmax - (gen)
            gen_labels_map.update({node: gen})

    values = [val_map.get(node, 0) for node in G.nodes()]

    nx.draw_networkx_nodes(G, pos=pos, cmap=plt.get_cmap('YlGnBu'), node_color=values, vmin=0, vmax=vmax)
    # nx.draw_networkx_nodes(G, pos=pos, node_color=)
    # gen_labels_map[1]=''
    # gen_labels_map[10]=''
    # gen_labels_map[20]=''
    # gen_labels_map[30]=''
    # nx.draw_networkx_labels(G, pos=pos, with_labels=True, labels=gen_labels_map)
    nx.draw_networkx_edges(G, pos=pos, edge_color='grey', lw=2)

    # V = nx.Graph()
    # V.add_nodes_from([1, 10])
    # V2 = nx.Graph()
    # V2.add_nodes_from([20,30])
    # nx.draw_networkx_nodes(V, pos=pos, node_color='red')
    # nx.draw_networkx_nodes(V2, pos=pos, node_color='orange')
    # nx.draw_networkx_labels(V, pos=pos, with_labels=False, labels={1:'V', 10:'V', 20:'V', 30:'V'})
    plt.axis(False)
    plt.tight_layout()
    plt.show()

