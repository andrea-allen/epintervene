import unittest
import numpy as np
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
import networkx as nx
import matplotlib.pyplot as plt


class TestRandomRolloutSimulation(unittest.TestCase):
    def test_vaccinated_nodes_are_not_done_twice(self):
        N = 10
        self.beta = 0.99
        self.gamma = 0.0001
        graph = nx.generators.complete_graph(N)
        N = len(graph.nodes())
        self.Beta = np.full((N, N), self.beta)
        self.Gamma = np.full(N, self.gamma)
        adjacency_matrix = np.array(nx.adjacency_matrix(graph).todense())
        sim = extended_simulation.RandomRolloutSimulation(adjacency_matrix)
        sim.add_infection_event_rates(self.Beta)
        sim.add_recover_event_rates(self.Gamma)
        sim.configure_intervention([1, 2], [0.0, 0.0], [.9, .1])

        first_round_vaccinated_nodes = sim.intervene(0)
        node_labels = [node._label for node in first_round_vaccinated_nodes]
        first_round_unique_nodes = np.unique(node_labels)
        self.assertEqual(len(first_round_unique_nodes), 9)

        second_round_vaccinated_nodes = sim.intervene(1)
        second_round_node_labels = [node._label for node in second_round_vaccinated_nodes]
        second_round_unique_nodes = np.unique(second_round_node_labels)
        self.assertEqual(len(second_round_unique_nodes), 1)
        for nl in second_round_unique_nodes:
            self.assertNotIn(nl, first_round_unique_nodes)



if __name__ == '__main__':
    unittest.main()