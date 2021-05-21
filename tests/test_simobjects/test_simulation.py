import unittest
import numpy as np
from epintervene.simobjects import simulation
import networkx as nx
import matplotlib.pyplot as plt

#TODO needs to be edited

class TestSimulation(unittest.TestCase):
    def test_sim_creates_patientzero(self):
        A = np.random.random_integers(0, 1, (10, 10))
        A = (A + A.T)/2
        sim = simulation.Simulation(A)
        self.assertRaises(AttributeError, sim._initialize_patient_zero)

        good_Beta = np.full((10,10), 0.2)
        good_Gamma = np.full(10, 0.5)
        sim.add_recover_event_rates(good_Gamma)
        sim.add_infection_event_rates(good_Beta)

        sim._initialize_patient_zero()
        self.assertEqual(len(sim._potential_recovery_events._event_list), 1)

    def test_beta_rate_matrix_throws_exception(self):
        A = np.random.random_integers(0, 1, (10, 10))
        A = (A + A.T)/2
        sim = simulation.Simulation(A)
        bad_Beta = np.zeros((2,2))
        self.assertRaises(ValueError, sim.add_infection_event_rates, bad_Beta)

    def test_beta_rate_matrix_alters_event_rates(self):
        A = np.random.random_integers(0, 1, (10, 10))
        A = (A + A.T)/2
        sim = simulation.Simulation(A)
        good_Beta = np.full((10,10), 0.2)
        good_Gamma = np.full(10, 0.5)
        for edge in sim._potential_IS_events._event_list:
            self.assertEqual(edge.event_rate, 0)

        sim.add_infection_event_rates(good_Beta)
        sim.add_recover_event_rates(good_Gamma)
        for edge in sim._potential_IS_events._event_list:
            self.assertEqual(edge.event_rate, 0.2)

    def test_IS_edges_are_updated_after_single_step(self):
        N = 10
        self.beta = 0.99
        self.gamma = 0.0001
        graph = nx.generators.complete_graph(N)
        N = len(graph.nodes())
        self.Beta = np.full((N, N), self.beta)
        self.Gamma = np.full(N, self.gamma)
        adjacency_matrix = np.array(nx.adjacency_matrix(graph).todense())
        sim = simulation.Simulation(adjacency_matrix)
        sim.add_infection_event_rates(self.Beta)
        sim.add_recover_event_rates(self.Gamma)
        sim._initialize_patient_zero()
        print('before single step')

        self.assertGreaterEqual(len(sim._potential_IS_events._event_list), 1)
        self.assertEqual(len(sim._potential_recovery_events._event_list), 1)
        sim._single_step()

        print('after single step')
        self.assertEqual(len(sim._potential_recovery_events._event_list), 2)
        infected_nodes = list(map(lambda node: node.get_label(), sim._potential_recovery_events._event_list))
        for edge in sim._potential_IS_events._event_list:
            self.assertIn(edge.get_left_node().get_label(), infected_nodes)
            self.assertNotIn(edge.get_right_node().get_label(), infected_nodes)
        if len(sim._potential_IS_events._event_list) > 0:
            sim._single_step()
            print('after second single step')
            self.assertEqual(len(sim._potential_recovery_events._event_list), 3)
            infected_nodes = list(map(lambda node: node.get_label(), sim._potential_recovery_events._event_list))
            for edge in sim._potential_IS_events._event_list:
                self.assertIn(edge.get_left_node().get_label(), infected_nodes)
                self.assertNotIn(edge.get_right_node().get_label(), infected_nodes)

    def test_custom_time_series_results(self):
        N = 100
        self.beta = 0.99
        self.gamma = 0.0001
        graph = nx.generators.erdos_renyi_graph(N, 0.02)
        N = len(graph.nodes())
        self.Beta = np.full((N, N), self.beta)
        self.Gamma = np.full(N, self.gamma)
        adjacency_matrix = np.array(nx.adjacency_matrix(graph).todense())
        sim = simulation.Simulation(adjacency_matrix)
        sim.add_infection_event_rates(self.Beta)
        sim.add_recover_event_rates(self.Gamma)
        sim._initialize_patient_zero()

        sim.run_sim()

        custom_time_results = sim.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=15)
        self.assertEqual(len(custom_time_results[0]), 1000)
        self.assertEqual(max(custom_time_results[0]), 15-0.015)

if __name__ == '__main__':
    unittest.main()
