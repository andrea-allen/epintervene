import unittest
import numpy as np
from epintervene.simobjects import simulation


class TestSimulation(unittest.TestCase):
    def test_sim_creates_patientzero(self):
        A = np.random.random_integers(0, 1, (10, 10))
        A = (A + A.T)/2
        sim = simulation.Simulation(A)
        sim.initialize_patient_zero()
        self.assertEqual(len(sim.current_infected), 1)

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
        for edge in sim.current_IS_edges:
            self.assertEqual(edge.event_rate, 0)

        sim.add_infection_event_rates(good_Beta)
        for edge in sim.current_IS_edges:
            self.assertEqual(edge.event_rate, 0.2)
