import unittest
import numpy as np
from epintervene import simulator


class TestSimulator(unittest.TestCase):
    def test_simulator(self):
        simulator.run()
        output_data = np.loadtxt('../data/sample_output.txt', delimiter='\n')
        self.assertListEqual(list(output_data), [1, 2, 3])


if __name__ == '__main__':
    unittest.main()
