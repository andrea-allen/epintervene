import unittest
from epintervene.simobjects import network


class TestNetwork(unittest.TestCase):
    def test_node(self):
        first_node = network.Node(label=1, generation=0, state=0, recover_rate=.1)
        second_node = network.Node(label=2, generation=1, state=1, recover_rate=.2)
        self.assertFalse(first_node.equals(second_node))

        first_node.infect()
        self.assertEqual(first_node.state, 1)
        first_node.recover()
        self.assertEqual(first_node.state, 2)

        third_node_same_label = network.Node(label=1, generation=0, state=0, recover_rate=0.1)
        self.assertTrue(third_node_same_label.equals(first_node))

    def test_edge(self):
        first_node = network.Node(label=1, generation=0, state=0, recover_rate=.1)
        second_node = network.Node(label=2, generation=1, state=1, recover_rate=.2)

        edge = network.Edge(first_node, second_node, infect_rate=0.5)
        should_be_same_edge = network.Edge(first_node, second_node, infect_rate=0.4)

        self.assertTrue(edge.equals(should_be_same_edge))

    def test_edge_infection(self):
        first_node = network.Node(label=1, generation=0, state=1, recover_rate=.1)
        second_node = network.Node(label=2, generation=-1, state=0, recover_rate=.2)

        self.assertEqual(first_node.state, 1)
        self.assertEqual(first_node.generation, 0)
        self.assertEqual(second_node.state, 0)
        self.assertEqual(second_node.generation, -1)

        edge = network.Edge(first_node, second_node, infect_rate=0.5)

        edge.infect()

        self.assertEqual(second_node.state, 1)
        self.assertEqual(second_node.generation, 1)


if __name__ == '__main__':
    unittest.main()
