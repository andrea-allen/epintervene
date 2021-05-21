import unittest
from epintervene.simobjects import network
from epintervene.simobjects import nodestate

#TODO: needs to be edited

class TestNetwork(unittest.TestCase):
    def test_node(self):
        first_node = network.Node(label=1, generation=0, state=nodestate.NodeState.SUSCEPTIBLE, event_rate=.1)
        second_node = network.Node(label=2, generation=1, state=nodestate.NodeState.INFECTED, event_rate=.2)
        self.assertFalse(first_node.equals(second_node))

        first_node.infect()
        self.assertEqual(first_node._state, nodestate.NodeState.INFECTED)
        first_node.recover()
        self.assertEqual(first_node._state, nodestate.NodeState.RECOVERED)

        third_node_same_label = network.Node(label=1, generation=0, state=nodestate.NodeState.SUSCEPTIBLE, event_rate=0.1)
        self.assertTrue(third_node_same_label.equals(first_node))

    def test_edge(self):
        first_node = network.Node(label=1, generation=0, state=nodestate.NodeState.SUSCEPTIBLE, event_rate=.1)
        second_node = network.Node(label=2, generation=1, state=nodestate.NodeState.INFECTED, event_rate=.2)

        edge = network.Edge(first_node, second_node, event_rate=0.5)
        should_be_same_edge = network.Edge(first_node, second_node, event_rate=0.4)

        self.assertTrue(edge.equals(should_be_same_edge))

    def test_edge_infection(self):
        first_node = network.Node(label=1, generation=0, state=nodestate.NodeState.INFECTED, event_rate=.1)
        second_node = network.Node(label=2, generation=-1, state=nodestate.NodeState.SUSCEPTIBLE, event_rate=.2)

        self.assertEqual(first_node.get_state(), nodestate.NodeState.INFECTED)
        self.assertEqual(first_node.get_generation(), 0)
        self.assertEqual(second_node.get_state(), nodestate.NodeState.SUSCEPTIBLE)
        self.assertEqual(second_node.get_generation(), -1)

        edge = network.Edge(first_node, second_node, event_rate=0.5)

        edge.infect()

        self.assertEqual(second_node.get_state(), nodestate.NodeState.INFECTED)
        self.assertEqual(second_node.get_generation(), 1)


if __name__ == '__main__':
    unittest.main()
