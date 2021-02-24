

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

    def display_info(self):
        print('Node index: ', self.label, ' state: ', self.state, ' event_rate: ', self.event_rate, ' gen: ', self.generation)

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