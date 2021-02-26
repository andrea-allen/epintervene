import time
import random
import numpy as np
from epintervene.simobjects import network


class Event:
    def __init__(self, event_name, event_rate):
        # create classes of events such as IS events (infected to susceptible event) and give it a rate parameter
        self.event_rate = event_rate
        self.event_name = event_name


# todo planning: for each type of event, you're gonna want a matrix or list of event rates
# be able to add different event classes and then build off of those, when you add an event class it also adds
# a time series tracker
# time series should increment at pre-determined steps or benchmarks

# Simulation base class either has nothing, or is SI
# options to choose from SEIR or whatever else or add classes as I see fit
class Simulation:
    def __init__(self, adj_matrix, max_unitless_sim_time=1000000):
        self.total_sim_time = max_unitless_sim_time
        self.current_sim_time = 0
        self.A = adj_matrix
        self.N = len(self.A[0])
        self.current_IS_edges = []
        self.current_infected = []
        self.has_been_infected_labels = []
        self.gen_collection = {} #I think this should be kept because it's what sets apart this sim framework as unique
        self.active_nodes = [] #add documentation
        self.total_num_timesteps = 0
        self.time_series = [0]
        self.real_time_srs_infc = []
        self.real_time_srs_rec = []
        self.generational_emergence = {0: 0}

    def initialize_patient_zero(self):
        N = len(self.A[0])
        p_zero_idx = random.randint(0, N - 1)
        patient_zero = network.Node(p_zero_idx, 0, 1, recover_rate=0)
        self.active_nodes.append(patient_zero)
        self.gen_collection[0] = [p_zero_idx]
        self.current_infected.append(patient_zero)
        self.has_been_infected_labels.append(p_zero_idx)
        for j in range(0, len(self.A[0])):
            if self.A[p_zero_idx, j] == 1:
                neighbor = network.Node(j, -1, 0, recover_rate=0)
                self.active_nodes.append(neighbor)
                edge_ij = network.Edge(patient_zero, neighbor, infect_rate=0)
                self.current_IS_edges.append(edge_ij)

    def add_infection_event_rates(self, Beta):
        if len(Beta) != self.N:
            raise ValueError(Beta, 'Beta matrix must be N by N to match network size')
        else:
            self.Beta = Beta #Document: Beta must be an N x N matrix with entries for the probability of infection between each pair of nodes
            for edge in self.current_IS_edges:
                edge.set_event_rate(self.A[edge.left_node.label, edge.right_node.label])