import random
import numpy as np
from epintervene.simobjects import network
from epintervene.simobjects import eventtype


class Event:
    def __init__(self, event_name, event_rate):
        # create classes of events such as IS events (infected to susceptible event) and give it a rate parameter
        self.event_rate = event_rate
        self.event_name = event_name


class EventList:
    def __init__(self, event_type, event_list):
        self.event_list = event_list
        self.event_type = event_type

    def add_to_event_list(self, event):
        self.event_list.append(event)

    def remove_from_event_list(self, event):
        self.event_list.remove(event)

    def get_type(self):
        return self.event_type.name


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
        self.Beta = None
        self.Gamma = None
        self.current_IS_edges = EventList(eventtype.EventType.INFECTEDSUSCEPTIBLE, [])
        self.current_infected = EventList(eventtype.EventType.RECOVER, [])
        self.recovered_nodes = []
        self.has_been_infected_labels = []
        self.highest_gen = 0
        self.gen_collection = {}  # I think this should be kept because it's what sets apart this sim framework as unique
        self.active_nodes = []  # add documentation
        self.total_num_timesteps = 0
        self.time_series = [0]
        self.real_time_srs_infc = []
        self.real_time_srs_rec = []
        self.generational_emergence = {0: 0}

    def initialize_patient_zero(self):
        N = len(self.A[0])
        p_zero_idx = random.randint(0, N - 1)
        if self.Gamma is None:
            raise AttributeError(self.Gamma,
                                 'Please provide a recovery vector Gamma to the simulation via the add_recover_event_rates method')
        if self.Beta is None:
            raise AttributeError(self.Beta,
                                 'Please provide an infection rate matrix Beta to the simulation via the add_infection_event_rates method')
        patient_zero = network.Node(p_zero_idx, 0, 1, recover_rate=self.Gamma[p_zero_idx])
        self.active_nodes.append(patient_zero)
        self.gen_collection[0] = [p_zero_idx]
        self.current_infected.add_to_event_list(patient_zero)
        self.has_been_infected_labels.append(p_zero_idx)
        for j in range(0, len(self.A[0])):
            if self.A[p_zero_idx, j] == 1:
                neighbor = network.Node(j, -1, 0, recover_rate=self.Gamma[j])
                self.active_nodes.append(neighbor)
                edge_ij = network.Edge(patient_zero, neighbor,
                                       infect_rate=self.Beta[patient_zero.label, neighbor.label])
                self.current_IS_edges.add_to_event_list(edge_ij)

    def add_infection_event_rates(self, Beta):
        if len(Beta) != self.N:
            raise ValueError(Beta, 'Beta matrix must be N by N to match network size')
        else:
            self.Beta = Beta  # Document: Beta must be an N x N matrix with entries for the probability of infection between each pair of nodes
            for edge in self.current_IS_edges.event_list:
                edge.set_event_rate(self.Beta[edge.left_node.label, edge.right_node.label])

    def add_recover_event_rates(self, Gamma):
        if len(Gamma) != self.N:
            raise ValueError(Gamma, 'Gamma vector must be length N to match network size')
        else:
            self.Gamma = Gamma  # Document: Beta must be an N x N matrix with entries for the probability of infection between each pair of nodes
            for node in self.current_infected.event_list:
                node.set_event_rate(self.Gamma[node.label])

    def run_sim(self):
        self.initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.current_IS_edges.event_list) == 0:
                break

    def single_step(self, visualize=False):
        if visualize:
            self.visualize_network()
        event_catolog = [self.current_IS_edges, self.current_infected]
        tau = draw_tau(event_catolog)

        self.time_series.append(self.time_series[-1] + tau)
        self.real_time_srs_infc.append(len(self.current_infected.event_list))
        self.real_time_srs_rec.append(len(self.recovered_nodes))

        event_class = draw_event_class(event_catolog)  # This is returning a whole list of events
        if event_class is not None:
            next_event = draw_event(
                event_class)  # Have a type that is an EventList class that contains Type, and the list of events
            if event_class.event_type == eventtype.EventType.INFECTEDSUSCEPTIBLE:  # Todo have events know how to do their own acrobatics
                infection_event = next_event
                infection_event.infect()
                self.current_IS_edges.remove_from_event_list(infection_event)
                self.current_infected.add_to_event_list(infection_event.right_node)

                try:
                    self.gen_collection[infection_event.right_node.generation].append(
                        infection_event.right_node.label)  # maybe move toward storing actual node objects? but also this could get huge. Could also append a vector that tracks what real time the node became infected too
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self.gen_collection[infection_event.right_node.generation] = [infection_event.right_node.label]
                    self.highest_gen += 1
                    self.generational_emergence[self.highest_gen] = self.current_sim_time
                self.has_been_infected_labels.append(infection_event.right_node.label)
                self.update_IS_edges()
                self.add_IS_edges(infection_event.right_node)
            if event_class.event_type == eventtype.EventType.RECOVER:
                recovery_event = next_event
                self.current_infected.remove_from_event_list(recovery_event)
                recovery_event.recover()
                self.update_IS_edges()
                self.recovered_nodes.append(recovery_event)
        self.current_sim_time += tau

    def update_IS_edges(self):
        updated_V_IS = []
        for edge in self.current_IS_edges.event_list:
            if (edge.left_node.state == 1) and (edge.right_node.state == 0):
                updated_V_IS.append(edge)
        self.current_IS_edges.event_list = updated_V_IS

    def add_IS_edges(self, infected_node):
        for j in range(0, len(self.A[infected_node.label])):
            if self.A[infected_node.label][j] == 1:
                candidate_node = network.Node(j, -1, 0, self.Gamma[j])
                neighbor_node = self.existing_node(candidate_node)
                if neighbor_node.state == 0:
                    edge_ij = network.Edge(infected_node, neighbor_node, self.Beta[infected_node.label][j])
                    if not self.edge_list_contains(edge_ij):
                        self.current_IS_edges.add_to_event_list(edge_ij)

    def existing_node(self, candidate_node):
        for node in self.active_nodes:
            if candidate_node.equals(node):
                return node
        self.active_nodes.append(candidate_node)
        return candidate_node

    def edge_list_contains(self, edge):
        for e in self.current_IS_edges.event_list:
            if e.equals(edge):
                return True
        return False

    def prune_IS_edges(self):
        for edge in self.current_IS_edges.event_list:
            try:
                edge_exists_in_network = (self.A[edge.left_node.label][edge.right_node.label] == 1)
                if not edge_exists_in_network:
                    self.current_IS_edges.remove_from_event_list(edge)
            except IndexError:
                self.current_IS_edges.remove_from_event_list(
                    edge)  # This happens if the new network no longer contains that node, can remove them

    def visualize_network(self):
        print('In progress')

    def tabulate_generation_results(self, max_gens):
        # todo would be good to also have a "real time" vector that gives the option to plot that version of a timeseries
        # do this next in order to compare time intervention vs not, time steps and real time? just real time probably
        gens = max_gens
        total_infected_time_srs = np.zeros(gens)
        s = 1
        m = 1
        s_max = 1
        for gen in range(max_gens):  # {0: 1, 1: 12, 14, 2: 16, 42, ....
            total_infected_time_srs[gen] = s
            try:
                m = len(self.gen_collection[gen + 1])  # num infected in gen g
                s += m
                s_max = s
            except KeyError:
                # make m=0 and s=the last s for the rest of the "time series"
                s = s_max
                m = 0
        return total_infected_time_srs

    def tabulate_continuous_time(self, time_buckets=100):
        max_time = max(self.time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        infection_time_series = np.zeros(len(time_partition))
        recover_time_series = np.zeros(len(time_partition))
        ts_length = len(self.time_series)

        try:
            for t in range(ts_length - 1):
                real_time_val = self.time_series[t]
                real_time_infected = self.real_time_srs_infc[t]
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        infection_time_series[idx] = real_time_infected
        except IndexError:
            print('No values for infection time series')

        try:
            for t in range(ts_length - 1):
                real_time_val = self.time_series[t]
                real_time_recovered = self.real_time_srs_rec[t]
                # determine what index it goes in
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        recover_time_series[idx] = real_time_recovered
        except IndexError:
            print('No values for recovery time series')

        return time_partition, infection_time_series, recover_time_series


def draw_tau(event_catalog):
    list_of_events = []
    for event_group in event_catalog:
        list_of_events.append(event_group.event_list)
    list_of_events = [item for sublist in list_of_events for item in sublist]
    sum_of_rates = np.sum(event.event_rate for event in list_of_events)
    if sum_of_rates == 0:
        sum_of_rates += .0001
    tau = np.random.exponential(1 / sum_of_rates)
    return tau


def draw_event_class(event_catalog):
    num_event_classes = len(event_catalog)
    partition_end_markers = {}
    total_combined_rate = 0
    for event_class in range(num_event_classes):
        events = event_catalog[event_class].event_list
        total_combined_rate += np.sum(event.event_rate for event in events)
        partition_end_markers[event_class] = total_combined_rate
    random_draw = random.uniform(0, total_combined_rate)
    for i in partition_end_markers.keys():
        if random_draw < partition_end_markers[i]:
            return event_catalog[i]  # Document: return the whole list of events in that class


def draw_event(event_list):
    possible_events = event_list.event_list
    max_rate = max(event.event_rate for event in possible_events)
    accepted = False
    random_event = None
    L = len(event_list.event_list)  # Document: drawing weighted with re-sampling (with replacement)
    while not accepted:
        random_idx = np.random.randint(0, L)
        random_event = event_list.event_list[random_idx]
        accept_rate = random_event.event_rate / max_rate
        random_draw = random.uniform(0, 1)
        if random_draw < accept_rate:
            accepted = True
    return random_event
