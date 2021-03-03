import random
import numpy as np
from epintervene.simobjects import network
from epintervene.simobjects import eventtype
from epintervene.simobjects import nodestate


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


class Simulation:
    def __init__(self, adj_matrix, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        self.total_sim_time = max_unitless_sim_time
        self.current_sim_time = 0
        self.A = adj_matrix
        self.N = len(self.A[0])
        self.Beta = None
        self.Gamma = None
        self.potential_IS_events = EventList(eventtype.EventType.INFECTEDSUSCEPTIBLE, [])
        self.potential_recovery_events = EventList(eventtype.EventType.RECOVER, [])
        self.recovered_nodes = []
        self.current_infected_nodes = []
        self.highest_gen = 0
        self.gen_collection = {}
        self.active_nodes = []
        self.total_num_timesteps = 0
        self.time_series = [0]
        self.real_time_srs_infc = []
        self.real_time_srs_rec = []
        self.generational_emergence = {0: 0}
        self.membership_groups = membership_groups
        self.membership_time_series_infc = {0: []}
        self.membership_time_series_rec = {0: []}
        self.track_memberships = False
        self.node_memberships = node_memberships

    def initialize_patient_zero(self):
        N = len(self.A[0])
        p_zero_idx = random.randint(0, N - 1)
        if self.Gamma is None:
            raise AttributeError(self.Gamma,
                                 'Please provide a recovery vector Gamma to the simulation via the '
                                 'add_recover_event_rates method')
        if self.Beta is None:
            raise AttributeError(self.Beta,
                                 'Please provide an infection rate matrix Beta to the simulation via the '
                                 'add_infection_event_rates method')
        patient_zero = network.Node(p_zero_idx, 0, nodestate.NodeState.INFECTED, event_rate=self.Gamma[p_zero_idx])
        if self.track_memberships:
            patient_zero.set_membership(self.node_memberships[p_zero_idx])
        self.active_nodes.append(patient_zero)
        self.current_infected_nodes.append(patient_zero)
        self.gen_collection[0] = [p_zero_idx]
        self.potential_recovery_events.add_to_event_list(patient_zero)
        for j in range(0, len(self.A[0])):
            if self.A[p_zero_idx, j] == 1:
                neighbor = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, event_rate=self.Gamma[j])
                self.active_nodes.append(neighbor)
                edge_ij = network.Edge(patient_zero, neighbor,
                                       event_rate=self.Beta[patient_zero.get_label(), neighbor.get_label()])
                self.potential_IS_events.add_to_event_list(edge_ij)

    def add_infection_event_rates(self, Beta):
        if len(Beta) != self.N:
            raise ValueError(Beta, 'Beta matrix must be N by N to match network size')
        else:
            self.Beta = Beta  # Document: Beta must be an N x N matrix with entries for the probability of infection
            # between each pair of nodes
            for edge in self.potential_IS_events.event_list:
                edge.set_event_rate(self.Beta[edge.get_left_node().get_label(), edge.get_right_node().get_label()])

    def add_recover_event_rates(self, Gamma):
        if len(Gamma) != self.N:
            raise ValueError(Gamma, 'Gamma vector must be length N to match network size')
        else:
            self.Gamma = Gamma  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self.potential_recovery_events.event_list:
                node.set_event_rate(self.Gamma[node.get_label()])

    def run_sim(self, with_memberships=False):
        # TODO function to clear!
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self.initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if len(self.potential_IS_events.event_list) == 0:
                break

    def single_step(self, visualize=False):
        if visualize:
            self.visualize_network()
        event_catolog = [self.potential_IS_events, self.potential_recovery_events]
        tau = draw_tau(event_catolog)

        self.time_series.append(self.time_series[-1] + tau)
        self.real_time_srs_infc.append(len(self.current_infected_nodes))
        self.real_time_srs_rec.append(len(self.recovered_nodes))
        if self.track_memberships:
            self.record_membership_states()

        event_class = draw_event_class(event_catolog)  # This is returning a whole list of events
        if event_class is not None:
            next_event = draw_event(event_class)
            if event_class.event_type == eventtype.EventType.INFECTEDSUSCEPTIBLE:
                infection_event = next_event
                infection_event.infect()
                self.potential_IS_events.remove_from_event_list(infection_event)
                self.potential_recovery_events.add_to_event_list(infection_event.get_right_node())
                self.current_infected_nodes.append(infection_event.get_right_node())

                try:
                    self.gen_collection[infection_event.get_right_node().get_generation()].append(
                        infection_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self.gen_collection[infection_event.get_right_node().get_generation()] = [
                        infection_event.get_right_node().get_label()]
                    self.highest_gen += 1
                    self.generational_emergence[self.highest_gen] = self.current_sim_time
                self.update_IS_events()
                self.add_IS_events(infection_event.get_right_node())
            if event_class.event_type == eventtype.EventType.RECOVER:
                recovery_event = next_event
                self.potential_recovery_events.remove_from_event_list(recovery_event)
                try:
                    self.current_infected_nodes.remove(recovery_event)
                except:
                    print(recovery_event.get_label())
                recovery_event.recover()
                self.update_IS_events()
                self.recovered_nodes.append(recovery_event)
        self.current_sim_time += tau

    def update_IS_events(self):
        updated_IS_events = []
        for edge in self.potential_IS_events.event_list:
            if (edge.get_left_node().get_state() == nodestate.NodeState.INFECTED) \
                    and (edge.get_right_node().get_state() == nodestate.NodeState.SUSCEPTIBLE):
                updated_IS_events.append(edge)
        self.potential_IS_events.event_list = updated_IS_events

    def add_IS_events(self, infected_node):
        for j in range(0, len(self.A[infected_node.get_label()])):
            if self.A[infected_node.get_label()][j] == 1:
                candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self.Gamma[j])
                if self.track_memberships:
                    candidate_node.set_membership(self.node_memberships[candidate_node.get_label()])
                neighbor_node = self.existing_node(candidate_node)
                if neighbor_node.get_state() == nodestate.NodeState.SUSCEPTIBLE:
                    edge_ij = network.Edge(infected_node, neighbor_node, self.Beta[infected_node.get_label()][j])
                    if not self.edge_list_contains(edge_ij):
                        self.potential_IS_events.add_to_event_list(edge_ij)

    def existing_node(self, candidate_node):
        for node in self.active_nodes:
            if candidate_node.equals(node):
                return node
        self.active_nodes.append(candidate_node)
        return candidate_node

    def edge_list_contains(self, edge):
        for e in self.potential_IS_events.event_list:
            if e.equals(edge):
                return True
        return False

    def prune_IS_edges(self):
        for edge in self.potential_IS_events.event_list:
            try:
                edge_exists_in_network = (
                        self.A[edge.get_left_node().get_label()][edge.get_right_node().get_label()] == 1)
                if not edge_exists_in_network:
                    self.potential_IS_events.remove_from_event_list(edge)
            except IndexError:
                self.potential_IS_events.remove_from_event_list(
                    edge)  # This happens if the new network no longer contains that node, can remove them

    def record_membership_states(self):
        for group in self.membership_groups:
            infected_in_group = list(node for node in self.current_infected_nodes if node.get_membership() == group)
            self.membership_time_series_infc[group].append(len(infected_in_group))

    def init_membership_state_time_series(self):
        self.membership_time_series_infc = {}
        for group in self.membership_groups:
            self.membership_time_series_infc[group] = []
            self.membership_time_series_rec[group] = []

    def visualize_network(self):
        print('In progress')

    def tabulate_generation_results(self, max_gens):
        gens = max_gens
        infection_time_srs_by_gen = np.zeros(gens)
        s = 1
        m = 1
        s_max = 1
        for gen in range(max_gens):  # {0: 1, 1: 12, 14, 2: 16, 42, ....
            infection_time_srs_by_gen[gen] = s
            try:
                m = len(self.gen_collection[gen + 1])  # num infected in gen g
                s += m
                s_max = s
            except KeyError:
                # make m=0 and s=the last s for the rest of the "time series"
                s = s_max
                m = 0
        return infection_time_srs_by_gen

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

    def tabulate_continuous_time_with_groups(self, time_buckets=100):
        max_time = max(self.time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        infection_time_series = {}
        for group in self.membership_groups:
            infection_time_series[group] = np.zeros(len(time_partition))
        # TODO tbd grouped recovery
        # recover_time_series = np.zeros(len(time_partition))
        ts_length = len(self.time_series)

        try:
            for t in range(ts_length - 1):
                for group in self.membership_groups:
                    real_time_val = self.time_series[t]
                    real_time_infected = self.membership_time_series_infc[group][t]
                    for idx in range(time_buckets):
                        upper_val = time_partition[idx]
                        if real_time_val < upper_val:
                            infection_time_series[group][idx] = real_time_infected
        except IndexError:
            print('No values for infection time series')

            # TODO tbd for recovery
            # try:
            #     for t in range(ts_length - 1):
            #         real_time_val = self.time_series[t]
            #         real_time_recovered = self.real_time_srs_rec[t]
            #         # determine what index it goes in
            #         for idx in range(time_buckets):
            #             upper_val = time_partition[idx]
            #             if real_time_val < upper_val:
            #                 recover_time_series[idx] = real_time_recovered
            # except IndexError:
            print('No values for recovery time series')

        return time_partition, infection_time_series


class SimulationSEIR(Simulation):
    def __init__(self, adjmatrix, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(adjmatrix, max_unitless_sim_time, membership_groups=membership_groups, node_memberships=node_memberships)
        self.potential_ES_events = EventList(eventtype.EventType.EXPOSEDSUSCEPTIBLE, [])
        self.potential_EI_events = EventList(eventtype.EventType.EXPOSEDINFECTED, [])
        self.current_exposed = []
        self.real_time_srs_exp = []
        self.Beta_ExposedSusceptible = None
        self.Theta_ExposedInfected = None
        self.cont_time_exposed_dict = {}
        self.membership_time_series_exp = {0: []}
        # TODO deal with generations of exposed emergence as well
        # should probably also determine in the model if there's completely asymptomatic cases as well (E->R events)

    def add_exposed_event_rates(self, Beta_Exp):
        if len(Beta_Exp) != self.N:
            raise ValueError(Beta_Exp, 'Beta_Exp matrix must be size N by N to match network size')
        else:
            self.Beta_ExposedSusceptible = Beta_Exp  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self.potential_ES_events.event_list:
                node.set_event_rate(self.Beta_ExposedSusceptible[node.get_label()])

    def add_exposed_infected_event_rates(self, Theta):
        if len(Theta) != self.N:
            raise ValueError(Theta, 'Theta ExposedInfectious vector must be length N to match network size')
        else:
            self.Theta_ExposedInfected = Theta  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self.potential_ES_events.event_list:
                node.set_event_rate(self.Beta_ExposedSusceptible[node.get_label()])

    def run_sim(self, with_memberships=False):
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self.initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            # Run one step
            self.single_step()

            self.total_num_timesteps += 1
            if (len(self.potential_IS_events.event_list) == 0) \
                    and (len(self.potential_ES_events.event_list) == 0) \
                    and (len(self.potential_EI_events.event_list) == 0):
                break

    def single_step(self, visualize=False):
        if visualize:
            self.visualize_network()
        event_catolog = [self.potential_IS_events, self.potential_recovery_events, self.potential_ES_events,
                         self.potential_EI_events]
        tau = draw_tau(event_catolog)

        self.time_series.append(self.time_series[-1] + tau)
        self.real_time_srs_infc.append(len(self.current_infected_nodes))
        self.real_time_srs_rec.append(len(self.recovered_nodes))
        self.real_time_srs_exp.append(len(self.current_exposed))
        if self.track_memberships:
            self.record_membership_states()

        event_class = draw_event_class(event_catolog)  # This is returning a whole list of events
        if event_class is not None:
            next_event = draw_event(event_class)
            if event_class.event_type == eventtype.EventType.INFECTEDSUSCEPTIBLE:
                infection_event = next_event
                infection_event.expose()
                self.potential_IS_events.remove_from_event_list(infection_event)
                self.current_exposed.append(infection_event.get_right_node())
                infection_event.get_right_node().set_event_rate(
                    self.Theta_ExposedInfected[infection_event.get_right_node().get_label()])
                self.potential_EI_events.add_to_event_list(infection_event.get_right_node())
                try:
                    self.gen_collection[infection_event.get_right_node().get_generation()].append(
                        infection_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self.gen_collection[infection_event.get_right_node().get_generation()] = [
                        infection_event.get_right_node().get_label()]
                    self.highest_gen += 1
                    self.generational_emergence[self.highest_gen] = self.current_sim_time
                self.update_IS_events()
                self.update_ES_events()
                self.add_ES_events(infection_event.get_right_node())  # todo this should be add ES edges actually
            elif event_class.event_type == eventtype.EventType.EXPOSEDSUSCEPTIBLE:
                exposure_event = next_event
                exposure_event.expose()
                self.potential_ES_events.remove_from_event_list(exposure_event)
                self.current_exposed.append(exposure_event.get_right_node())
                exposure_event.get_right_node().set_event_rate(
                    self.Theta_ExposedInfected[exposure_event.get_right_node().get_label()])
                self.potential_EI_events.add_to_event_list(
                    exposure_event.get_right_node())  # TODO modify rate above here
                try:
                    self.gen_collection[exposure_event.get_right_node().get_generation()].append(
                        exposure_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self.gen_collection[exposure_event.get_right_node().get_generation()] = [
                        exposure_event.get_right_node().get_label()]
                    self.highest_gen += 1
                    self.generational_emergence[self.highest_gen] = self.current_sim_time
                self.update_IS_events()
                self.update_ES_events()
                self.add_ES_events(exposure_event.get_right_node())

            elif event_class.event_type == eventtype.EventType.EXPOSEDINFECTED:
                exposed_infected_event = next_event
                self.potential_EI_events.remove_from_event_list(exposed_infected_event)
                exposed_infected_event.infect()
                self.update_IS_events()
                self.update_ES_events()
                self.current_infected_nodes.append(exposed_infected_event)
                self.current_exposed.remove(exposed_infected_event)
                self.add_IS_events(exposed_infected_event)

            elif event_class.event_type == eventtype.EventType.RECOVER:
                recovery_event = next_event
                self.potential_recovery_events.remove_from_event_list(recovery_event)
                recovery_event.recover()
                self.update_IS_events()
                self.recovered_nodes.append(recovery_event)
                self.current_infected_nodes.remove(recovery_event)

        self.current_sim_time += tau

    def update_ES_events(self):
        updated_V_ES = []
        for edge in self.potential_ES_events.event_list:
            if (edge.get_left_node().get_state() == nodestate.NodeState.EXPOSED) \
                    and (edge.get_right_node().get_state() == nodestate.NodeState.SUSCEPTIBLE):
                updated_V_ES.append(edge)
        self.potential_ES_events.event_list = updated_V_ES

    def add_ES_events(self, infected_node):
        for j in range(0, len(self.A[infected_node.get_label()])):
            if self.A[infected_node.get_label()][j] == 1:
                candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self.Theta_ExposedInfected[j])
                if self.track_memberships:
                    candidate_node.set_membership(self.node_memberships[candidate_node.get_label()])
                neighbor_node = self.existing_node(candidate_node)
                if neighbor_node.get_state() == nodestate.NodeState.SUSCEPTIBLE:
                    edge_ij = network.Edge(infected_node, neighbor_node,
                                           self.Beta_ExposedSusceptible[infected_node.get_label()][j])
                    if not self.edge_list_contains(edge_ij):
                        self.potential_ES_events.add_to_event_list(edge_ij)
    def record_membership_states(self):
        # TODO assign a time series vector for number of groups for membership
        # ex. current_infected_group2.append(len(current_infected_nodes where membership==group2)
        for group in self.membership_groups:
            infected_in_group = list(node for node in self.current_infected_nodes if node.get_membership() == group)
            self.membership_time_series_infc[group].append(len(infected_in_group))
            exposed_in_group = list(node for node in self.current_exposed if node.get_membership() == group)
            self.membership_time_series_exp[group].append(len(exposed_in_group))
        # maybe don't worry about recovery? need to change the list of recovered nodes to be whole Node objects,
        # not just labels
        # self.real_time_srs_rec.append(len(self.recovered_nodes))

    def init_membership_state_time_series(self):
        self.membership_time_series_infc = {}
        for group in self.membership_groups:
            self.membership_time_series_infc[group] = []
            self.membership_time_series_rec[group] = []
            self.membership_time_series_exp[group] = []

    def tabulate_continuous_time(self, time_buckets=100):
        max_time = max(self.time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        exposed_time_series = np.zeros(len(time_partition))
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

        try:
            for t in range(ts_length - 1):
                real_time_val = self.time_series[t]
                real_time_exposed = self.real_time_srs_exp[t]
                # determine what index it goes in
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        exposed_time_series[idx] = real_time_exposed
        except IndexError:
            print('No values for exposed time series')

        return time_partition, infection_time_series, recover_time_series, exposed_time_series

    def tabulate_continuous_time_with_groups(self, time_buckets=100):
        max_time = max(self.time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        infection_time_series = {}
        exposed_time_series = {}
        for group in self.membership_groups:
            infection_time_series[group] = np.zeros(len(time_partition))
            exposed_time_series[group] = np.zeros(len(time_partition))
        ts_length = len(self.time_series)

        try:
            for t in range(ts_length - 1):
                for group in self.membership_groups:
                    real_time_val = self.time_series[t]
                    real_time_infected = self.membership_time_series_infc[group][t]
                    real_time_exposed = self.membership_time_series_exp[group][t]
                    for idx in range(time_buckets):
                        upper_val = time_partition[idx]
                        if real_time_val < upper_val:
                            infection_time_series[group][idx] = real_time_infected
                            exposed_time_series[group][idx] = real_time_exposed
        except IndexError:
            print('No values for infection time series')

            # TODO tbd for recovery

        return time_partition, infection_time_series, exposed_time_series


def draw_tau(event_catalog):
    list_of_events = []
    for event_group in event_catalog:
        list_of_events.append(event_group.event_list)
    list_of_events = [item for sublist in list_of_events for item in sublist]
    sum_of_rates = np.sum(event.get_event_rate() for event in list_of_events)
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
        total_combined_rate += np.sum(event.get_event_rate() for event in events)
        partition_end_markers[event_class] = total_combined_rate
    random_draw = random.uniform(0, total_combined_rate)
    for i in partition_end_markers.keys():
        if random_draw < partition_end_markers[i]:
            return event_catalog[i]  # Document: return the whole list of events in that class


def draw_event(event_list):
    possible_events = event_list.event_list
    max_rate = max(event.get_event_rate() for event in possible_events)
    accepted = False
    random_event = None
    L = len(event_list.event_list)  # Document: drawing weighted with re-sampling (with replacement)
    while not accepted:
        random_idx = np.random.randint(0, L)
        random_event = event_list.event_list[random_idx]
        accept_rate = random_event.get_event_rate() / max_rate
        random_draw = random.uniform(0, 1)
        if random_draw < accept_rate:
            accepted = True
    return random_event
