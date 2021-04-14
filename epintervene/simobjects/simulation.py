import random
import numpy as np
from epintervene.simobjects import network
from epintervene.simobjects import eventtype
from epintervene.simobjects import nodestate


class EventList:
    def __init__(self, event_type, event_list):
        self._event_list = event_list
        self._event_type = event_type

    def add_to_event_list(self, event):
        self._event_list.append(event)

    def remove_from_event_list(self, event):
        try:
            self._event_list.remove(event)
        except ValueError:
            pass

    def get_type(self):
        return self._event_type.name

    def __getattribute__(self, item):
        return object.__getattribute__(self, item)


class Simulation:
    """
    An SIR simulation object storing network information, rate parameters, and simulation run results.

    """
    def __init__(self, N, adj_matrix=None, adj_list=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        """
        Create a single instance of an SIR Simulation object.

        :param adj_matrix: Input array consisting of the adjacency matrix for the underlying network.
        :param max_unitless_sim_time: Max value the event-driven simulation will cease running at when time draw reaches this value.
        :param membership_groups: Tags for membership groups or meta-populations in underlying network model. Cqn be any type.
        :param node_memberships: Ordered list of membership group values where list index corresponds to node label in the adjacency matrix/list.
        """
        self.total_sim_time = max_unitless_sim_time
        self._current_sim_time = 0
        self._A = adj_matrix
        self._adjlist = adj_list
        self._N = N
        self._uniform_beta = None
        self._uniform_gamma = None
        self._Beta = None
        self._Gamma = None
        self._potential_IS_events = EventList(eventtype.EventType.INFECTEDSUSCEPTIBLE, [])
        self._out_degree_IS_events = {}
        self._in_degree_IS_events = {}
        self._potential_recovery_events = EventList(eventtype.EventType.RECOVER, [])
        self._recovered_nodes = []
        self._current_infected_nodes = []
        self._highest_gen = 0
        self._gen_collection = {}
        self._active_nodes = []
        self._active_node_dict = {}
        self._total_num_timesteps = 0
        self._time_series = [0]
        self._real_time_srs_infc = []
        self._real_time_srs_rec = []
        self._generational_emergence = {0: 0}
        self._membership_groups = membership_groups
        self._membership_time_series_infc = {0: []}
        self._membership_time_series_rec = {0: []}
        self.track_memberships = False
        self.use_uniform_rate = True
        self._node_memberships = node_memberships


    def get_Beta(self):
        return self._Beta

    def get_Gamma(self):
        return self._Gamma

    def get_gen_collection(self):
        return self._gen_collection

    def get_generational_emergence(self):
        return self._generational_emergence

    def get_time_series_vals(self):
        return self._time_series

    def get_real_time_series_infected(self):
        return self._real_time_srs_infc

    def get_real_time_series_recovered(self):
        return self._real_time_srs_rec

    def get_membership_time_series_infected(self):
        return self._membership_time_series_infc

    def get_membership_time_series_recovered(self):
        return self._membership_time_series_rec

    def set_adjlist(self, adjlist):
        self._adjlist = adjlist

    def set_uniform_beta(self, beta):
        self._uniform_beta = beta

    def set_uniform_gamma(self, gamma):
        self._uniform_gamma = gamma

    def _initialize_patient_zero(self):
        if self._N==0:
            self._N = len(self._A[0])
        p_zero_idx = random.randint(0, self._N - 1)
        if self._Gamma is None and self._uniform_gamma is None:
            raise AttributeError(self._Gamma,
                                 'Please provide a recovery vector Gamma to the simulation via the '
                                 'add_recover_event_rates method, or a uniform gamma with the set_uniform_gamma method')
        if self._Beta is None and self._uniform_beta is None:
            raise AttributeError(self._Beta,
                                 'Please provide an infection rate matrix Beta to the simulation via the '
                                 'add_infection_event_rates method, or set a uniform beta with the set_uniform_beta method')
        if self.use_uniform_rate:
            patient_zero = network.Node(p_zero_idx, 0, nodestate.NodeState.INFECTED, event_rate=self._uniform_gamma)
        else:
            patient_zero = network.Node(p_zero_idx, 0, nodestate.NodeState.INFECTED, event_rate=self._Gamma[p_zero_idx])
        if self.track_memberships:
            patient_zero.set_membership(self._node_memberships[p_zero_idx])
        self._active_nodes.append(patient_zero)
        self._active_node_dict[p_zero_idx] = patient_zero
        self._current_infected_nodes.append(patient_zero)
        self._gen_collection[0] = [p_zero_idx]
        self._potential_recovery_events.add_to_event_list(patient_zero)
        self._out_degree_IS_events[p_zero_idx] = []
        self._add_IS_events(patient_zero)

    def add_infection_event_rates(self, Beta):
        """
        Add matrix that is same dimensions of adj_matrix where each entry is transmission rate
        between nodes

        :param Beta: Matrix where each entry (i,j) contains float between 0 and 1 encoding probability of single transmission event from infected node (i) to susceptible node (j)
        """
        if len(Beta) != self._N:
            raise ValueError(Beta, 'Beta matrix must be N by N to match network size')
        else:
            self._Beta = Beta  # Document: Beta must be an N x N matrix with entries for the probability of infection
            # between each pair of nodes
            for edge in self._potential_IS_events._event_list:
                edge.set_event_rate(self._Beta[edge.get_left_node().get_label(), edge.get_right_node().get_label()])

    def add_recover_event_rates(self, Gamma):
        """
        Add a vector of length the size of the adj_matrix where each entry is the recovery rate for each node

        :param Gamma: Vector where each entry [i] is the recovery rate for node i
        """
        if len(Gamma) != self._N:
            raise ValueError(Gamma, 'Gamma vector must be length N to match network size')
        else:
            self._Gamma = Gamma  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self._potential_recovery_events._event_list:
                node.set_event_rate(self._Gamma[node.get_label()])

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False, viz_graph=None, viz_pos=None):
        """
        Main method for running a single realization of the epidemic simulation.

        :param with_memberships: Set to True if specified memberships in Simulation configuration
        :param uniform_rate: Will use a uniform transmission and recovery event probabilites instead of a customized matrix. Increases simulation speed, use if possible.
        :param wait_for_recovery: Will run simulation until all possible nodes have recovered.
        """
        self.use_uniform_rate = uniform_rate
        if self.use_uniform_rate and (
                self._uniform_beta is None or self._uniform_gamma is None):
            raise Exception(
                'Since uniform_rate=True, must set all transmission and recovery parameter values via set_uniform_{param} method}')
        if with_memberships:
            self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            # Run one step
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos)

            self._total_num_timesteps += 1
            if len(self._potential_IS_events._event_list) == 0:
                if not wait_for_recovery:
                    break
                elif len(self._current_infected_nodes) == 0:
                    break

    def _single_step(self, visualize=False, uniform_rate=True, viz_graph=None, viz_pos=None):
        self.use_uniform_rate = uniform_rate
        if visualize:
            self._visualize_network(viz_graph, viz_pos)
        event_catolog = [self._potential_IS_events, self._potential_recovery_events]
        # note: to optimize, going to stop updating IS edges before this step. it will throw off tau by a tiny bit, but usually only by 1 or 2 events and shouldn't affect the whole distribution
        tau = draw_tau(event_catolog, uniform_rate=self.use_uniform_rate)
        self._current_sim_time += tau

        self._time_series.append(self._time_series[-1] + tau)
        self._real_time_srs_infc.append(len(self._current_infected_nodes))
        self._real_time_srs_rec.append(len(self._recovered_nodes))
        if self.track_memberships:
            self._record_membership_states()
        event_class = draw_event_class(event_catolog, uniform_rate=self.use_uniform_rate)
        if event_class is not None:
            next_event = draw_event(event_class, self.use_uniform_rate)
            if event_class._event_type == eventtype.EventType.INFECTEDSUSCEPTIBLE:
                infection_event = next_event
                infection_event.infect()
                self._potential_IS_events.remove_from_event_list(infection_event)
                self._potential_recovery_events.add_to_event_list(infection_event.get_right_node())
                self._current_infected_nodes.append(infection_event.get_right_node())

                try:
                    self._gen_collection[infection_event.get_right_node().get_generation()].append(
                        infection_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self._gen_collection[infection_event.get_right_node().get_generation()] = [
                        infection_event.get_right_node().get_label()]
                    self._highest_gen += 1
                    self._generational_emergence[self._highest_gen] = self._current_sim_time
                self._update_IS_events(infection_IS_event=infection_event)
                self._add_IS_events(infection_event.get_right_node())
            if event_class._event_type == eventtype.EventType.RECOVER:
                recovery_event = next_event
                self._potential_recovery_events.remove_from_event_list(recovery_event)
                try:
                    self._current_infected_nodes.remove(recovery_event)
                except:
                    print(recovery_event.get_label())
                recovery_event.recover()
                self._update_IS_events(recovery_event=recovery_event)
                self._recovered_nodes.append(recovery_event)

    def _update_IS_events(self, infection_IS_event=None, recovery_event=None):
        if infection_IS_event is not None:
            try:
                in_degree_events = self._in_degree_IS_events[infection_IS_event.get_right_node().get_label()]
                for event in in_degree_events:
                    self._potential_IS_events.remove_from_event_list(event)
                    left_node_idx = event.get_left_node().get_label()
                    try:
                        self._out_degree_IS_events[left_node_idx].remove(event)
                    except ValueError:
                        pass
                self._in_degree_IS_events[infection_IS_event.get_right_node().get_label()] = []
            except KeyError:
                pass
        elif recovery_event is not None:
            try:
                out_degree_events = self._out_degree_IS_events[recovery_event.get_label()]
                for event in out_degree_events:
                    self._potential_IS_events.remove_from_event_list(event)
                self._out_degree_IS_events[recovery_event.get_label()] = []
            except KeyError:
                pass
            try:
                # note this will only do anything during vaccination, removes any edges pointing to the now
                #  vaccinated node
                in_degree_events = self._in_degree_IS_events[recovery_event.get_label()]
                for event in in_degree_events:
                    self._potential_IS_events.remove_from_event_list(event)
                self._in_degree_IS_events[recovery_event.get_label()] = []
            except KeyError:
                pass

    def _add_IS_events(self, infected_node):
        infection_adjlist = self._adjlist[infected_node.get_label()]
        if len(infection_adjlist)>1:
            infected_label = infected_node.get_label()
            length = len(infection_adjlist)
            for n in range(1, length):
                j = infection_adjlist[n]
                if self.use_uniform_rate:
                    candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self._uniform_gamma)
                else:
                    candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self._Gamma[j])
                if self.track_memberships:
                    candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                neighbor_node = self._existing_node(candidate_node)
                if neighbor_node.get_state() == nodestate.NodeState.SUSCEPTIBLE:
                    if self.use_uniform_rate:
                        edge_ij = network.Edge(infected_node, neighbor_node, self._uniform_beta)
                    else:
                        edge_ij = network.Edge(infected_node, neighbor_node, self._Beta[infected_label][j])
                    self._potential_IS_events.add_to_event_list(edge_ij)
                    try:
                        self._out_degree_IS_events[infected_node.get_label()].append(edge_ij)
                    except KeyError:
                        self._out_degree_IS_events[infected_node.get_label()] = [edge_ij]
                    try:
                        self._in_degree_IS_events[neighbor_node.get_label()].append(edge_ij)
                    except KeyError:
                        self._in_degree_IS_events[neighbor_node.get_label()] = [edge_ij]

    def _existing_node(self, candidate_node):
        if candidate_node.get_label() in self._active_node_dict.keys():
            return self._active_node_dict[candidate_node.get_label()]
        else:
            self._active_node_dict[candidate_node.get_label()] = candidate_node
        return candidate_node

    def _edge_list_contains(self, edge):
        for e in self._potential_IS_events._event_list:
            if e.equals(edge):
                return True
        return False

    def _prune_IS_edges(self):
        #TODO need to work on to use only adjlist, new adjlist
        for edge in self._potential_IS_events._event_list:
            try:
                edge_exists_in_network = (
                        self._A[edge.get_left_node().get_label()][edge.get_right_node().get_label()] == 1)
                if not edge_exists_in_network:
                    self._potential_IS_events.remove_from_event_list(edge)
            except IndexError:
                self._potential_IS_events.remove_from_event_list(
                    edge)  # This happens if the new network no longer contains that node, can remove them

    def _record_membership_states(self):
        for group in self._membership_groups:
            infected_in_group = list(node for node in self._current_infected_nodes if node.get_membership() == group)
            self._membership_time_series_infc[group].append(len(infected_in_group))

    def _init_membership_state_time_series(self):
        self._membership_time_series_infc = {}
        for group in self._membership_groups:
            self._membership_time_series_infc[group] = []
            self._membership_time_series_rec[group] = []

    def _visualize_network(self, viz_graph, viz_pos):
        network.visualize(self._N, viz_graph, viz_pos, self._gen_collection)
        print('In progress')

    def tabulate_generation_results(self, max_gens):
        """
        Compile the results from the simulation in terms of epidemic generations.

        :param max_gens: Maximum generations desired to tabulate and return information for
        :return: Array indexed by generations from 0 to max_gens of how many cumulative infected nodes there are in generations less than or equal to the index
        """
        gens = max_gens
        infection_time_srs_by_gen = np.zeros(gens)
        s = 1
        m = 1
        s_max = 1
        for gen in range(max_gens):  # {0: 1, 1: 12, 14, 2: 16, 42, ....
            infection_time_srs_by_gen[gen] = s
            try:
                m = len(self._gen_collection[gen + 1])  # num infected in gen g
                s += m
                s_max = s
            except KeyError:
                # make m=0 and s=the last s for the rest of the "time series"
                s = s_max
                m = 0
        return infection_time_srs_by_gen

    def tabulate_continuous_time(self, time_buckets=100, custom_range=False, custom_t_lim=100):
        """
        Compile the results from the simulation as a time series.

        :param time_buckets: How many ticks to partition the total time into
        :param custom_range: Defaults to False, set to True if desire a set maximum time value (suggested if collecting an ensemble of results)
        :param custom_t_lim: If custom_range is True, provide max time value for the time series
        :return: three arrays: time series values, number of infections time series, number recovered time series
        """
        max_time = max(self._time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        if custom_range:
            time_partition = np.linspace(0, custom_t_lim, time_buckets, endpoint=False)
        infection_time_series = np.zeros(len(time_partition))
        recover_time_series = np.zeros(len(time_partition))
        ts_length = len(self._time_series)

        try:
            for t in range(ts_length - 1):
                real_time_val = self._time_series[t]
                real_time_infected = self._real_time_srs_infc[t]
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        infection_time_series[idx] = real_time_infected
        except IndexError:
            print('No values for infection time series')

        try:
            for t in range(ts_length - 1):
                real_time_val = self._time_series[t]
                real_time_recovered = self._real_time_srs_rec[t]
                # determine what index it goes in
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        recover_time_series[idx] = real_time_recovered
        except IndexError:
            print('No values for recovery time series')

        return time_partition, infection_time_series, recover_time_series

    def tabulate_continuous_time_with_groups(self, time_buckets=100, custom_range=False, custom_t_lim=100):
        """
        Compile the results of the simulation as a dictionary of time series, yielding a distinct time series for
        each membership group.

        :param time_buckets: How many ticks to partition the total time into
        :param custom_range: Defaults to False, set to True if desire a set maximum time value (suggested if collecting an ensemble of results)
        :param custom_t_lim: If custom_range is True, provide max time value for the time series
        :return: two objects: the time series values, and a dict of infection time series keyed by membership group
        """
        max_time = max(self._time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        if custom_range:
            time_partition = np.linspace(0, custom_t_lim, time_buckets, endpoint=False)
        infection_time_series = {}
        for group in self._membership_groups:
            infection_time_series[group] = np.zeros(len(time_partition))
        # TODO tbd grouped recovery
        # recover_time_series = np.zeros(len(time_partition))
        ts_length = len(self._time_series)

        try:
            for t in range(ts_length - 1):
                for group in self._membership_groups:
                    real_time_val = self._time_series[t]
                    real_time_infected = self._membership_time_series_infc[group][t]
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
    """
    An SEIR simulation object storing network information, rate parameters, and simulation run results.

    """
    def __init__(self, N, adjmatrix=None, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        """
        Create a single instance of an SEIR Simulation object.

        :param N: Number of nodes in the network
        :param adjmatrix: Adjacency matrix, numpy array, can be None
        :param adjlist: Adjacency list, must be symmetric. Nested list should have entries of form e.g. [i j k] where node i has neighbors j and k. New source node is always the left-most node of a new line.
        :param max_unitless_sim_time: Max value the event-driven simulation will cease running at when time draw reaches this value.
        :param membership_groups:  Tags for membership groups or meta-populations in underlying network model. Cqn be any type.
        :param node_memberships: Ordered list of membership group values where list index corresponds to node label in the adjacency matrix/list.
        """
        super().__init__(adj_matrix=adjmatrix, adj_list=adjlist, N=N, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups, node_memberships=node_memberships)
        self._potential_ES_events = EventList(eventtype.EventType.EXPOSEDSUSCEPTIBLE, [])
        self._potential_EI_events = EventList(eventtype.EventType.EXPOSEDINFECTED, [])
        self._current_exposed = []
        self._real_time_srs_exp = []
        self._Beta_ExposedSusceptible = None
        self._Theta_ExposedInfected = None
        self._cont_time_exposed_dict = {}
        self._membership_time_series_exp = {0: []}
        self._in_degree_ES_events = {}
        self._out_degree_ES_events = {}
        self._uniform_gamma_ei = None
        self._uniform_beta_es = None
        # TODO deal with generations of exposed emergence as well
        # should probably also determine in the model if there's completely asymptomatic cases as well (E->R events)

    def set_uniform_beta_es(self, beta):
        self._uniform_beta_es = beta

    def set_uniform_gamma_ei(self, gamma):
        self._uniform_gamma_ei = gamma

    def add_exposed_event_rates(self, Beta_Exp):
        if len(Beta_Exp) != self._N:
            raise ValueError(Beta_Exp, 'Beta_Exp matrix must be size N by N to match network size')
        else:
            self._Beta_ExposedSusceptible = Beta_Exp  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self._potential_ES_events._event_list:
                node.set_event_rate(self._Beta_ExposedSusceptible[node.get_label()])

    def add_exposed_infected_event_rates(self, Theta):
        if len(Theta) != self._N:
            raise ValueError(Theta, 'Theta ExposedInfectious vector must be length N to match network size')
        else:
            self._Theta_ExposedInfected = Theta  # Document: Beta must be an N x N matrix with entries for the probability of
            # infection between each pair of nodes
            for node in self._potential_ES_events._event_list:
                node.set_event_rate(self._Beta_ExposedSusceptible[node.get_label()])

    # TODO number recovered is going over 100%, ugh.
    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False):
        """
        Main method for running a single realization of the epidemic simulation.

        :param with_memberships: Set to True if specified memberships in Simulation configuration.
        :param uniform_rate: Will use a uniform transmission and recovery event probabilites instead of a customized matrix. Increases simulation speed, use if possible.
        :param wait_for_recovery: Will run simulation until all possible nodes have recovered.
        """
        self.use_uniform_rate = uniform_rate
        if self.use_uniform_rate and (self._uniform_beta_es is None or self._uniform_gamma_ei is None or self._uniform_gamma is None or self._uniform_beta is None):
            raise Exception('Since uniform_rate=True, must set all transmission and recovery parameter values via set_uniform_{param} method}')
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            # Run one step
            self._single_step(uniform_rate=uniform_rate)

            self._total_num_timesteps += 1
            if (len(self._potential_IS_events._event_list) == 0) \
                    and (len(self._potential_ES_events._event_list) == 0) \
                    and (len(self._potential_EI_events._event_list) == 0):
                if not wait_for_recovery:
                    break
                elif len(self._current_infected_nodes) == 0:
                    break

    def _single_step(self, visualize=False, uniform_rate=True):
        self.use_uniform_rate = uniform_rate
        if visualize:
            self._visualize_network()
        event_catolog = [self._potential_IS_events, self._potential_recovery_events, self._potential_ES_events,
                         self._potential_EI_events]
        tau = draw_tau(event_catolog, uniform_rate=self.use_uniform_rate)

        self._time_series.append(self._time_series[-1] + tau)
        self._real_time_srs_infc.append(len(self._current_infected_nodes))
        self._real_time_srs_rec.append(len(self._recovered_nodes))
        self._real_time_srs_exp.append(len(self._current_exposed))
        if self.track_memberships:
            self._record_membership_states()

        event_class = draw_event_class(event_catolog, uniform_rate=self.use_uniform_rate)  # This is returning a whole list of events
        if event_class is not None:
            next_event = draw_event(event_class, self.use_uniform_rate)
            if event_class._event_type == eventtype.EventType.INFECTEDSUSCEPTIBLE:
                infection_event = next_event
                infection_event.expose()
                self._potential_IS_events.remove_from_event_list(infection_event)
                self._current_exposed.append(infection_event.get_right_node())
                if self.use_uniform_rate:
                    infection_event.get_right_node().set_event_rate(
                        self._uniform_gamma_ei)
                else:
                    infection_event.get_right_node().set_event_rate(
                        self._Theta_ExposedInfected[infection_event.get_right_node().get_label()])
                self._potential_EI_events.add_to_event_list(infection_event.get_right_node())
                try:
                    self._gen_collection[infection_event.get_right_node().get_generation()].append(
                        infection_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self._gen_collection[infection_event.get_right_node().get_generation()] = [
                        infection_event.get_right_node().get_label()]
                    self._highest_gen += 1
                    self._generational_emergence[self._highest_gen] = self._current_sim_time
                self._update_IS_events(infection_IS_event = infection_event)
                self._update_ES_events(infection_ES_event = infection_event)
                self._add_ES_events(infection_event.get_right_node())
            elif event_class._event_type == eventtype.EventType.EXPOSEDSUSCEPTIBLE:
                exposure_event = next_event
                exposure_event.expose()
                self._potential_ES_events.remove_from_event_list(exposure_event)
                self._current_exposed.append(exposure_event.get_right_node())
                if self.use_uniform_rate:
                    exposure_event.get_right_node().set_event_rate(
                        self._uniform_gamma_ei)
                else:
                    exposure_event.get_right_node().set_event_rate(
                        self._Theta_ExposedInfected[exposure_event.get_right_node().get_label()])
                self._potential_EI_events.add_to_event_list(
                    exposure_event.get_right_node())
                try:
                    self._gen_collection[exposure_event.get_right_node().get_generation()].append(
                        exposure_event.get_right_node().get_label())
                except KeyError:  # Need a better way than KeyError to catch a new generation
                    self._gen_collection[exposure_event.get_right_node().get_generation()] = [
                        exposure_event.get_right_node().get_label()]
                    self._highest_gen += 1
                    self._generational_emergence[self._highest_gen] = self._current_sim_time
                self._update_IS_events(infection_IS_event=exposure_event)
                self._update_ES_events(infection_ES_event=exposure_event)
                self._add_ES_events(exposure_event.get_right_node())

            elif event_class._event_type == eventtype.EventType.EXPOSEDINFECTED:
                exposed_infected_event = next_event
                self._potential_EI_events.remove_from_event_list(exposed_infected_event)
                exposed_infected_event.infect()
                if self.use_uniform_rate:
                    exposed_infected_event.set_event_rate(
                        self._uniform_gamma)
                else:
                    exposed_infected_event.set_event_rate(
                        self._Gamma[exposed_infected_event.get_label()])
                self._potential_recovery_events.add_to_event_list(exposed_infected_event)
                self._update_ES_events(exposed_infected_event=exposed_infected_event)
                self._current_infected_nodes.append(exposed_infected_event)
                self._current_exposed.remove(exposed_infected_event)
                self._add_IS_events(exposed_infected_event)

            elif event_class._event_type == eventtype.EventType.RECOVER:
                recovery_event = next_event
                self._potential_recovery_events.remove_from_event_list(recovery_event)
                recovery_event.recover()
                self._update_IS_events(recovery_event=recovery_event)
                self._recovered_nodes.append(recovery_event)
                self._current_infected_nodes.remove(recovery_event)

        self._current_sim_time += tau

    def _update_ES_events(self, infection_ES_event=None, exposed_infected_event=None):
        if infection_ES_event is not None:
            try:
                in_degree_events = self._in_degree_ES_events[infection_ES_event.get_right_node().get_label()]
                for event in in_degree_events:
                    self._potential_ES_events.remove_from_event_list(event)
                    left_node_idx = event.get_left_node().get_label()
                    try:
                        self._out_degree_ES_events[left_node_idx].remove(event)
                    except ValueError:
                        pass
                self._in_degree_ES_events[infection_ES_event.get_right_node().get_label()] = []
            except KeyError:
                pass
        elif exposed_infected_event is not None:
            try:
                out_degree_events = self._out_degree_ES_events[exposed_infected_event.get_label()]
                for event in out_degree_events:
                    self._potential_ES_events.remove_from_event_list(event)
                self._out_degree_ES_events[exposed_infected_event.get_label()] = []
            except KeyError:
                pass
            try:
                in_degree_events = self._in_degree_ES_events[exposed_infected_event.get_label()]
                for event in in_degree_events:
                    self._potential_ES_events.remove_from_event_list(event)
                self._in_degree_ES_events[exposed_infected_event.get_label()] = []
            except KeyError:
                pass

    def _add_ES_events(self, infected_node):
        infection_adjlist = self._adjlist[infected_node.get_label()]
        if len(infection_adjlist)>1:
            infected_label = infected_node.get_label()
            length = len(infection_adjlist)
            for n in range(1, length):
                j = infection_adjlist[n]
                if self.use_uniform_rate:
                    candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self._uniform_gamma_ei)
                else:
                    candidate_node = network.Node(j, -1, nodestate.NodeState.SUSCEPTIBLE, self._Theta_ExposedInfected[j])
                if self.track_memberships:
                    candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                neighbor_node = self._existing_node(candidate_node) # 18 seconds
                if neighbor_node.get_state() == nodestate.NodeState.SUSCEPTIBLE:
                    if self.use_uniform_rate:
                        edge_ij = network.Edge(infected_node, neighbor_node, self._uniform_beta)
                    else:
                        edge_ij = network.Edge(infected_node, neighbor_node, self._Beta_ExposedSusceptible[infected_node.get_label()][j])
                    self._potential_ES_events.add_to_event_list(edge_ij)
                    try:
                        self._out_degree_ES_events[infected_node.get_label()].append(edge_ij)
                    except KeyError:
                        self._out_degree_ES_events[infected_node.get_label()] = [edge_ij]
                    try:
                        self._in_degree_ES_events[neighbor_node.get_label()].append(edge_ij)
                    except KeyError:
                        self._in_degree_ES_events[neighbor_node.get_label()] = [edge_ij]

    def _record_membership_states(self):
        # TODO assign a time series vector for number of groups for membership
        # ex. current_infected_group2.append(len(current_infected_nodes where membership==group2)
        for group in self._membership_groups:
            infected_in_group = list(node for node in self._current_infected_nodes if node.get_membership() == group)
            self.membership_time_series_infc[group].append(len(infected_in_group))
            exposed_in_group = list(node for node in self._current_exposed if node.get_membership() == group)
            self._membership_time_series_exp[group].append(len(exposed_in_group))
        # maybe don't worry about recovery? need to change the list of recovered nodes to be whole Node objects,
        # not just labels
        # self.real_time_srs_rec.append(len(self.recovered_nodes))

    def _init_membership_state_time_series(self):
        self.membership_time_series_infc = {}
        for group in self._membership_groups:
            self.membership_time_series_infc[group] = []
            self._membership_time_series_rec[group] = []
            self._membership_time_series_exp[group] = []

    def tabulate_continuous_time(self, time_buckets=100, custom_range=False, custom_t_lim=100):
        # TODO add doc strings since return type differs
        max_time = max(self._time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        if custom_range:
            time_partition = np.linspace(0, custom_t_lim, time_buckets, endpoint=False)
        exposed_time_series = np.zeros(len(time_partition))
        infection_time_series = np.zeros(len(time_partition))
        recover_time_series = np.zeros(len(time_partition))
        ts_length = len(self._time_series)

        try:
            for t in range(ts_length - 1):
                real_time_val = self._time_series[t]
                real_time_infected = self._real_time_srs_infc[t]
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        infection_time_series[idx] = real_time_infected
        except IndexError:
            print('No values for infection time series')

        try:
            for t in range(ts_length - 1):
                real_time_val = self._time_series[t]
                real_time_recovered = self._real_time_srs_rec[t]
                # determine what index it goes in
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        recover_time_series[idx] = real_time_recovered
        except IndexError:
            print('No values for recovery time series')

        try:
            for t in range(ts_length - 1):
                real_time_val = self._time_series[t]
                real_time_exposed = self._real_time_srs_exp[t]
                # determine what index it goes in
                for idx in range(time_buckets):
                    upper_val = time_partition[idx]
                    if real_time_val < upper_val:
                        exposed_time_series[idx] = real_time_exposed
        except IndexError:
            print('No values for exposed time series')

        return time_partition, infection_time_series, recover_time_series, exposed_time_series

    def tabulate_continuous_time_with_groups(self, time_buckets=100, custom_range=False, custom_t_lim=100):
        max_time = max(self._time_series)
        time_partition = np.linspace(0, max_time + 1, time_buckets, endpoint=False)
        if custom_range:
            time_partition = np.linspace(0, custom_t_lim, time_buckets, endpoint=False)
        infection_time_series = {}
        exposed_time_series = {}
        for group in self._membership_groups:
            infection_time_series[group] = np.zeros(len(time_partition))
            exposed_time_series[group] = np.zeros(len(time_partition))
        ts_length = len(self._time_series)

        try:
            for t in range(ts_length - 1):
                for group in self._membership_groups:
                    real_time_val = self._time_series[t]
                    real_time_infected = self.membership_time_series_infc[group][t]
                    real_time_exposed = self._membership_time_series_exp[group][t]
                    for idx in range(time_buckets):
                        upper_val = time_partition[idx]
                        if real_time_val < upper_val:
                            infection_time_series[group][idx] = real_time_infected
                            exposed_time_series[group][idx] = real_time_exposed
        except IndexError:
            print('No values for infection time series')

            # TODO tbd for recovery

        return time_partition, infection_time_series, exposed_time_series


def draw_tau(event_catalog, uniform_rate=False):
    list_of_events = []
    for event_group in event_catalog:
        list_of_events.append(event_group._event_list)
    if uniform_rate:
        sum_of_rates = 0
        for i in range(len(list_of_events)):
            current_list = list_of_events[i]
            try:
                rate = current_list[0].get_event_rate()
                sum_of_rates += rate * len(current_list)
            except IndexError:
                continue

    else:
        list_of_events = [item for sublist in list_of_events for item in sublist]
        sum_of_rates = np.sum(event.get_event_rate() for event in list_of_events)
    if sum_of_rates == 0:
        sum_of_rates += .0001
    tau = np.random.exponential(1 / sum_of_rates)
    return tau


def draw_event_class(event_catalog, uniform_rate=False):
    num_event_classes = len(event_catalog)
    partition_end_markers = {}
    total_combined_rate = 0
    for event_class in range(num_event_classes):
        events = event_catalog[event_class]._event_list
        if uniform_rate:
            if len(events) > 0:
                idx = 0
                first_event_rate = events[idx].get_event_rate()
                if first_event_rate==0.0:
                    print(first_event_rate)
                    while first_event_rate==0:
                        idx += 1
                        first_event_rate = events[idx].get_event_rate()
                        print(first_event_rate)
                total_combined_rate += events[idx].get_event_rate() * len(events)
            partition_end_markers[event_class] = total_combined_rate
        else:
            total_combined_rate += np.sum(event.get_event_rate() for event in events)
            partition_end_markers[event_class] = total_combined_rate
    random_draw = random.uniform(0, total_combined_rate)
    for i in partition_end_markers.keys():
        if random_draw < partition_end_markers[i]:
            return event_catalog[i]  # Document: return the whole list of events in that class


def draw_event(event_list, use_uniform_rate=False):
    possible_events = event_list._event_list
    if use_uniform_rate:
        max_rate = possible_events[0].get_event_rate()
    else:
        max_rate = max(event.get_event_rate() for event in possible_events)
    accepted = False
    random_event = None
    L = len(event_list._event_list)  # Document: drawing weighted with re-sampling (with replacement)
    while not accepted:
        random_idx = np.random.randint(0, L)
        random_event = event_list._event_list[random_idx]
        if not use_uniform_rate:
            accept_rate = random_event.get_event_rate() / max_rate
            random_draw = random.uniform(0, 1)
            if random_draw < accept_rate:
                accepted = True
        else:
            return random_event
    return random_event
