from epintervene.simobjects import simulation
import numpy as np
from epintervene.simobjects import network
from epintervene.simobjects import nodestate


class UniversalInterventionSim(simulation.Simulation):
    def __init__(self, N, A=None, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_matrix=A, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen = None
        self.beta_redux = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time
        self.use_uniform_rate = True

    def simtype(self):
        print('I am a simulation class of type universal intervention')

    def configure_intervention(self, intervention_gen, beta_redux):
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False):
        self.use_uniform_rate = uniform_rate
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self._highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self._current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step(uniform_rate=self.use_uniform_rate)

            self._total_num_timesteps += 1
            if len(self._potential_IS_events._event_list) == 0:
                if not wait_for_recovery:
                    break
                elif len(self._current_infected_nodes) == 0:
                    break

    def intervene(self, reduce_current_edges=False):
        print('intervening')
        N = self._N
        if self.use_uniform_rate:
            self._uniform_beta = self.beta_redux
        else:
            new_Beta = np.full((N, N), self.beta_redux)
            self.Beta = new_Beta
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self._potential_IS_events._event_list:
                if self.use_uniform_rate:
                    edge.set_event_rate(self.beta_redux)
                else:
                    edge.set_event_rate(self.Beta[edge.i.label][edge.j.label])


class RandomInterventionSim(simulation.Simulation):
    def __init__(self, A, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(adj_matrix=A, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen = None
        self.beta_redux = None
        self.proportion_reduced = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time
        self.use_uniform_rate = True

    def simtype(self):
        print('I am a simulation class of type Random Intervention')

    def configure_intervention(self, intervention_gen, beta_redux, proportion_reduced):
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux
        self.proportion_reduced = proportion_reduced

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False):
        self.use_uniform_rate = uniform_rate
        if with_memberships: self.track_memberships = True
        self.use_uniform_rate = uniform_rate
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self._highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self._current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step(uniform_rate=self.use_uniform_rate)

            self._total_num_timesteps += 1
            if len(self._potential_IS_events._event_list) == 0:
                if not wait_for_recovery:
                    break
                elif len(self._current_infected_nodes) == 0:
                    break

    def intervene(self, reduce_current_edges=False):
        print('intervening')
        N = len(self._A[0])
        frac_of_network = self.proportion_reduced * N
        how_many = 1
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))
        vaccinated_nodes = []
        vax_labels = []
        random_set = np.unique(np.random.randint(0, N, how_many))
        for node_label in random_set:
            if len(vax_labels) < how_many:
                if node_label not in vax_labels:
                    candidate_node = network.Node(node_label, -1, None, self._Gamma[node_label])
                    if self.track_memberships:
                        candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                    existing_node = self._existing_node(candidate_node)
                    if existing_node.get_state() != nodestate.NodeState.VACCINATED:
                        existing_node.vaccinate()
                        vaccinated_nodes.append(existing_node)
                        vax_labels.append(node_label)
                        self._Beta[node_label] = np.full(N, self.beta_redux)
                        self._Beta[:, node_label] = np.full(N, self.beta_redux).T
                    self._update_IS_events(recovery_event=existing_node)
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self._potential_IS_events._event_list:
                edge.set_event_rate(self._Beta[edge.get_left_node().get_label()][edge.get_right_node().get_label()])


class RandomRolloutSimulation(simulation.Simulation):
    def __init__(self, N, adjmatrix=None, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_matrix=adjmatrix, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen_list = None
        self.beta_redux_list = None
        self.proportion_reduced_list = None
        self.intervened_status_list = []
        self.time_of_intervention_list = []
        self.next_up_intervention_entry = 0
        self.use_uniform_rate = True

    def simtype(self):
        print('I am a simulation class of type Multi intervention Random Intervention')

    def configure_intervention(self, intervention_gen_list, beta_redux_list, proportion_reduced_list):
        self.intervention_gen_list = intervention_gen_list
        self.beta_redux_list = beta_redux_list
        self.proportion_reduced_list = proportion_reduced_list
        for i in range(len(intervention_gen_list)):
            self.intervened_status_list.append(False)

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False):
        self.use_uniform_rate = uniform_rate
        if with_memberships: self.track_memberships = True
        self.use_uniform_rate = uniform_rate
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            if self.next_up_intervention_entry < len(self.intervention_gen_list):
                if not self.intervened_status_list[self.next_up_intervention_entry]:
                    if self._highest_gen >= self.intervention_gen_list[self.next_up_intervention_entry]:
                        self.intervene(self.next_up_intervention_entry)
                        self.time_of_intervention_list.append(self._current_sim_time)
                        self.intervened_status_list[self.next_up_intervention_entry] = True
                        self.next_up_intervention_entry += 1
                # Run one step
            self._single_step(uniform_rate=self.use_uniform_rate)

            self._total_num_timesteps += 1
            if len(self._potential_IS_events._event_list) == 0:
                if not wait_for_recovery:
                    break
                elif len(self._current_infected_nodes) == 0:
                    break

    def intervene(self, intervention_entry, reduce_current_edges=False):
        print('intervening')
        if self._N == 0:
            self._N = len(self._A[0])
        frac_of_network = self.proportion_reduced_list[intervention_entry] * self._N
        how_many = 0
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))
        vaccinated_nodes = []
        vax_labels = []
        while len(vaccinated_nodes) < how_many:
            random_set = np.unique(np.random.randint(0, self._N, how_many))
            for node_label in random_set:
                if len(vax_labels) < how_many:
                    if node_label not in vax_labels:
                        if self.use_uniform_rate:
                            candidate_node = network.Node(node_label, -1, None, self._uniform_gamma)
                        else:
                            candidate_node = network.Node(node_label, -1, None, self._Gamma[node_label])
                        if self.track_memberships:
                            candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                        existing_node = self._existing_node(candidate_node)
                        # all the nodes become active
                        if existing_node.get_state() != nodestate.NodeState.VACCINATED:
                            # these should also be added to the active nodes
                            existing_node.vaccinate()
                            vaccinated_nodes.append(existing_node)
                            vax_labels.append(node_label)
                            # self._Beta[node_label] = np.full(self._N, self.beta_redux_list[intervention_entry])
                            # self._Beta[:, node_label] = np.full(self._N, self.beta_redux_list[intervention_entry]).T
                            self._update_IS_events(recovery_event=existing_node)
                            if reduce_current_edges:
                                for edge in self._potential_IS_events._event_list:
                                    edge.set_event_rate(self._Beta[edge.get_left_node().get_label()][edge.get_right_node().get_label()])
        return vaccinated_nodes


class TargetedInterventionSim(simulation.Simulation):
    # TODO
    def simtype(self):
        print('I am a simulation class of type Targeted Intervention')


class RingInterventionSim(simulation.Simulation):
    # TODO
    def simtype(self):
        print('I am a simulation class of type Ring Intervention')


class AbsoluteTimeNetworkSwitchSim(simulation.Simulation):
    def __init__(self, A, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(adj_matrix=A, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.A_modified = None
        self.intervention_time = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time
        self.new_Beta = None
        self.new_Gamma = None

    def configure_intervention(self, new_adjacency_matrix, intervention_time, Beta, Gamma):
        self.A_modified = new_adjacency_matrix
        self.intervention_time = intervention_time
        self.new_Beta = Beta
        self.new_Gamma = Gamma

    def simtype(self):
        print('I am a simulation class of type Absolute Time Network Intervention')

    def run_sim(self, with_memberships=False, use_uniform_rate=True):
        self.use_uniform_rate = use_uniform_rate
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero()
        while self._current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self._current_sim_time > self.intervention_time:
                    self.intervene()
                    self.time_of_intervention = self._current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step()

            self._total_num_timesteps += 1
            if len(self._potential_IS_events._event_list) == 0:
                break

    def intervene(self):
        # need to re-size Beta and Gamma in case the network size has changed
        self.A = self.A_modified
        super().add_infection_event_rates(self.new_Beta)
        super().add_recover_event_rates(self.new_Gamma)
        self._prune_IS_edges()
        for node in self._potential_IS_events._event_list:
            self._add_IS_events(node)
        self._update_IS_events()
        print('Modifying network')
