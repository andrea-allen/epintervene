from epintervene.simobjects import simulation
import numpy as np
from epintervene.simobjects import network
from epintervene.simobjects import nodestate
import random


class UniversalInterventionSim(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
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

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
        self.use_uniform_rate = uniform_rate
        if self.use_uniform_rate and (
                self._uniform_beta is None or self._uniform_gamma is None):
            raise Exception(
                'Since uniform_rate=True, must set all transmission and recovery parameter values via set_uniform_{param} method}')
        if with_memberships:
            self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero(label=p_zero)
        self._gens_size_over_time = []
        self._gens_size_over_time.append(list(self._current_active_gen_sizes))
        while self._current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self._highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self._current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos,
                              record_active_gen_sizes=record_active_gen_sizes)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self, reduce_current_edges=False):
        print('intervening')
        N = self._N
        if self.use_uniform_rate:
            self._uniform_beta = self.beta_redux
        else:
            new_Beta = np.full((N, N), self.beta_redux)
            self.Beta = new_Beta
        # # change event rate for each existing edge pair
        # # NO LONGER NECESSARY WITH REQUIRED UNIFORM RATE
        # if reduce_current_edges:
        #     for edge in self._potential_IS_events._event_list:
        #         if self.use_uniform_rate:
        #             edge.set_event_rate(self.beta_redux)
        #         else:
        #             edge.set_event_rate(self.Beta[edge.i.label][edge.j.label])


class RandomInterventionSim(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
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

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
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
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos,
                              record_active_gen_sizes=record_active_gen_sizes)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self):
        print('intervening')
        if self._N == 0:
            self._N = len(self._adjlist)
        vaccinated_nodes = []
        vax_labels = []
        random_set = list(np.where(np.random.uniform(0, 1, self._N) <= self.proportion_reduced)[0])
        for node_label in random_set:
            if self.use_uniform_rate:
                candidate_node = network.Node(node_label, -1, None, self._uniform_gamma)
            else:
                candidate_node = network.Node(node_label, -1, None, self._Gamma[node_label])
            if self.track_memberships:
                candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
            existing_node = self._existing_node(candidate_node)
            if existing_node.get_state() != nodestate.NodeState.VACCINATED:
                if existing_node.get_state() != nodestate.NodeState.INFECTED:
                    existing_node.vaccinate()
                    vaccinated_nodes.append(existing_node)
                    vax_labels.append(node_label)
                # self._Beta[node_label] = np.full(N, self.beta_redux)
                # self._Beta[:, node_label] = np.full(N, self.beta_redux).T
                    self._update_IS_events(recovery_event=existing_node)


class RandomRolloutSimulation(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
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

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
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
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos,
                              record_active_gen_sizes=record_active_gen_sizes)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self, intervention_entry):
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
                            if existing_node.get_state() != nodestate.NodeState.INFECTED:
                                existing_node.vaccinate()
                                vaccinated_nodes.append(existing_node)
                                vax_labels.append(node_label)
                            # self._Beta[node_label] = np.full(self._N, self.beta_redux_list[intervention_entry])
                            # self._Beta[:, node_label] = np.full(self._N, self.beta_redux_list[intervention_entry]).T
                                self._update_IS_events(recovery_event=existing_node)
                            # else:
                            #     print("I was already infected")
        return vaccinated_nodes


class TargetedInterventionSim(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen = None
        self.beta_redux = None
        self.proportion_reduced = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time
        self.use_uniform_rate = True

    def simtype(self):
        print('I am a simulation class of type Targeted Intervention')

    def configure_intervention(self, intervention_gen, beta_redux, proportion_reduced):
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux
        self.proportion_reduced = proportion_reduced

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
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
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos,
                              record_active_gen_sizes=record_active_gen_sizes)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self):
        print('intervening')
        if self._N == 0:
            self._N = len(self._adjlist)
        frac_of_network = self.proportion_reduced * self._N
        how_many = 0
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))

        degrees = list([len(neighbor) - 1 for neighbor in self._adjlist])
        degree_classes = {}
        for i in range(len(degrees)):
            degree = degrees[i]
            try:
                degree_classes[degree].append(i)
            except KeyError:
                degree_classes[degree] = [i]

        keys_decreasing = np.sort(list(degree_classes.keys()))[::-1]

        num_vaccinated = 0
        infct_hd_nodes = 0
        for degree in keys_decreasing:
            current_node_group = degree_classes[degree] # To match analytical prediction, once a degree class is bumped up to meet the quota, the whole class is vaccinated
            if num_vaccinated < how_many:
                for node in current_node_group:
                    if self.use_uniform_rate:
                        candidate_node = network.Node(node, -1, None, self._uniform_gamma)
                    else:
                        candidate_node = network.Node(node, -1, None, self._Gamma[node])
                    if self.track_memberships:
                        candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                    existing_node = self._existing_node(candidate_node)
                    if existing_node.get_state() != nodestate.NodeState.VACCINATED:
                        if existing_node.get_state() != nodestate.NodeState.INFECTED: #Do not vaccinate currently infected nodes
                            existing_node.vaccinate()
                            num_vaccinated += 1 #Only counted if high degree nodes are successfully vaccinated
                        else:
                            infct_hd_nodes += 1 #Counts high-degree currently infected nodes
                    self._update_IS_events(recovery_event=existing_node)
        # print(f"Number of already infected nodes: {infct_hd_nodes}, which is {(infct_hd_nodes/how_many)} percent of the projected vacine proportion")
        # print(f"In total, {num_vaccinated} were vaccinated, out of a planned {how_many}")

class TargetedRolloutSimulation(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen_list = None
        self.beta_redux_list = None
        self.proportion_reduced_list = None
        self.intervened_status_list = []
        self.time_of_intervention_list = []
        self.next_up_intervention_entry = 0
        self.use_uniform_rate = True

    def simtype(self):
        print('Simulation class of type Multi intervention Targeted Intervention')

    def configure_intervention(self, intervention_gen_list, beta_redux_list, proportion_reduced_list):
        self.intervention_gen_list = intervention_gen_list
        self.beta_redux_list = beta_redux_list
        self.proportion_reduced_list = proportion_reduced_list
        for i in range(len(intervention_gen_list)):
            self.intervened_status_list.append(False)

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
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
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph, viz_pos=viz_pos,
                              record_active_gen_sizes=record_active_gen_sizes)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self, intervention_entry):
        if self._N == 0:
            self._N = len(self._adjlist)
        frac_of_network = self.proportion_reduced_list[intervention_entry] * self._N
        how_many = 0
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))

        degrees = list([len(neighbor) - 1 for neighbor in self._adjlist])
        degree_classes = {}
        for i in range(len(degrees)):
            degree = degrees[i]
            try:
                degree_classes[degree].append(i)
            except KeyError:
                degree_classes[degree] = [i]

        keys_decreasing = np.sort(list(degree_classes.keys()))[::-1]

        num_vaccinated = 0
        infct_hd_nodes = 0
        for degree in keys_decreasing:
            current_node_group = degree_classes[degree] # To match analytical prediction, once a degree class is bumped up to meet the quota, the whole class is vaccinated
            if num_vaccinated < how_many:
                for node in current_node_group:
                    if self.use_uniform_rate:
                        candidate_node = network.Node(node, -1, None, self._uniform_gamma)
                    else:
                        candidate_node = network.Node(node, -1, None, self._Gamma[node])
                    if self.track_memberships:
                        candidate_node.set_membership(self._node_memberships[candidate_node.get_label()])
                    existing_node = self._existing_node(candidate_node)
                    if existing_node.get_state() != nodestate.NodeState.VACCINATED:
                        if existing_node.get_state() != nodestate.NodeState.INFECTED: #Do not vaccinate currently infected nodes
                            existing_node.vaccinate()
                            num_vaccinated += 1 #Only counted if high degree nodes are successfully vaccinated
                        else:
                            infct_hd_nodes += 1 #Counts high-degree currently infected nodes
                    self._update_IS_events(recovery_event=existing_node)

class RingInterventionSim(simulation.Simulation):
    # TODO
    def simtype(self):
        print('I am a simulation class of type Ring Intervention')


class AbsoluteTimeNetworkSwitchSim(simulation.Simulation):
    def __init__(self, N, adjlist=None, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(N=N, adj_list=adjlist, max_unitless_sim_time=max_unitless_sim_time,
                         membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.A_modified = None
        self.adjlist_modified = None
        self.intervention_time = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time
        self.new_Beta = None
        self.new_Gamma = None
        self.use_uniform_rate = True

    def configure_intervention(self, intervention_time, new_adjlist=None, new_beta=None, new_gamma=None):
        self.adjlist_modified = new_adjlist
        # self.A_modified = new_adjacency_matrix
        self.intervention_time = intervention_time
        self.new_Beta = new_beta
        self.new_Gamma = new_gamma

    def simtype(self):
        print('I am a simulation class of type Absolute Time Network Intervention')

    def run_sim(self, with_memberships=False, uniform_rate=True, wait_for_recovery=False, visualize=False,
                viz_graph=None, viz_pos=None, p_zero=None, kill_by=None, record_active_gen_sizes=False):
        self.use_uniform_rate = uniform_rate
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self._init_membership_state_time_series()
        self._initialize_patient_zero(label=p_zero)
        while self._current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self._current_sim_time > self.intervention_time:
                    self.intervene()
                    self.time_of_intervention = self._current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step(uniform_rate=uniform_rate, visualize=visualize, viz_graph=viz_graph)

            self._total_num_timesteps += 1
            if self._current_number_IS_events == 0:
                if not wait_for_recovery:
                    break
                elif len(self._recovery_events) == 0:
                    break
            if kill_by is not None and self._highest_gen >= kill_by:
                should_quit = self._check_active_gens(kill_by)
                if should_quit:
                    break

    def intervene(self):
        if self.adjlist_modified is not None:
            self.set_adjlist(self.adjlist_modified)
        if self.new_Gamma is not None and self.new_Beta is not None:
            super().set_uniform_beta(self.new_Beta)
            super().set_uniform_gamma(self.new_Gamma)
        self._reset_node_states()
        for node_key in list(self._active_node_dict.keys()):
            node = self._active_node_dict[node_key]
            if node.get_state() == nodestate.NodeState.INFECTED:
                self._add_IS_events(node)
