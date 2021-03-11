from epintervene.simobjects import simulation
import numpy as np


class UniversalInterventionSim(simulation.Simulation):
    def __init__(self, A, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(adj_matrix=A, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen = None
        self.beta_redux = None
        self.intervened = False
        self.time_of_intervention = max_unitless_sim_time

    def simtype(self):
        print('I am a simulation class of type universal intervention')

    def configure_intervention(self, intervention_gen, beta_redux):
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux

    def run_sim(self, with_memberships=False):
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self._initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self.highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step()

            self.total_num_timesteps += 1
            if len(self.potential_IS_events.event_list) == 0:
                break

    def intervene(self, reduce_current_edges=False):
        N = len(self.Beta[0])
        new_Beta = np.full((N, N), self.beta_redux)
        self.Beta = new_Beta
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self.potential_IS_events.event_list:
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

    def simtype(self):
        print('I am a simulation class of type Random Intervention')

    def configure_intervention(self, intervention_gen, beta_redux, proportion_reduced):
        self.intervention_gen = intervention_gen
        self.beta_redux = beta_redux
        self.proportion_reduced = proportion_reduced

    def run_sim(self, with_memberships=False):
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self._initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self.highest_gen >= self.intervention_gen:
                    self.intervene()
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step()

            self.total_num_timesteps += 1
            if len(self.potential_IS_events.event_list) == 0:
                break

    def intervene(self, reduce_current_edges=False):
        print('intervening')
        N = len(self.A[0])
        frac_of_network = self.proportion_reduced * N
        how_many = 1
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))
        random_set = np.random.randint(0, N, how_many)
        for node in random_set:
            self.Beta[node] = np.full(N, self.beta_redux)
            # TODO column as well?
            self.Beta[:, node] = np.full(N, self.beta_redux).T
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self.potential_IS_events.event_list:
                edge.set_event_rate(self.Beta[edge.get_left_node().get_label()][edge.get_right_node().get_label()])


class MultiInterventionSim(simulation.Simulation):
    def __init__(self, A, max_unitless_sim_time=1000000, membership_groups=None, node_memberships=None):
        super().__init__(adj_matrix=A, max_unitless_sim_time=max_unitless_sim_time, membership_groups=membership_groups,
                         node_memberships=node_memberships)
        self.intervention_gen_list = None
        self.beta_redux_list = None
        self.proportion_reduced_list = None
        self.intervened_status_list = []
        self.time_of_intervention_list = []
        self.next_up_intervention_entry = 0

    def simtype(self):
        print('I am a simulation class of type Multi intervention Random Intervention')

    def configure_intervention(self, intervention_gen_list, beta_redux_list, proportion_reduced_list):
        self.intervention_gen_list = intervention_gen_list
        self.beta_redux_list = beta_redux_list
        self.proportion_reduced_list = proportion_reduced_list
        for i in range(len(intervention_gen_list)):
            self.intervened_status_list.append(False)

    def run_sim(self, with_memberships=False):
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self._initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            if self.next_up_intervention_entry < len(self.intervention_gen_list):
                if not self.intervened_status_list[self.next_up_intervention_entry]:
                    if self.highest_gen >= self.intervention_gen_list[self.next_up_intervention_entry]:
                        self.intervene(self.next_up_intervention_entry)
                        self.time_of_intervention_list.append(self.current_sim_time)
                        self.intervened_status_list[self.next_up_intervention_entry] = True
                        self.next_up_intervention_entry += 1 #TODO this isn't incrementing
                # Run one step
            self._single_step()

            self.total_num_timesteps += 1
            if len(self.potential_IS_events.event_list) == 0:
                break

    def intervene(self, intervention_entry, reduce_current_edges=False):
        print('intervening')
        N = len(self.A[0])
        frac_of_network = self.proportion_reduced_list[intervention_entry] * N
        how_many = 1
        if frac_of_network > 1:
            how_many = int(np.round(frac_of_network, 0))
        # TODO this random set needs to be nodes who are not currently vaccinated, also should they be in infected and/or recovered?
        random_set = np.random.randint(0, N, how_many)
        for node in random_set:
            self.Beta[node] = np.full(N, self.beta_redux_list[intervention_entry])
            # TODO column as well?
            self.Beta[:, node] = np.full(N, self.beta_redux_list[intervention_entry]).T
        # change event rate for each existing edge pair
        if reduce_current_edges:
            for edge in self.potential_IS_events.event_list:
                edge.set_event_rate(self.Beta[edge.get_left_node().get_label()][edge.get_right_node().get_label()])


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

    def run_sim(self, with_memberships=False):
        if with_memberships: self.track_memberships = True
        if self.track_memberships:
            self.init_membership_state_time_series()
        self._initialize_patient_zero()
        while self.current_sim_time < self.total_sim_time:
            if not self.intervened:
                if self.current_sim_time > self.intervention_time:
                    self.intervene()
                    self.time_of_intervention = self.current_sim_time
                    self.intervened = True
            # Run one step
            self._single_step()

            self.total_num_timesteps += 1
            if len(self.potential_IS_events.event_list) == 0:
                break

    def intervene(self):
        # need to re-size Beta and Gamma in case the network size has changed
        self.A = self.A_modified
        super().add_infection_event_rates(self.new_Beta)
        super().add_recover_event_rates(self.new_Gamma)
        self.prune_IS_edges()
        for node in self.potential_IS_events.event_list:
            self.add_IS_events(node)
        self._update_IS_events()
        print('Modifying network')
