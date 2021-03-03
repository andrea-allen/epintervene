import time
import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
import networkx as nx
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
import math

def random_vaccination():
    print('Manually testing random vaccination')
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb()
    G, pos = nb.from_degree_distribution(10000, powerlaw)
    A = np.array(nx.adjacency_matrix(G).todense())
    # TODO need to make sure this isn't messing with the simulation run time
    sim = extended_simulation.RandomInterventionSim(A)
    Beta = np.full((len(A), len(A)), 0.9)
    Gamma = np.full(len(A), 0.001)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    sim.configure_intervention(intervention_gen=4, beta_redux=0.0, proportion_reduced=0.2)
    sim.run_sim()

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.show()

def sbm_membership():
    G, pos, A = create_stochastic_block_model()
    # SIR model sandbox, want to have different curves for different group memberships

    staff_pop = int(len(A)/3)
    res_pop = int(len(A)/3)
    county_pop = int(len(A)/3)

    node_membership_vector = []
    for i in range(staff_pop):
        node_membership_vector.append('staff')
    for j in range(res_pop):
        node_membership_vector.append('residents')
    for k in range(county_pop):
        node_membership_vector.append('county')
    sim = simulation.Simulation(A, membership_groups=['staff', 'residents', 'county'], node_memberships=node_membership_vector)
    Beta = np.full((len(A), len(A)), 0.5)
    Gamma = np.full(len(A), 0.9)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    sim.run_sim(with_memberships=True)

    ts, membership_ts_infc = sim.tabulate_continuous_time_with_groups(time_buckets=1000)
    plt.figure(0)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group)
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.show()

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(time_buckets=1000)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.show()


    # SEIR model sandbox
    # todo add membership feature for infected and exposed to SEIR model
    sim = simulation.SimulationSEIR(A, node_memberships=node_membership_vector, membership_groups=['staff', 'residents', 'county'])
    Beta_IS = np.full((len(A), len(A)), 0.25)
    Gamma = np.full(len(A), 1.0)
    Beta_ES = np.full((len(A), len(A)), 0.25)
    Theta_EI = np.full(len(A), 1.0)
    sim.add_infection_event_rates(Beta_IS)
    sim.add_exposed_event_rates(Beta_ES)
    sim.add_recover_event_rates(Gamma)
    sim.add_exposed_infected_event_rates(Theta_EI)

    sim.run_sim(with_memberships=True)

    ts, membership_ts_infc, membership_ts_exp = sim.tabulate_continuous_time_with_groups(1000)
    plt.figure(0)
    colors = {'staff':'blue', 'residents':'red', 'county':'orange'}
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group, color=colors[group])
        plt.plot(ts, membership_ts_exp[group], label=group, ls='--', color=colors[group])
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.show()



    ts, infect_ts, recover_ts, exposed_ts = sim.tabulate_continuous_time(1000)

    plt.figure(1)
    plt.plot(ts, infect_ts, color='red', label='Infected')
    plt.plot(ts, exposed_ts, color='orange', label='Exposed')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(max_gens=20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.show()

def run():
    print('Running experiments in sandbox')
    # Stochastic block model code
    A = np.zeros((60, 60))
    block_1_p = 0.5
    block_2_p = 0.15
    block_3_p = 0.15
    block_12_p = 0.02
    block_13_p = 0.8
    block_23_p = 0.15
    N = len(A[0])
    mid_point = int(N/2)
    for i in range(0, mid_point):
        for j in range(i+1, mid_point):
            A[i,j] = 1 if np.random.rand() < block_1_p else 0
            A[j,i] = A[i,j]
        for j in range(mid_point, N):
            A[i,j] = 1 if np.random.rand() < block_12_p else 0
            A[j,i] = A[i,j]
    for i in range(mid_point, N):
        for j in range(i+1, N):
            A[i,j] = 1 if np.random.rand() < block_2_p else 0
            A[j,i] = A[i,j]


    G, pos = network.NetworkBuilder.from_adjacency_matrix(A, return_pos=True)
    node_colors = {}
    for n in G.nodes():
        if n < mid_point:
            node_colors[n] = 'blue'
        else:
            node_colors[n] = 'red'
    nx.draw_networkx_nodes(G, pos, G.nodes(), node_color=node_colors.values(), node_size=20)
    nx.draw_networkx_edges(G, pos)
    plt.show()

    A = np.random.random_integers(0, 1, (100, 100))
    A = (A + A.T) / 2
    np.fill_diagonal(A, 0)

    # SIR model sandbox

    sim = simulation.Simulation(A)
    Beta = np.full((len(A), len(A)), 0.5)
    Gamma = np.full(len(A), 0.9)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    sim.run_sim()
    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)

    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.show()

    # SEIR model sandbox
    sim = simulation.SimulationSEIR(A)
    Beta_IS = np.full((len(A), len(A)), 0.25)
    Gamma = np.full(len(A), 1.0)
    Beta_ES = np.full((len(A), len(A)), 0.25)
    Theta_EI = np.full(len(A), 1.0)
    sim.add_infection_event_rates(Beta_IS)
    sim.add_exposed_event_rates(Beta_ES)
    sim.add_recover_event_rates(Gamma)
    sim.add_exposed_infected_event_rates(Theta_EI)
    sim.run_sim()
    ts, infect_ts, recover_ts, exposed_ts = sim.tabulate_continuous_time(1000)

    plt.figure(1)
    plt.plot(ts, infect_ts, color='red', label='Infected')
    plt.plot(ts, exposed_ts, color='orange', label='Exposed')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.show()

    return A

def create_stochastic_block_model():
    A = np.zeros((300, 300))
    staff_block = 0.8
    residents_block = 0.9
    county_block = 0.1
    staff_res_block = 0.7
    res_county_block = 0.001
    staff_county_block = 0.1 #TODO how to determine this rate

    staff_pop = int(len(A)/3)
    res_pop = int(len(A)/3)
    county_pop = int(len(A)/3)

    staff_range = np.arange(0, staff_pop)
    res_range = np.arange(staff_pop, staff_pop+res_pop)
    county_range = np.arange(staff_pop+res_pop, staff_pop+res_pop+county_pop)

    for i in staff_range:
        for j in staff_range[i:]:
            A[i,j] = 1 if np.random.rand() < staff_block else 0
            A[j,i] = A[i,j]
        for j in res_range:
            A[i,j] = 1 if np.random.rand() < staff_res_block else 0
            A[j,i] = A[i,j]
        for j in county_range:
            A[i,j] = 1 if np.random.rand() < staff_county_block else 0
            A[j,i] = A[i,j]
    for i in res_range:
        for j in res_range[i:]:
            A[i,j] = 1 if np.random.rand() < residents_block else 0
            A[j,i] = A[i,j]
        for j in county_range:
            A[i,j] = 1 if np.random.rand() < res_county_block else 0
            A[j,i] = A[i,j]
    for i in county_range:
        for j in county_range[i:]:
            A[i,j] = 1 if np.random.rand() < county_block else 0
            A[j,i] = A[i,j]

    np.fill_diagonal(A, 0)
    G, pos = network.NetworkBuilder.from_adjacency_matrix(A, return_pos=True)
    G.remove_nodes_from(list(nx.isolates(G)))

    node_colors = {}
    for n in G.nodes():
        if n in staff_range:
            node_colors[n] = 'blue'
        elif n in res_range:
            node_colors[n] = 'red'
        elif n in county_range:
            node_colors[n] = 'orange'
    nx.draw_networkx_nodes(G, pos, G.nodes(), node_color=node_colors.values(), node_size=20)
    nx.draw_networkx_edges(G, pos)
    plt.show()
    return G, pos, A

def power_law_degree_distrb(maxk=40, alpha=2, mu=5):
    p_k = np.empty(maxk)
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-alpha)) * (math.e ** (-k / mu))
    p_k = p_k / np.sum(p_k)
    return p_k

if __name__=='__main__':
    sbm_membership()
    run()
    random_vaccination()




