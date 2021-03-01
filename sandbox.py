
import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
import networkx as nx
from epintervene.simobjects import simulation

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
    # todo make sure all the lists sum to total nodes
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


