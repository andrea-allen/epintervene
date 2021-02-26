
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

    A = np.random.random_integers(0, 1, (10, 10))
    A = (A + A.T) / 2
    sim = simulation.Simulation(A)
    sim.add_infection_event_rates(np.zeros((2, 2)))
    return A


