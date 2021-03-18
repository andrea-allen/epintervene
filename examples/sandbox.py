import time
import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
import networkx as nx
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
import math
import time

def optimizing():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb(mu=100)
    start_time=time.time()
    degree_distrb = binomial_degree_distb(400, 2.5)
    degree_distrb = powerlaw
    print(f'net work time {time.time()-start_time}')


    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, degree_distrb)
    adjlist = nb.create_adjacency_list(G)

    A = np.array(nx.adjacency_matrix(G).todense())
    sim = simulation.Simulation(adj_matrix=A, adj_list=adjlist, N=len(A))
    # sim.set_adjlist(adjlist)
    Beta = np.full((len(A), len(A)), 0.1)
    Gamma = np.full(len(A), 0.001)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    # 46 seconds for a major sim with 10000 nodes, 1.5 mean degree, and beta of 0.9
    # commenting out the update_IS_events method, sim takes same time, seemed to have no effect
    start_time = time.time()
    # sim.run_sim(wait_for_recovery=True)
    sim.run_sim(wait_for_recovery=False, uniform_rate=True)
    print(f'Total time for a single sim took {time.time()-start_time}')


    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous Time Results for Random Intervention Simulation')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Generational Cumulative Results for Random Intervention Simulation')
    plt.show()


def random_vaccination():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb(mu=100)

    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, powerlaw)
    adjlist = nb.create_adjacency_list(G)
    A = np.array(nx.adjacency_matrix(G).todense())
    sim = extended_simulation.RandomRolloutSimulation(N=len(A), adjmatrix=A, adjlist=adjlist)
    # sim = extended_simulation.RandomInterventionSim(N=len(A), adjmatrix=A, adjlist=adjlist)
    Beta = np.full((len(A), len(A)), 0.9)
    Gamma = np.full(len(A), 0.001)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    sim.configure_intervention(intervention_gen_list=[3,4], beta_redux_list=[0,0], proportion_reduced_list=[0.01, 0.03])
    # sim.configure_intervention(3, 0, .9) #TODO ready to commit, experiment w intervention rollouts
    start = time.time()
    sim.run_sim(wait_for_recovery=False, uniform_rate=True)
    print(f'total sim time {time.time()-start}')

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous Time Results for Random Intervention Simulation')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Generational Cumulative Results for Random Intervention Simulation')
    plt.show()

def membership():
    # Creating a network from a Stochastic Block Model
    nb = network.NetworkBuilder
    G, pos, A = create_zoo_stochastic_block_model()
    adjlist = nb.create_adjacency_list(G)
    # SIR model sandbox, want to have different curves for different group memberships

    tiger_population = int(len(A)/3)
    bird_population = int(len(A)/3)
    elephant_population = int(len(A)/3)

    node_membership_vector = []
    for i in range(tiger_population):
        node_membership_vector.append('tiger')
    for j in range(bird_population):
        node_membership_vector.append('bird')
    for k in range(elephant_population):
        node_membership_vector.append('elephant')
    sim = simulation.Simulation(N=len(A), adj_matrix=A, adj_list=adjlist, membership_groups=['tiger', 'bird', 'elephant'], node_memberships=node_membership_vector)
    Beta = np.full((len(A), len(A)), 0.5)
    Gamma = np.full(len(A), 0.9)
    sim.add_infection_event_rates(Beta)
    sim.add_recover_event_rates(Gamma)
    sim.run_sim(with_memberships=True, wait_for_recovery=True)

    ts, membership_ts_infc = sim.tabulate_continuous_time_with_groups(time_buckets=1000, custom_range=True, custom_t_lim=15)
    plt.figure(0)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group)
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous time results with group membership')
    plt.show()

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=10)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.title('SIR Continuous time results (total population, without showing group membership)')
    plt.legend(loc='upper left')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Cumulative Epidemic Generation size results')
    plt.show()


    # SEIR model sandbox
    sim = simulation.SimulationSEIR(N=len(A), adjmatrix=A, adjlist=adjlist, node_memberships=node_membership_vector, membership_groups=['tiger', 'bird', 'elephant'])
    Beta_IS = np.full((len(A), len(A)), 0.25)
    Gamma = np.full(len(A), 1.0)
    Beta_ES = np.full((len(A), len(A)), 0.25)
    Theta_EI = np.full(len(A), 1.0)
    sim.add_infection_event_rates(Beta_IS)
    sim.add_exposed_event_rates(Beta_ES)
    sim.add_recover_event_rates(Gamma)
    sim.add_exposed_infected_event_rates(Theta_EI)

    adjlist = nb.create_adjacency_list(G)
    sim.set_adjlist(adjlist)

    sim.run_sim(with_memberships=True, wait_for_recovery=True)

    ts, membership_ts_infc, membership_ts_exp = sim.tabulate_continuous_time_with_groups(1000, custom_range=True, custom_t_lim=10)
    plt.figure(0)
    colors = {'tiger':'blue', 'bird':'red', 'elephant':'orange'}
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group, color=colors[group])
        plt.plot(ts, membership_ts_exp[group], label=group, ls='--', color=colors[group])
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.title('SEIR Continuous time results with group membership')
    plt.show()



    ts, infect_ts, recover_ts, exposed_ts = sim.tabulate_continuous_time(1000, custom_range=True, custom_t_lim=10)

    plt.figure(1)
    plt.plot(ts, infect_ts, color='red', label='Infected')
    plt.plot(ts, exposed_ts, color='orange', label='Exposed')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    plt.title('SEIR Continuous time results (total population, without showing group membership)')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(max_gens=20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SEIR Cumulative Epidemic Generation size results')
    plt.show()

def create_zoo_stochastic_block_model():
    A = np.zeros((300, 300))
    tiger_block = 0.15
    bird_block = 0.15
    elephant_block = 0.2
    tiger_bird_block = 0.01
    bird_elephant_block = 0.02
    tiger_elephant_block = 0.005

    tiger_pop = int(len(A)/3)
    bird_pop = int(len(A)/3)
    elephant_pop = int(len(A)/3)

    tiger_range = np.arange(0, tiger_pop)
    bird_range = np.arange(tiger_pop, tiger_pop+bird_pop)
    elephant_range = np.arange(tiger_pop+bird_pop, tiger_pop+bird_pop+elephant_pop)

    for i in tiger_range:
        for j in range(i, tiger_range[-1]):
            A[i,j] = 1 if np.random.rand() < tiger_block else 0
            A[j,i] = A[i,j]
        for j in bird_range:
            A[i,j] = 1 if np.random.rand() < tiger_bird_block else 0
            A[j,i] = A[i,j]
        for j in elephant_range:
            A[i,j] = 1 if np.random.rand() < tiger_elephant_block else 0
            A[j,i] = A[i,j]
    for i in bird_range:
        for j in range(i, bird_range[-1]):
            A[i,j] = 1 if np.random.rand() < bird_block else 0
            A[j,i] = A[i,j]
        for j in elephant_range:
            A[i,j] = 1 if np.random.rand() < bird_elephant_block else 0
            A[j,i] = A[i,j]
    for i in elephant_range:
        for j in range(i, elephant_range[-1]):
            A[i,j] = 1 if np.random.rand() < elephant_block else 0
            A[j,i] = A[i,j]

    np.fill_diagonal(A, 0)
    G, pos = network.NetworkBuilder.from_adjacency_matrix(A, return_pos=True)
    G.remove_nodes_from(list(nx.isolates(G)))

    node_colors = {}
    for n in G.nodes():
        if n in tiger_range:
            node_colors[n] = 'blue'
        elif n in bird_range:
            node_colors[n] = 'red'
        elif n in elephant_range:
            node_colors[n] = 'orange'
    nx.draw_networkx_nodes(G, pos, G.nodes(), node_color=node_colors.values(), node_size=20)
    nx.draw_networkx_edges(G, pos)
    plt.title('Tigers (blue), Bird (red), Elephant (orange) network')
    plt.show()
    return G, pos, A

def power_law_degree_distrb(maxk=40, alpha=2, mu=5):
    p_k = np.empty(maxk)
    p_k[0] = 0
    for k in range(1, maxk):
        p_k[k] = (k ** (-alpha)) * (math.e ** (-k / mu))
    p_k = p_k / np.sum(p_k)
    return p_k

def binomial_degree_distb(N, lam=6):
    p_k = np.empty(N)
    p = lam / N
    for k in range(0, len(p_k)):
        p_k[k] = (p ** k) * ((1 - p) ** (N - k)) * math.comb(N, k)
    return p_k

if __name__=='__main__':
    optimizing()
    membership()
    random_vaccination()




