import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
import networkx as nx
from epintervene.simobjects import simulation
import math


def basic_SIR_simulation():
    print('Basic SIR simulation tutorial')
    # # Create a NetworkBuilder object, to help with creating and configuring the network for the simulation:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 3)
    # # Generating a network from the above degree distribution, with 100 nodes
    G, pos = nb.from_degree_distribution(100, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist = nb.create_adjacency_list(G)

    # # Constructing the simulation object
    my_simulation = simulation.Simulation(adj_list=adjlist, N=len(adjlist))

    # # Setting required configurations
    my_simulation.set_uniform_beta(beta=0.9) # .9 people per day per infected person
    my_simulation.set_uniform_gamma(gamma=0.2) # 5 days to recover

    # # Running the simulation
    my_simulation.run_sim()

    # # Obtaining time series results
    time_series_vals, infected_time_series, recovered_time_series = my_simulation.tabulate_continuous_time()
    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected')
    plt.plot(time_series_vals, recovered_time_series, label='recovered')
    plt.legend(loc='upper left')
    plt.title('Auto generated time series values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining time series results with a custom time series end point.
    # # ****TIP: Run at least one simulation without specifying a custom limit, in order to get a sense of what range
    # # of continuous time
    # # the simulation takes to fully run. Then on a second run after experimentation, you can specify a limit.
    time_series_vals, infected_time_series, recovered_time_series = my_simulation.tabulate_continuous_time(
        time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected')
    plt.plot(time_series_vals, recovered_time_series, label='recovered')
    plt.legend(loc='upper left')
    plt.title('Custom increments of time values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining cumulative infections by generations of infection:
    gen_results = my_simulation.tabulate_generation_results(max_gens=30)

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results, label='cumulative infections')
    plt.plot(np.arange(30), gen_results)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.show()



def basic_SEIR_simulation():
    print('Basic SEIR simulation')
    # # Create a NetworkBuilder object, to help with creating and configuring the network for the simulation:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 3)
    # # Generating a network from the above degree distribution
    G, pos = nb.from_degree_distribution(100, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist = nb.create_adjacency_list(G)

    # # Constructing the simulation object
    my_simulation = simulation.SimulationSEIR(N=len(adjlist), adjlist=adjlist)

    # # Setting required configurations
    my_simulation.set_uniform_beta(beta=0.9)
    my_simulation.set_uniform_beta_es(beta=0.45) # Exposed nodes half as infectious as fully infected
    my_simulation.set_uniform_gamma(gamma=0.1) # 10 days to recover
    my_simulation.set_uniform_gamma_ei(gamma=0.2) # 5 days to develop symptoms

    # # Running the simulation
    my_simulation.run_sim(wait_for_recovery=True) # Setting wait for recovery to True, to observe the full epidemic curve

    # # Obtaining time series results
    time_series_vals, infected_time_series, recovered_time_series, exposed_time_series = my_simulation.tabulate_continuous_time()

    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected')
    plt.plot(time_series_vals, recovered_time_series, label='recovered')
    plt.plot(time_series_vals, exposed_time_series, label='exposed')
    plt.legend(loc='upper left')
    plt.title('Auto generated time series values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining time series results with a custom time series end point.
    # # ****TIP: Run at least one simulation without specifying a custom limit, in order to get a sense of what range
    # # of continuous time
    # # the simulation takes to fully run. Then on a second run after experimentation, you can specify a limit.
    time_series_vals, infected_time_series, recovered_time_series, exposed_time_series = my_simulation.tabulate_continuous_time(
        time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected')
    plt.plot(time_series_vals, recovered_time_series, label='recovered')
    plt.plot(time_series_vals, exposed_time_series, label='exposed')
    plt.legend(loc='upper left')
    plt.title('Custom increments of time values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining cumulative infections by generations of infection:
    gen_results = my_simulation.tabulate_generation_results(max_gens=30)

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results, label='cumulative infections')
    plt.plot(np.arange(30), gen_results)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.show()


def basic_SIR_groups_simulation():
    print('Basic SIR simulation with membership groups')
    # Creating a network from a Stochastic Block Model
    nb = network.NetworkBuilder
    G, pos, _ = create_zoo_stochastic_block_model(tiger_elephant_block=.05, tiger_bird_block=.03,
                                                  bird_elephant_block=.02)
    adjlist = nb.create_adjacency_list(G)

    # # Assigning indices for the split population
    tiger_population = int(len(adjlist) / 3)
    bird_population = int(len(adjlist) / 3)
    elephant_population = int(len(adjlist) / 3)

    # # Creating a vector of node memberships, indexed by node label (index from above)
    node_membership_vector = []
    for i in range(tiger_population):
        node_membership_vector.append('tiger')
    for j in range(bird_population):
        node_membership_vector.append('bird')
    for k in range(elephant_population):
        node_membership_vector.append('elephant')

    # # Setting up the simulation with memberships
    sim = simulation.Simulation(N=len(adjlist), adj_list=adjlist, membership_groups=['tiger', 'bird', 'elephant'],
                                node_memberships=node_membership_vector)

    # # Configuring the disease parameters
    sim.set_uniform_beta(0.7)
    sim.set_uniform_gamma(0.21)

    # # Running the simulation, must specify to track memberships
    sim.run_sim(with_memberships=True, wait_for_recovery=True, uniform_rate=True)

    # # Recording results with group memberships
    ts, membership_ts_infc, membership_ts_rec = sim.tabulate_continuous_time_with_groups(time_buckets=1000, custom_range=True,
                                                                      custom_t_lim=30)

    # # Plotting and interpreting the results
    plt.figure(0)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=f'{group} infected')
        plt.plot(ts, membership_ts_rec[group], label=f'{group} recovered', ls='--')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous time results with group membership')
    plt.box(on=False)
    plt.show()

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=30)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.title('SIR Continuous time results \n(total population, without showing group membership)')
    plt.legend(loc='upper left')
    plt.box(on=False)
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Cumulative Epidemic Generation size results')
    plt.box(on=False)
    plt.show()


def basic_SEIR_groups_simulation():
    print('Basic SEIR simulation with membership groups')
    # Creating a network from a Stochastic Block Model
    nb = network.NetworkBuilder
    G, pos, _ = create_zoo_stochastic_block_model(tiger_elephant_block=.05, tiger_bird_block=.03,
                                                  bird_elephant_block=.02)
    adjlist = nb.create_adjacency_list(G)

    # # Assigning indices for the split population
    tiger_population = int(len(adjlist) / 3)
    bird_population = int(len(adjlist) / 3)
    elephant_population = int(len(adjlist) / 3)

    # # Creating a vector of node memberships, indexed by node label (index from above)
    node_membership_vector = []
    for i in range(tiger_population):
        node_membership_vector.append('tiger')
    for j in range(bird_population):
        node_membership_vector.append('bird')
    for k in range(elephant_population):
        node_membership_vector.append('elephant')

    # # Setting up the simulation with memberships
    sim = simulation.SimulationSEIR(N=len(adjlist), adjlist=adjlist, node_memberships=node_membership_vector,
                                    membership_groups=['tiger', 'bird', 'elephant'])

    # # Configuring the disease parameters
    sim.set_uniform_gamma(0.1)
    sim.set_uniform_beta(0.5)
    sim.set_uniform_beta_es(0.25)
    sim.set_uniform_gamma_ei(0.2)

    # # Running the simulation
    sim.run_sim(with_memberships=True, wait_for_recovery=True, uniform_rate=True)

    # # Recording the results
    ts, membership_ts_infc, membership_ts_exp, membership_ts_rec = sim.tabulate_continuous_time_with_groups()

    # # Plotting and interpreting the results
    plt.figure(0)
    colors = {'tiger': 'blue', 'bird': 'red', 'elephant': 'orange'}
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=f'{group} infected', color=colors[group])
        plt.plot(ts, membership_ts_exp[group], label=f'{group} exposed', ls='--', color=colors[group])
        plt.plot(ts, membership_ts_rec[group], label=f'{group} recovered', ls=':', color=colors[group])
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper right')
    plt.title('SEIR Continuous time results with group membership')
    plt.show()

    ts, infect_ts, recover_ts, exposed_ts = sim.tabulate_continuous_time(1000, custom_range=True, custom_t_lim=30)

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


def create_zoo_stochastic_block_model(N=300, tiger_bird_block=0.01, bird_elephant_block=0.02,
                                      tiger_elephant_block=0.005):
    A = np.zeros((N, N))
    tiger_block = 10 / N
    bird_block = 10 / N
    elephant_block = 20 / N
    tiger_bird_block = tiger_bird_block / N
    bird_elephant_block = bird_elephant_block / N
    tiger_elephant_block = tiger_elephant_block / N

    tiger_pop = int(len(A) / 3)
    bird_pop = int(len(A) / 3)
    elephant_pop = int(len(A) / 3)

    tiger_range = np.arange(0, tiger_pop)
    bird_range = np.arange(tiger_pop, tiger_pop + bird_pop)
    elephant_range = np.arange(tiger_pop + bird_pop, tiger_pop + bird_pop + elephant_pop)

    for i in tiger_range:
        for j in range(i, tiger_range[-1]):
            A[i, j] = 1 if np.random.rand() < tiger_block else 0
            A[j, i] = A[i, j]
        for j in bird_range:
            A[i, j] = 1 if np.random.rand() < tiger_bird_block else 0
            A[j, i] = A[i, j]
        for j in elephant_range:
            A[i, j] = 1 if np.random.rand() < tiger_elephant_block else 0
            A[j, i] = A[i, j]
    for i in bird_range:
        for j in range(i, bird_range[-1]):
            A[i, j] = 1 if np.random.rand() < bird_block else 0
            A[j, i] = A[i, j]
        for j in elephant_range:
            A[i, j] = 1 if np.random.rand() < bird_elephant_block else 0
            A[j, i] = A[i, j]
    for i in elephant_range:
        for j in range(i, elephant_range[-1]):
            A[i, j] = 1 if np.random.rand() < elephant_block else 0
            A[j, i] = A[i, j]

    np.fill_diagonal(A, 0)
    G, pos = network.NetworkBuilder.from_adjacency_matrix(A, return_pos=True)

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
    plt.box(on=False)
    # plt.show()
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


if __name__ == '__main__':
    basic_SIR_simulation()
    basic_SEIR_simulation()
    basic_SIR_groups_simulation()
    basic_SEIR_groups_simulation()
