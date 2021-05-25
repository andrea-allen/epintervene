import time
import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
import networkx as nx
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
import math
import time
import random

def basic_simulation():
    # TODO
    print('Basic simulation')

def visualize_network():
    nb = network.NetworkBuilder
    # powerlaw = power_law_degree_distrb(mu=100)
    degree_distrb = binomial_degree_distb(400, 6)

    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(45, degree_distrb, True)
    adjlist = nb.create_adjacency_list(G)
    N = len(G.nodes())
    sim = simulation.Simulation(adj_list=adjlist, N=N)
    sim.set_uniform_beta(0.9)
    sim.set_uniform_gamma(0.001)
    sim.run_sim(wait_for_recovery=False, uniform_rate=True, visualize=True, viz_graph=G, viz_pos=pos)

def chain_degree_dist(N):
    p_k = np.zeros(N)
    p_k[0] = 0
    p_k[1] = .01
    p_k[2] = .99
    return p_k
def chain_network():
    nb = network.NetworkBuilder
    chain_net = chain_degree_dist(10)
    plaw = power_law_degree_distrb(400)
    start_time=time.time()
    total_counting_time = 0
    total_sim_time = 0
    degree_distrb = binomial_degree_distb(400, 2.5)
    # degree_distrb = plaw
    print(f'net work time {time.time()-start_time}')


    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, degree_distrb)
    # G, pos = nb.from_degree_distribution(10000, plaw)
    # G = nx.generators.balanced_tree(2, 5)
    # G = nx.Graph()
    # G.add_nodes_from(list(np.arange(40)))
    # for i in range(39):
    #     G.add_edge(i, i+1)
    # G = nx.generators.balanced_tree(2, 5)
    # G.add_edges_from([(0,63), (63,64), (63,65), (64,66), (64,67), (65,68), (68,69)])
    # G = nx.generators.complete_graph(10)
    # nx.draw(G, with_labels=True)
    # plt.show()

    ### STAR NETWORK
    # G = nx.Graph()
    # G.add_nodes_from([0])
    #
    # arm_length = 5
    # range1 = np.arange(arm_length, arm_length+arm_length)
    # range2 = np.arange(arm_length+arm_length, arm_length+arm_length*2)
    # range3 = np.arange(arm_length+arm_length*2, arm_length+arm_length*3)
    # range4 = np.arange(arm_length+arm_length*3, arm_length+arm_length*4)
    # G.add_edge(0,1)
    # G.add_edge(0,2)
    # G.add_edge(0,3)
    # G.add_edge(0,4)
    # G.add_edge(1, range1[0])
    # G.add_edge(2, range2[0])
    # G.add_edge(3, range3[0])
    # G.add_edge(4, range4[0])
    # for j in range(len(range1)-1):
    #     G.add_edge(range1[j], range1[j+1])
    #     G.add_edge(range2[j], range2[j+1])
    #     G.add_edge(range3[j], range3[j+1])
    #     G.add_edge(range4[j], range4[j+1])
    # nx.draw(G, with_labels=True)
    # plt.show()
    #
    #
    # pos = nx.spring_layout(G)
    adjlist = nb.create_adjacency_list(G)

    A = np.array(nx.adjacency_matrix(G).todense())

    # dd_g = np.array(nx.degree_histogram(G)) / (
    #     sum(nx.degree_histogram(G)))
    # avg = 0
    # for k in range(len(dd_g)):
    #     avg += k * dd_g[k]
    #
    # Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    # G0 = G.subgraph(Gcc[0])
    # dd_gcc = np.array(nx.degree_histogram(G0)) / (
    #     sum(nx.degree_histogram(G0)))
    # Gcc_avg = 0
    # for k in range(len(dd_gcc)):
    #     Gcc_avg += k * dd_gcc[k]
    # dd_g1 = g1_of(dd_gcc)
    # Gcc_avg_excess = z1_of(dd_g1)


    start_time = time.time()
    total_infectious_degree_counter = 0
    number_took_off = 0
    averaging_active = None
    average_active_sizes = None
    averaging_total = None
    averaging_count = 0
    gen_results_average = np.zeros(100)
    full_ts = None
    full_ts_infc = None
    full_ts_rec = None
    generational_distribution = np.zeros((100, len(A)))

    num_sims = 30

    for i in range(num_sims):
        if i % 30 == 0:
            print(f'sim number{i}')
        sim = simulation.Simulation(adj_list=adjlist, N=len(A))
        sim.set_uniform_beta(beta=0.004)
        sim.set_uniform_gamma(gamma=0.001)

        # time_s = time.time()
        sim.run_sim(wait_for_recovery=False, uniform_rate=True, viz_pos=pos, viz_graph=G, visualize=False, p_zero=None,
                    kill_by=16, record_active_gen_sizes=False)
        # actual_run_time = time.time()-time_s
        # total_sim_time += time.time() - time_s
        # print(f'sim time was {actual_run_time}')
        # current_counter = sim.infectious_degree_counter
        # if len(current_counter)>0:
        #     print(f'len of counter {len(current_counter)}')
        #     if len(current_counter)>100:
        #         total_infectious_degree_counter += sum(sim.infectious_degree_counter)/(len(sim.infectious_degree_counter))
        #         number_took_off += 1
        # time_s = time.time()
        ts, infect_ts, recover_ts, active_gen_ts, total_gen_ts = sim.tabulate_continuous_time(1000, custom_range=True, custom_t_lim=5000, active_gen_info=True)
        # ts, infect_ts, recover_ts, active_gen_ts, total_gen_ts, active_sizes_ts = sim.tabulate_continuous_time(1000, custom_range=True, custom_t_lim=5000, active_gen_info=True, active_gen_sizes=True)
        # ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000, custom_range=True, custom_t_lim=5000)
        if full_ts is None:
            full_ts = ts
            full_ts_infc = infect_ts
            full_ts_rec = recover_ts
        else:
            full_ts_infc += infect_ts
            full_ts_rec += recover_ts
        # print(f'result counting time was {time.time()-time_s}')
        # total_counting_time += time.time()-time_s
        # plt.plot(ts, infect_ts, color='red')
        # plt.plot(ts, recover_ts, color='green')
        # plt.show()
        # time_s = time.time()
        gen_results = sim.tabulate_generation_results(100)
        # print(f'generation counting time was {time.time()-time_s}')
        if averaging_active is None:
            averaging_active = np.array(active_gen_ts)
            averaging_total = np.array(total_gen_ts)
            # average_active_sizes = active_sizes_ts
        else:
            # print(np.array(total_gen_ts)[-1])
            if total_gen_ts[-1] > 2:
                averaging_active += np.array(active_gen_ts)
                averaging_total += np.array(total_gen_ts)
                # average_active_sizes += active_sizes_ts
                averaging_count += 1
        gen_results_average += gen_results
        for gen in range(len(gen_results)):
            num_total_infctd = int(gen_results[gen])
            try:
                generational_distribution[gen][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', gen, ', gen_s: ', num_total_infctd)
                continue
        # plt.figure(1, frameon=True)
        # plt.plot(ts, infect_ts, color='blue', lw=2, label='Infected')
        # # plt.plot(ts, recover_ts, color='green', label='Recovered')
        # plt.xlabel('Time', fontsize=12)
        # plt.ylabel('Infections', fontsize=12)
        # plt.ylabel('Number of nodes in class')
        # plt.legend(loc='upper left')
        # plt.xticks([])
        # plt.yticks([])
        # plt.title('SIR Continuous Time Results for Random Intervention Simulation')
        # plt.show()

        # ts_by_gen = sim.tabulate_generation_results(20)
        # plt.figure(2)
        # plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
        # plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
        # plt.xlabel('Generation number')
        # plt.ylabel('Cumulative infections by generation')
        # plt.title('SIR Generational Cumulative Results for Random Intervention Simulation')
        # plt.show()
    # print(active_gen_ts, total_gen_ts)
    # average_active_sizes = average_active_sizes / averaging_count
    # for g in range(1, 12, 3):
    #     plt.plot(average_active_sizes.T[g], label=f'generation {g}')
    plt.legend(loc='upper right')
    plt.xlabel('time')
    plt.ylabel('average number active nodes')
    plt.show()

    plt.plot(ts, averaging_active / averaging_count, label='active gens')
    plt.plot(ts, averaging_total / averaging_count, label='total gens')
    plt.legend(loc='upper left')
    plt.show()
    plt.title('gen results')
    plt.scatter(np.arange(100), gen_results_average/num_sims)
    plt.show()

    plt.plot(full_ts, full_ts_infc/num_sims)
    plt.plot(full_ts, full_ts_rec/num_sims)
    plt.show()
    # averaging results:
    for gen in range(100):
        gen_time_series = generational_distribution[gen]
        gen_time_series = gen_time_series / num_sims
        generational_distribution[gen] = gen_time_series
    plt.title('Generation distributions')
    plt.title('Generation distributions')
    plt.plot(np.arange(390), generational_distribution[3][2:392])
    plt.plot(np.arange(390), generational_distribution[6][2:392])
    plt.plot(np.arange(390), generational_distribution[10][2:392])
    plt.plot(np.arange(390), generational_distribution[15][2:392])
    plt.semilogy()
    plt.ylim(.00005, .1)
    plt.show()
    # print(f'total infectious degree avg {total_infectious_degree_counter/number_took_off}')
    print(f'Total time for all sim took {time.time()-start_time}')
    print(f'total sim time {total_sim_time}')
    print(f'total counting time {total_counting_time}')


def sim_testing():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb()
    start_time=time.time()
    degree_distrb = binomial_degree_distb(400, 3)
    # degree_distrb = powerlaw
    print(f'net work time {time.time()-start_time}')


    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(30, degree_distrb)
    adjlist = nb.create_adjacency_list(G)

    for i in range(4):
        print(i)
        sim = simulation.Simulation(adj_list=adjlist, N=len(G.nodes()))
        # sim = extended_simulation.RandomRolloutSimulation(adjmatrix=A, adjlist=adjlist, N=len(A))
        # sim.set_adjlist(adjlist)
        # Beta = np.full((len(A), len(A)), 0.0015)
        # Gamma = np.full(len(A), 0.001)
        # sim.add_infection_event_rates(Beta)
        # sim.add_recover_event_rates(Gamma)
        sim.set_uniform_beta(0.004)
        sim.set_uniform_gamma(0.001)
        sim.run_sim(wait_for_recovery=False, uniform_rate=True, kill_by=None)
        ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
        plt.figure(1, frameon=True)
        plt.plot(ts, infect_ts, color='blue', lw=3, label='Infected')
        # plt.plot(ts, recover_ts, color='green', label='Recovered')
        plt.xlabel('Time', fontsize=18)
        plt.ylabel('Infections', fontsize=18)
        # plt.ylabel('Number of nodes in class')
        # plt.legend(loc='upper left')
        plt.xticks([])
        plt.yticks([])
        # plt.title('SIR Continuous Time Results for Random Intervention Simulation')
        # plt.show()

        ts_by_gen = sim.tabulate_generation_results(20)
        plt.figure(2)
        plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
        plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
        plt.xlabel('Generation number')
        plt.ylabel('Cumulative infections by generation')
        plt.title('SIR Generational Cumulative Results for Random Intervention Simulation')
        plt.show()


def optimizing():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb()
    start_time=time.time()
    degree_distrb = binomial_degree_distb(400, 2.5)
    degree_distrb = powerlaw
    print(f'net work time {time.time()-start_time}')


    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, degree_distrb)
    adjlist = nb.create_adjacency_list(G)

    A = np.array(nx.adjacency_matrix(G).todense())
    start_ensemble_time = time.time()
    total_total_time = 0
    for i in range(1000):
        # print(i)
        sim = simulation.Simulation(adj_list=adjlist, N=len(A))
        # sim = extended_simulation.RandomRolloutSimulation(adjmatrix=A, adjlist=adjlist, N=len(A))
        # sim.set_adjlist(adjlist)
        # Beta = np.full((len(A), len(A)), 0.0015)
        # Gamma = np.full(len(A), 0.001)
        # sim.add_infection_event_rates(Beta)
        # sim.add_recover_event_rates(Gamma)
        sim.set_uniform_beta(0.004)
        sim.set_uniform_gamma(0.001)
        # sim.configure_intervention(intervention_gen_list=[5,6,7], beta_redux_list=[0, 0, 0], proportion_reduced_list=[0.01,0.05,0.10])
        # 46 seconds for a major sim with 10000 nodes, 1.5 mean degree, and beta of 0.9
        # commenting out the update_IS_events method, sim takes same time, seemed to have no effect
        # sim.run_sim(wait_for_recovery=True)
        start_time = time.time()
        sim.run_sim(wait_for_recovery=False, uniform_rate=True, kill_by=16)
        ## TODO how much is because of that random exponential vs the uniforms
        total_total_time += time.time()-start_time
        # print(sum([len(sim.get_gen_collection()[key]) for key in list(sim.get_gen_collection().keys())]))
        # if max(sim._real_time_srs_infc) > 100:
        #     plt.plot(sim._time_series[1:], sim._real_time_srs_infc)
        #     plt.plot(sim._time_series[1:], sim._real_time_srs_rec)
        #     plt.plot(sim._time_series[1:], np.array(sim._real_time_srs_infc)+np.array(sim._real_time_srs_rec))
        #     plt.show()
        # TODO: turns out 2.5 seconds (35 percent) is taken by not single step. What's it coming from?
        # TODO 44 percent of time from getting the length of the events in the dictionary
    print(f'Total time for all sim start to end took {time.time()-start_ensemble_time}')


    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
    plt.figure(1, frameon=True)
    plt.plot(ts, infect_ts, color='blue', lw=3, label='Infected')
    # plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time', fontsize=18)
    plt.ylabel('Infections', fontsize=18)
    # plt.ylabel('Number of nodes in class')
    # plt.legend(loc='upper left')
    plt.xticks([])
    plt.yticks([])
    # plt.title('SIR Continuous Time Results for Random Intervention Simulation')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Generational Cumulative Results for Random Intervention Simulation')
    plt.show()

def uniform_reduction():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb(mu=100)

    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, powerlaw)
    adjlist = nb.create_adjacency_list(G)
    A = np.array(nx.adjacency_matrix(G).todense())
    sim = extended_simulation.UniversalInterventionSim(N=len(A), adjlist=adjlist)
    sim.configure_intervention(intervention_gen=4, beta_redux=0.6)
    sim.set_uniform_beta(0.9)
    sim.set_uniform_gamma(0.1)

    start = time.time()
    sim.run_sim(wait_for_recovery=False, uniform_rate=True)
    print(f'total sim time {time.time()-start}')

    reg_sim = simulation.Simulation(N=len(A), adj_list=adjlist)
    reg_sim.set_uniform_beta(0.9)
    reg_sim.set_uniform_gamma(0.1)

    start = time.time()
    reg_sim.run_sim(wait_for_recovery=False, uniform_rate=True)
    print(f'total sim time {time.time()-start}')


    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(1000)
    reg_ts, reg_infect_ts, reg_recover_ts = reg_sim.tabulate_continuous_time(1000)
    plt.figure(1)
    plt.plot(ts, infect_ts, color='blue', label='Intv. Infected', ls='--')
    plt.plot(ts, reg_infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Intv. Recovered', ls='--')
    plt.plot(ts, reg_recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous Time Results for Universal Intervention Simulation')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    ts_by_gen_reg = reg_sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.plot(np.arange(len(ts_by_gen_reg)), ts_by_gen_reg)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen_reg)), ts_by_gen_reg)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Generational Cumulative Results for Universal Intervention Simulation')
    plt.show()


def random_rollout_vaccination():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb(mu=100)

    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, powerlaw)
    adjlist = nb.create_adjacency_list(G)
    A = np.array(nx.adjacency_matrix(G).todense())
    sim = extended_simulation.RandomRolloutSimulation(N=len(A), adjmatrix=A, adjlist=adjlist)
    # sim = extended_simulation.RandomInterventionSim(N=len(A), adjmatrix=A, adjlist=adjlist)
    # Beta = np.full((len(A), len(A)), 0.9)
    # Gamma = np.full(len(A), 0.001)
    # sim.add_infection_event_rates(Beta)
    # sim.add_recover_event_rates(Gamma)
    sim.set_uniform_beta(0.9)
    sim.set_uniform_gamma(0.001)
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
    plt.title('SIR Continuous Time Results for Random Rollout Intervention Simulation')
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Generational Cumulative Results for Random Rollout Intervention Simulation')
    plt.show()

def random_vaccination():
    nb = network.NetworkBuilder
    powerlaw = power_law_degree_distrb(mu=100)

    # Creating a network from a power law degree distribution
    G, pos = nb.from_degree_distribution(10000, powerlaw)
    adjlist = nb.create_adjacency_list(G)
    A = np.array(nx.adjacency_matrix(G).todense())
    sim = extended_simulation.RandomInterventionSim(N=len(A), adjmatrix=A, adjlist=adjlist)
    # sim = extended_simulation.RandomInterventionSim(N=len(A), adjmatrix=A, adjlist=adjlist)
    # Beta = np.full((len(A), len(A)), 0.9)
    # Gamma = np.full(len(A), 0.001)
    # sim.add_infection_event_rates(Beta)
    # sim.add_recover_event_rates(Gamma)
    sim.set_uniform_beta(0.9)
    sim.set_uniform_gamma(0.001)
    sim.configure_intervention(intervention_gen=4, beta_redux=0.0, proportion_reduced=0.10)
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
    G, pos, A = create_zoo_stochastic_block_model(tiger_elephant_block=.05, tiger_bird_block=.03, bird_elephant_block=.02)
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
    sim = simulation.Simulation(N=len(A), adj_list=adjlist, membership_groups=['tiger', 'bird', 'elephant'],
                                node_memberships=node_membership_vector)
    sim.set_uniform_beta(0.7)
    sim.set_uniform_gamma(0.21)
    sim.run_sim(with_memberships=True, wait_for_recovery=True, uniform_rate=True)

    ts, membership_ts_infc = sim.tabulate_continuous_time_with_groups(time_buckets=1000, custom_range=True,
                                                                      custom_t_lim=30)
    plt.figure(0)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group)
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous time results with group membership')
    plt.box(on=False)
    plt.show()

    ts, infect_ts, recover_ts = sim.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=30)
    plt.figure(1, frameon=False)
    plt.plot(ts, infect_ts, color='blue', label='Infected')
    plt.plot(ts, recover_ts, color='green', label='Recovered')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes in class')
    plt.title('SIR Continuous time results \n(total population, without showing group membership)')
    plt.legend(loc='upper left')
    plt.box(on=False)
    # plt.show()

    ts_by_gen = sim.tabulate_generation_results(20)
    plt.figure(2, frameon=False)
    plt.plot(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.scatter(np.arange(len(ts_by_gen)), ts_by_gen)
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative infections by generation')
    plt.title('SIR Cumulative Epidemic Generation size results')
    plt.box(on=False)
    plt.show()


    # SEIR model sandbox
    sim = simulation.SimulationSEIR(N=len(A), adjlist=adjlist, node_memberships=node_membership_vector, membership_groups=['tiger', 'bird', 'elephant'])
    # Idea: To set non uniform rates, you could do a regular SIR model with membership groups,
    # and have the membership groups have their own beta/gamma values (i.e. if it spreads more quickly
    # in certain communities, not just because of the network) but then you don't have to track memberships, just use them


    sim.set_uniform_gamma(0.00001)
    sim.set_uniform_beta(0.5)
    sim.set_uniform_beta_es(0.5)
    sim.set_uniform_gamma_ei(1.0)

    # should fix intervention code in the extended simulation package.
    # then need to push and fix all the references. (do after weekend)
    sim.run_sim(with_memberships=True, wait_for_recovery=True, uniform_rate=True)

    ts, membership_ts_infc, membership_ts_exp = sim.tabulate_continuous_time_with_groups(1000, custom_range=True, custom_t_lim=30)
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

def create_zoo_stochastic_block_model(N=300, tiger_bird_block = 0.01, bird_elephant_block = 0.02, tiger_elephant_block = 0.005):
    A = np.zeros((N, N))
    tiger_block = 10/N
    bird_block = 10/N
    elephant_block = 20/N
    tiger_bird_block = tiger_bird_block/N
    bird_elephant_block = bird_elephant_block/N
    tiger_elephant_block = tiger_elephant_block/N

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
    # G.remove_nodes_from(list(nx.isolates(G)))

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

def temporal_intervention():
    # tabulate a bunch of ensemble results
    # Creating a network from a Stochastic Block Model
    intervention_time = 0.5
    nb = network.NetworkBuilder
    G, pos, A = create_zoo_stochastic_block_model(N=3000)
    G_switch, pos, A_switch = create_zoo_stochastic_block_model(N=3000, tiger_bird_block=0,bird_elephant_block=0,tiger_elephant_block=0)
    adjlist = nb.create_adjacency_list(G)
    adjlist_switch = nb.create_adjacency_list(G_switch)
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
    sim = extended_simulation.AbsoluteTimeNetworkSwitchSim(N=len(A), A=A, adjlist=adjlist, membership_groups=['tiger', 'bird', 'elephant'], node_memberships=node_membership_vector)
    sim.configure_intervention(intervention_time=intervention_time, new_adjacency_matrix=A_switch, new_adjlist=adjlist_switch)
    sim.set_uniform_beta(0.5)
    sim.set_uniform_gamma(0.9)
    sim.run_sim(with_memberships=True, wait_for_recovery=False, use_uniform_rate=True)

    sim_regular = simulation.Simulation(N=len(A), adj_matrix=A, adj_list=adjlist, membership_groups=['tiger', 'bird', 'elephant'], node_memberships=node_membership_vector)
    sim_regular.set_uniform_beta(0.5)
    sim_regular.set_uniform_gamma(0.9)
    sim_regular.run_sim(with_memberships=True, wait_for_recovery=False, uniform_rate=True)

    inft_reg_results = np.zeros(1000)
    rec_reg_results = np.zeros(1000)
    inft_results = np.zeros(1000)
    rec_results = np.zeros(1000)
    start_t = time.time()
    for i in range(100):
        print(i)
        sim_regular = simulation.Simulation(N=len(A), adj_matrix=A, adj_list=adjlist,
                                            membership_groups=['tiger', 'bird', 'elephant'],
                                            node_memberships=node_membership_vector)
        sim_regular.set_uniform_beta(0.5)
        sim_regular.set_uniform_gamma(0.9)
        sim_regular.run_sim(with_memberships=True, wait_for_recovery=True, uniform_rate=True)
        ts, infct, rec = sim_regular.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=6)
        inft_reg_results += infct
        rec_reg_results += rec

        sim = extended_simulation.AbsoluteTimeNetworkSwitchSim(N=len(A), A=A, adjlist=adjlist,
                                                               membership_groups=['tiger', 'bird', 'elephant'],
                                                               node_memberships=node_membership_vector)
        sim.configure_intervention(intervention_time=intervention_time, new_adjacency_matrix=A_switch, new_adjlist=adjlist_switch)
        sim.set_uniform_beta(0.5)
        sim.set_uniform_gamma(0.9)
        sim.run_sim(with_memberships=True, wait_for_recovery=True, use_uniform_rate=True)
        ts, infct, rec = sim.tabulate_continuous_time(time_buckets=1000, custom_range=True, custom_t_lim=6)
        inft_results += infct
        rec_results += rec
    print(f'total time was {time.time()-start_t}')

    inft_reg_results = inft_reg_results/100
    rec_reg_results = rec_reg_results/100
    inft_results = inft_results/100
    rec_results = rec_results/100

    plt.figure(1)
    plt.plot(ts, inft_reg_results, label='Regular Infected')
    plt.plot(ts, rec_reg_results, label='Regular Recovered')
    # plt.legend(loc='uppwer right')
    #
    # plt.figure(2)
    plt.plot(ts, inft_results, label='Net Switch Infected')
    plt.plot(ts, rec_results, label='Net Switch Recovered')
    plt.vlines(intervention_time, ymin=0, ymax=np.max(rec_reg_results), ls='--', color='red', label='Network switch time')

    plt.legend(loc='upper right')
    # plt.show()

    plt.figure('difference')
    plt.plot(ts, np.abs(inft_reg_results-inft_results), label='Difference in Infections')
    plt.vlines(intervention_time, ymin=0, ymax=np.max(np.abs(inft_reg_results-inft_results)), ls='--', color='red', label='Network switch time')

    plt.legend(loc='upper left')
    plt.show()

    ts, membership_ts_infc = sim_regular.tabulate_continuous_time_with_groups(time_buckets=1000, custom_range=True, custom_t_lim=6)
    plt.figure(0)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group+' Regular')
    ts, membership_ts_infc = sim.tabulate_continuous_time_with_groups(time_buckets=1000, custom_range=True, custom_t_lim=6)
    for group in membership_ts_infc.keys():
        plt.plot(ts, membership_ts_infc[group], label=group+' Net Switch')
    plt.xlabel('Time t')
    plt.ylabel('Number of nodes infected in network group')
    plt.legend(loc='upper left')
    plt.title('SIR Continuous time results with group membership')
    plt.show()



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

def z1_of(g_0):
    z1 = 0
    for k in range(len(g_0)):
        z1 += (k * g_0[k])
    return z1

def g1_of(g_0):
    g_1 = np.zeros(len(g_0))
    for k in range(len(g_0) - 1):
        g_1[k] = (k + 1) * g_0[k + 1]
    return g_1 / (z1_of(g_0))

def speed_random():
    start_time = time.time()
    for i in range(1000000):
        number = np.random.randint(0, 100)
    numpy_time = time.time() - start_time
    start_time = time.time()
    for i in range(1000000):
        random_num = random.randint(0, 100)
    randint_time = time.time() - start_time
    start_time = time.time()
    for i in range(1000000):
        random_uni = random.uniform(0, 100)
    rand_uni_time = time.time()-start_time
    start_time = time.time()
    for i in range(1000000):
        random_uni = random.uniform(0, 100)
        random_rounded = int(random_uni)
    rand_uni_round_time = time.time()-start_time
    print(f'Numpy random time {numpy_time}')
    print(f'Randint random time {randint_time}')
    print(f'Random Unfirom time {rand_uni_time}')
    print(f'Random Unfirom Rounded time {rand_uni_round_time}')
    start_time = time.time()
    timing_uni = 0
    for i in range(1000000):
        t_s = time.time()
        random_uni = random.uniform(0, 100)
        timing_uni += time.time()-t_s
    time_with_timing = time.time()-start_time
    print(f'Timing Unfirom time {time_with_timing}')
    print(f'Timing time {timing_uni}')

def expovariate_versions():
    distribution = np.zeros(1000000)
    start_time = time.time()
    for i in range(1000000):
        distribution[i] = random.expovariate((10*.004+110*.001))
    total_random_expo = time.time()-start_time

    distribution_np = np.zeros(1000000)
    start_time = time.time()
    for i in range(1000000):
        distribution_np[i] = np.random.exponential(1/(10*.004+110*.001))
    total_np_random_expo = time.time()-start_time

    print(f'Random time {total_random_expo}')
    print(f'Random numpy time {total_np_random_expo}')

    plt.figure('random')
    plt.hist(distribution, bins=100, label='random', alpha=0.5)

    # plt.figure('numpy')
    plt.hist(distribution_np, bins=100, label='numpy random', alpha=0.5)
    plt.legend(loc='upper right')

    plt.show()



if __name__=='__main__':
    # visualize_network()
    # uniform_reduction()
    # chain_network()
    # optimizing()
    # sim_testing()
    membership()
    # speed_random()
    # expovariate_versions()
    # membership()
    # random_vaccination()
    # random_rollout_vaccination()
    # temporal_intervention()




