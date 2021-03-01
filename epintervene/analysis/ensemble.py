import networkx as nx
import numpy as np
from epintervene.simobjects import simulation
from epintervene.simobjects import extended_simulation
from epintervene.simobjects import network
import time


def simulate_intervention_effects(degree_distrb, base_file_name='sim_results',
                                  num_sims=10000, num_nodes=1000, init_T=0.8,
                                  gen_intervene=3, T_intervene=0.4, recover_rate=.001, prop_reduced=0.0):
    # Ensemble run with no intervention
    start_time = time.time()
    size_distrb_per_gen_no_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                            -1, 0.0, init_T, recover_rate, prop_reduced)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt('./../data/' + base_file_name + '_size_distrb_per_gen_no_interv.txt', size_distrb_per_gen_no_intervention,
               delimiter=',')

    # WITH intervention
    start_time = time.time()
    size_distrb_per_gen_intervention = simulate_ensemble(degree_distrb, num_sims, num_nodes,
                                                         gen_intervene, T_intervene, init_T,
                                                         recover_rate, prop_reduced)
    print('Total time for ensemble was ', time.time() - start_time, ' with N=', num_nodes, ' nodes and ', num_sims,
          ' number of simulations.')
    np.savetxt('./../data/' + base_file_name + '_size_distrb_per_gen_with_interv.txt', size_distrb_per_gen_intervention, delimiter=',')


# Simulate ensemble saves results for generational time series, as opposed to real time results. Those are TBD
def simulate_ensemble(degree_distrb, num_sims=10, N=1000, intervention_gen=-1, intervention_T=0.0, initial_T=0.8,
                      gamma=0.1, prop_reduced=0.0):

    # If intervention_gen is set to -1, no intervention will occur
    # Generates the empirical s-slice for results of proportion of simulations that resulted in exactly s nodes
    # infected at generation g
    # Includes steps for doing so with an intervention step
    beta_init = -(gamma * initial_T) / (initial_T - 1)
    beta_interv = -(gamma * intervention_T) / (intervention_T - 1)
    outbreak_size_distrb_per_gen_matrix = np.zeros((100, N))
    G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
    A = np.array(nx.adjacency_matrix(G).todense())
    num_nodes_in_net = len(A[0])
    Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
    Gamma = np.full(num_nodes_in_net, gamma)
    for i in range(num_sims):
        if i % 500 == 0:
            G, pos = network.NetworkBuilder.from_degree_distribution(N, degree_distrb)
            A = np.array(nx.adjacency_matrix(G).todense())
            num_nodes_in_net = len(A[0])
            Beta = np.full((num_nodes_in_net, num_nodes_in_net), beta_init)
            Gamma = np.full(num_nodes_in_net, gamma)
        results = simulate(A, Beta, Gamma, i, intervention_gen, beta_interv, prop_reduced)
        for gen in range(len(results)):
            num_total_infctd = int(results[gen])
            try:
                outbreak_size_distrb_per_gen_matrix[gen][num_total_infctd] += 1
            except IndexError:
                print('Index error for g: ', gen, ', gen_s: ', num_total_infctd)
                continue
    # averaging:
    for gen in range(100):
        gen_time_series = outbreak_size_distrb_per_gen_matrix[gen]
        gen_time_series = gen_time_series / num_sims
        outbreak_size_distrb_per_gen_matrix[gen] = gen_time_series
    return outbreak_size_distrb_per_gen_matrix


def simulate(A, Beta, Gamma, current, intervention_gen=-1, beta_interv=-1.0, prop_reduced=0.0):
    start_time = time.time()
    # With intervention into the simulation code
    if intervention_gen > 0:
        # Example of random vaccination:
        sim = extended_simulation.RandomInterventionSim(A)
        sim.add_infection_event_rates(Beta)
        sim.add_recover_event_rates(Gamma)
        sim.configure_intervention(intervention_gen=intervention_gen, beta_redux=beta_interv,
                                   proportion_reduced=prop_reduced)
    else:
        sim = simulation.Simulation(A)
        sim.add_infection_event_rates(Beta)
        sim.add_recover_event_rates(Gamma)
    sim.run_sim()
    if current % 50 == 0:
        print('current sim ' + str(current))
        print("--- %s seconds to run simulation---" % (time.time() - start_time))
    results = sim.tabulate_generation_results(100)
    return results
