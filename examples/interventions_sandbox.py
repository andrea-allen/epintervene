import numpy as np
from epintervene.simobjects import network
import matplotlib.pyplot as plt
from epintervene.simobjects import simulation, extended_simulation
import math

def uniform_intervention_example():
    print('Uniform intervention simulation tutorial')
    # # **** TIP ***** For best demonstration, if time series results show a simulation in which only a single node
    # # was infected, try running the example again to see another result.

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
    my_simulation = extended_simulation.UniversalInterventionSim(adjlist=adjlist, N=len(adjlist))

    # # Setting required configurations
    my_simulation.set_uniform_beta(beta=0.9) # .9 people per day per infected person
    my_simulation.set_uniform_gamma(gamma=0.2) # 5 days to recover

    # # CONFIGURE THE INTERVENTION SPECIFICATION
    # # Here, this configuration means that at generation 4 (i.e. when the first person belonging to generation 4
    # # becomes infected), the entire population's rate of infection is reduced to 0.2
    my_simulation.configure_intervention(intervention_gen=4, beta_redux=0.2)

    # # Setting up a simulation object with the original parameters and no intervention:
    regular_simulation = simulation.Simulation(adj_list=adjlist, N=len(adjlist))
    regular_simulation.set_uniform_beta(beta=0.9)
    regular_simulation.set_uniform_gamma(gamma=0.2)

    # # Running the simulations
    my_simulation.run_sim()
    regular_simulation.run_sim()

    # # Obtaining time series results
    time_series_vals, infected_time_series, recovered_time_series = my_simulation.tabulate_continuous_time(time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    time_series_vals_reg, infected_time_series_reg, recovered_time_series_reg = regular_simulation.tabulate_continuous_time(time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected w/ intervention')
    plt.plot(time_series_vals_reg, infected_time_series_reg, label='infected normally')
    plt.plot(time_series_vals, recovered_time_series, label='recovered w/ intervention')
    plt.plot(time_series_vals_reg, recovered_time_series_reg, label='recovered normally')
    plt.legend(loc='upper left')
    plt.title('Auto generated time series values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining cumulative infections by generations of infection:
    gen_results = my_simulation.tabulate_generation_results(max_gens=30)
    gen_results_reg = regular_simulation.tabulate_generation_results(max_gens=30)

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results, label='cumulative infections w/ intervention')
    plt.scatter(np.arange(30), gen_results_reg, label='normal cumulative infections')
    plt.plot(np.arange(30), gen_results)
    plt.plot(np.arange(30), gen_results_reg)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.legend(loc='lower right')
    plt.show()

def random_intervention_example():
    print('Random intervention simulation tutorial')
    # # **** TIP ***** For best demonstration, if time series results show a simulation in which only a single node
    # # was infected, try running the example again to see another result.

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
    my_simulation = extended_simulation.RandomInterventionSim(adjlist=adjlist, N=len(adjlist))

    # # Setting required configurations
    my_simulation.set_uniform_beta(beta=0.9)  # .9 people per day per infected person
    my_simulation.set_uniform_gamma(gamma=0.2)  # 5 days to recover

    # # CONFIGURE THE INTERVENTION SPECIFICATION
    # # Here, this configuration means that at generation 4 (i.e. when the first person belonging to generation 4
    # # becomes infected), a random set of nodes making up 30% of the population will be reduced to 0 transmissibility
    # # For now, only beta=0 is supported for Random vaccination.
    my_simulation.configure_intervention(intervention_gen=4, proportion_reduced=0.3, beta_redux=0)

    # # Setting up a simulation object with the original parameters and no intervention:
    regular_simulation = simulation.Simulation(adj_list=adjlist, N=len(adjlist))
    regular_simulation.set_uniform_beta(beta=0.9)
    regular_simulation.set_uniform_gamma(gamma=0.2)

    # # Running the simulations
    my_simulation.run_sim()
    regular_simulation.run_sim()

    # # Obtaining time series results
    time_series_vals, infected_time_series, recovered_time_series = my_simulation.tabulate_continuous_time(
        time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    time_series_vals_reg, infected_time_series_reg, recovered_time_series_reg = regular_simulation.tabulate_continuous_time(
        time_buckets=1000,
        custom_range=True,
        custom_t_lim=15)
    # # Plotting the results:
    plt.plot(time_series_vals, infected_time_series, label='infected w/ intervention')
    plt.plot(time_series_vals_reg, infected_time_series_reg, label='infected normally')
    plt.plot(time_series_vals, recovered_time_series, label='recovered w/ intervention')
    plt.plot(time_series_vals_reg, recovered_time_series_reg, label='recovered normally')
    plt.legend(loc='upper left')
    plt.title('Auto generated time series values')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining cumulative infections by generations of infection:
    gen_results = my_simulation.tabulate_generation_results(max_gens=30)
    gen_results_reg = regular_simulation.tabulate_generation_results(max_gens=30)

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results, label='cumulative infections w/ intervention')
    plt.scatter(np.arange(30), gen_results_reg, label='normal cumulative infections')
    plt.plot(np.arange(30), gen_results)
    plt.plot(np.arange(30), gen_results_reg)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.legend(loc='lower right')
    plt.show()

def targeted_intervention_example():
    print('Targeted intervention simulation tutorial')

    # # In this example, we will run each simulation 10 times and take the ensemble average of the effects, to obtain
    # # a smoother picture of the type of results available.

    # # Create a NetworkBuilder object, to help with creating and configuring the network for the simulation:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 3)
    # # Generating a network from the above degree distribution, with 100 nodes
    G, pos = nb.from_degree_distribution(100, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist = nb.create_adjacency_list(G)

    # # Initializing the data structures for the ensemble of results:
    time_series_vals = None
    infected_ts_normal = None
    infected_ts_random = None
    infected_ts_targeted = None
    recovered_ts_normal = None
    recovered_ts_random = None
    recovered_ts_targeted = None

    # # Looping as many times as the number of simulations we want to run, we will initialize a NEW simulation object
    # # every time, and collate the results:
    for n in range(40):
        # # Constructing the simulation object
        my_simulation = extended_simulation.TargetedInterventionSim(adjlist=adjlist, N=len(adjlist))

        # # Setting required configurations
        my_simulation.set_uniform_beta(beta=0.9)  # .9 people per day per infected person
        my_simulation.set_uniform_gamma(gamma=0.2)  # 5 days to recover

        # # CONFIGURE THE INTERVENTION SPECIFICATION
        # # Here, this configuration means that at generation 4 (i.e. when the first person belonging to generation 4
        # # becomes infected), a set of nodes making up 30% of the population will be reduced to 0 transmissibility.
        # # The nodes are chosen in order of degree class, highest degree nodes included first, until 30% quota is reached.
        # # For now, only beta=0 is supported for Targeted vaccination.
        my_simulation.configure_intervention(intervention_gen=4, proportion_reduced=0.3, beta_redux=0)

        # # Setting up a simulation object with the original parameters and no intervention:
        regular_simulation = simulation.Simulation(adj_list=adjlist, N=len(adjlist))
        regular_simulation.set_uniform_beta(beta=0.9)
        regular_simulation.set_uniform_gamma(gamma=0.2)

        # # Setting up a simulation with Random vaccination (from example above) for comparison purposes.
        random_simulation = extended_simulation.RandomInterventionSim(adjlist=adjlist, N=len(adjlist))
        random_simulation.set_uniform_beta(beta=0.9)
        random_simulation.set_uniform_gamma(gamma=0.2)
        random_simulation.configure_intervention(intervention_gen=4, proportion_reduced=0.3, beta_redux=0)

        # # Running the simulations
        my_simulation.run_sim()
        regular_simulation.run_sim()
        random_simulation.run_sim()

        # # Obtaining time series results
        if time_series_vals is None:
            time_series_vals, infected_ts_targeted, recovered_ts_targeted = my_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            time_series_vals, infected_ts_normal, recovered_ts_normal = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            time_series_vals, infected_ts_random, recovered_ts_random = random_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
        else:
            _, i_t_t, r_t_t = my_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_targeted += i_t_t
            recovered_ts_targeted += r_t_t
            _, i_t_n, r_t_n = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_normal += i_t_n
            recovered_ts_normal += r_t_n
            _, i_t_r, r_t_r = random_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_random += i_t_r
            recovered_ts_random += r_t_r

    # # Normalizing for the 10 simulations
    infected_ts_normal = infected_ts_normal / 10
    infected_ts_random = infected_ts_random / 10
    infected_ts_targeted = infected_ts_targeted / 10
    recovered_ts_normal = recovered_ts_normal / 10
    recovered_ts_random = recovered_ts_random / 10
    recovered_ts_targeted = recovered_ts_targeted / 10

    # # Plotting the results:
    plt.plot(time_series_vals, infected_ts_targeted, label='infected w/ targeted intervention', color='blue')
    plt.plot(time_series_vals, infected_ts_random, label='infected w/ random intervention', color='red')
    plt.plot(time_series_vals, infected_ts_normal, label='infected normally', color='orange')
    plt.plot(time_series_vals, recovered_ts_targeted, label='recovered w/ targeted intervention', color='blue', ls='--')
    plt.plot(time_series_vals, recovered_ts_normal, label='recovered normally', color='orange', ls='--')
    plt.plot(time_series_vals, recovered_ts_random, label='recovered w/ random intervention', color='red', ls='--')
    plt.legend(loc='upper left')
    plt.title('Effects of Targeted Vaccination compared \n to normal and random')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Obtaining cumulative infections by generations of infection:
    # # These results in this example are not the ensemble result, just the results from a single (the latest) run
    gen_results = my_simulation.tabulate_generation_results(max_gens=30)
    gen_results_reg = regular_simulation.tabulate_generation_results(max_gens=30)
    gen_results_rand = random_simulation.tabulate_generation_results(max_gens=30)

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results, label='cumulative infections w/ targeted intervention')
    plt.scatter(np.arange(30), gen_results_reg, label='normal cumulative infections')
    plt.scatter(np.arange(30), gen_results_rand, label='cumulative infections w/ random intervention')
    plt.plot(np.arange(30), gen_results)
    plt.plot(np.arange(30), gen_results_reg)
    plt.plot(np.arange(30), gen_results_rand)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.legend(loc='lower right')
    plt.show()

def absolute_time_intervention_example():
    # # In this example, we will run each simulation 100 times and take the ensemble average of the effects, to obtain
    # # a smoother picture of the type of results available.

    # # This intervention allows you to switch the population network at a designated time during the run.
    # It is recommended that you first get familiar with the dynamics of the simulation on both networks, and the
    # relative time dynamics of t, so that you can provide a proper time value t for the intervention.

    # # Create a NetworkBuilder object, to help with creating and configuring the network for the simulation:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 3)
    # # Generating a network from the above degree distribution, with 100 nodes
    G, pos = nb.from_degree_distribution(100, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist_1 = nb.create_adjacency_list(G)

    # # Create a second network, for the intervention where we will switch networks:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 7)
    # # Generating a network from the above degree distribution, with 100 nodes
    G, pos = nb.from_degree_distribution(100, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist_2 = nb.create_adjacency_list(G)

    # # Initializing the data structures for the ensemble of results:
    time_series_vals = None
    infected_ts_normal = None
    infected_ts_netswitch = None
    recovered_ts_normal = None
    recovered_ts_netswitch = None

    # # Looping as many times as the number of simulations we want to run, we will initialize a NEW simulation object
    # # every time, and collate the results:
    num_runs = 100
    for n in range(num_runs):
        # # Constructing the simulation object
        my_simulation = extended_simulation.AbsoluteTimeNetworkSwitchSim(adjlist=adjlist_2, N=len(adjlist_1))

        # # Setting required configurations
        my_simulation.set_uniform_beta(beta=0.1)  # .7 people per day per infected person
        my_simulation.set_uniform_gamma(gamma=0.005)  # 5 days to recover

        # # CONFIGURE THE INTERVENTION SPECIFICATION
        # # Here, this configuration means that at time t=7, the population network will switch to that given by
        # adjlist_2.
        my_simulation.configure_intervention(intervention_time=7, new_adjlist=adjlist_1)

        # # Setting up a simulation object with the original parameters and no intervention, that only runs on
        # adjlist_1's network
        regular_simulation = simulation.Simulation(adj_list=adjlist_2, N=len(adjlist_1))
        regular_simulation.set_uniform_beta(beta=0.1)
        regular_simulation.set_uniform_gamma(gamma=0.005)

        # # Running the simulations
        my_simulation.run_sim()
        regular_simulation.run_sim()

        # # Obtaining time series results
        if time_series_vals is None:
            time_series_vals, infected_ts_netswitch, recovered_ts_netswitch = my_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            time_series_vals, infected_ts_normal, recovered_ts_normal = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)

        else:
            _, i_t_t, r_t_t = my_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_netswitch += i_t_t
            recovered_ts_netswitch += r_t_t
            _, i_t_n, r_t_n = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_normal += i_t_n
            recovered_ts_normal += r_t_n

    # # Normalizing for the 10 simulations
    infected_ts_normal = infected_ts_normal / num_runs
    infected_ts_netswitch = infected_ts_netswitch / num_runs
    recovered_ts_normal = recovered_ts_normal / num_runs
    recovered_ts_netswitch = recovered_ts_netswitch / num_runs

    # # Plotting the results:
    plt.plot(time_series_vals, infected_ts_netswitch, label='infected w/ net switch intervention', color='blue')
    plt.plot(time_series_vals, infected_ts_normal, label='infected normally', color='orange')
    plt.plot(time_series_vals, recovered_ts_netswitch, label='recovered w/ net switch intervention', color='blue', ls='--')
    plt.plot(time_series_vals, recovered_ts_normal, label='recovered normally', color='orange', ls='--')
    plt.legend(loc='upper left')
    plt.title('Effects of Network Switching at t=7 compared \n to normal run on Network 1')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

def rollout_intervention_examples():
    # # If you want to simulate a series of vaccinations, in which a progressive fraction of the population is
    #   vaccinated in either the Random or Targeted schemes, follow the below examples:
    # # In this example, we will run each simulation 10 times and take the ensemble average of the effects, to obtain
    # # a smoother picture of the type of results available.

    # # Create a NetworkBuilder object, to help with creating and configuring the network for the simulation:
    nb = network.NetworkBuilder
    # # Specify a degree distribution to create a configuration model network if desired:
    degree_distrb = binomial_degree_distb(400, 3)
    # # Generating a network from the above degree distribution, with 100 nodes
    G, pos = nb.from_degree_distribution(300, degree_distrb)
    # # If you already have a network, feel free to skip the above steps and generate your own NetworkX object,
    # # or skip this next step and provide a symmetric adjacency list.
    adjlist = nb.create_adjacency_list(G)

    # # Initializing the data structures for the ensemble of results:
    time_series_vals = None
    infected_ts_normal = None
    infected_ts_random_rollout = None
    infected_ts_targeted_rollout = None

    gen_results_reg = np.zeros(30)
    gen_results_targ = np.zeros(30)
    gen_results_rand = np.zeros(30)

    # # Looping as many times as the number of simulations we want to run, we will initialize a NEW simulation object
    # # every time, and collate the results:
    num_sims = 100
    for n in range(num_sims):
        # # Constructing the simulation object, a rollout with random vaccination
        random_rollout_simulation = extended_simulation.RandomRolloutSimulation(adjlist=adjlist, N=len(adjlist))

        # # Setting required configurations, for starting
        random_rollout_simulation.set_uniform_beta(beta=0.9)  # .9 people per day per infected person
        random_rollout_simulation.set_uniform_gamma(gamma=0.2)  # 5 days to recover

        # # Constructing a targeted rollout intervention object
        targeted_rollout_simulation = extended_simulation.TargetedRolloutSimulation(adjlist=adjlist, N=len(adjlist))
        targeted_rollout_simulation.set_uniform_beta(beta=0.9)
        targeted_rollout_simulation.set_uniform_gamma(gamma=0.2)

        # # CONFIGURE THE INTERVENTION SPECIFICATIONS
        # # Here we will configure the interventions for both random and targeted rollouts.
        # # We will define a list of which generations to intervene at, and what proportion of the population to
        # vaccinate (with 100% efficacy) at each successive generation.
        # # The ONLY DIFFERENCE between the random vs targeted scheme is the vaccination STRATEGY, for
        # Random Rollout, a random set of the designated percentage of individuals is selected.
        # # For Targeted Rollout, the designated percentage is selected by decreasing order of highest degree nodes.

        random_rollout_simulation.configure_intervention(intervention_gen_list=[3, 5, 7], beta_redux_list=[0,0,0], proportion_reduced_list=[.02, .03, .05])
        targeted_rollout_simulation.configure_intervention(intervention_gen_list=[3, 5, 7], beta_redux_list=[0,0,0], proportion_reduced_list=[.02, .03, .05])
        # # Notice how in this example, we are configuring the random and targeted interventions with the same
        # configurations, in order to observe the difference that the two strategies have compared to one another.

        # # Setting up a simulation object with the original parameters and no intervention:
        regular_simulation = simulation.Simulation(adj_list=adjlist, N=len(adjlist))
        regular_simulation.set_uniform_beta(beta=0.9)
        regular_simulation.set_uniform_gamma(gamma=0.2)

        # # Running the simulations
        regular_simulation.run_sim()
        random_rollout_simulation.run_sim()
        targeted_rollout_simulation.run_sim()

        # # Obtaining time series results
        if time_series_vals is None:
            time_series_vals, infected_ts_normal, _ = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            time_series_vals, infected_ts_random_rollout, _ = random_rollout_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            time_series_vals, infected_ts_targeted_rollout, _ = targeted_rollout_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            # # Obtaining cumulative infections by generations of infection:
            gen_results_reg += regular_simulation.tabulate_generation_results(max_gens=30)
            gen_results_targ += targeted_rollout_simulation.tabulate_generation_results(max_gens=30)
            gen_results_rand += random_rollout_simulation.tabulate_generation_results(max_gens=30)
        else:
            _, i_t_t, _ = regular_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_normal += i_t_t

            _, i_t_n, _ = random_rollout_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_random_rollout += i_t_n

            _, i_t_r, _ = targeted_rollout_simulation.tabulate_continuous_time(
                time_buckets=1000,
                custom_range=True,
                custom_t_lim=15)
            infected_ts_targeted_rollout += i_t_r

            # # Obtaining cumulative infections by generations of infection:
            gen_results_reg += regular_simulation.tabulate_generation_results(max_gens=30)
            gen_results_targ += targeted_rollout_simulation.tabulate_generation_results(max_gens=30)
            gen_results_rand += random_rollout_simulation.tabulate_generation_results(max_gens=30)

    # # Normalizing for the 10 simulations
    infected_ts_normal = infected_ts_normal / num_sims
    infected_ts_random_rollout = infected_ts_random_rollout / num_sims
    infected_ts_targeted_rollout = infected_ts_targeted_rollout / num_sims

    # # Plotting the results:
    plt.plot(time_series_vals, infected_ts_targeted_rollout, label='infected w/ targeted rollout intervention', color='blue')
    plt.plot(time_series_vals, infected_ts_random_rollout, label='infected w/ random rollout intervention', color='red')
    plt.plot(time_series_vals, infected_ts_normal, label='infected normally', color='orange')

    plt.legend(loc='upper left')
    plt.title('Effects of Targeted Vaccination compared \n to normal and random')
    plt.xlabel('Time')
    plt.ylabel('Number of nodes')
    plt.show()

    # # Plotting generation results
    plt.scatter(np.arange(30), gen_results_targ/num_sims, label='cumulative infections w/ targeted rollout intervention')
    plt.scatter(np.arange(30), gen_results_reg/num_sims, label='normal cumulative infections')
    plt.scatter(np.arange(30), gen_results_rand/num_sims, label='cumulative infections w/ random rollout intervention')
    plt.plot(np.arange(30), gen_results_targ/num_sims)
    plt.plot(np.arange(30), gen_results_reg/num_sims)
    plt.plot(np.arange(30), gen_results_rand/num_sims)
    plt.title('Cumulative infections by epidemic generations')
    plt.xlabel('Generation number')
    plt.ylabel('Cumulative Infections')
    plt.legend(loc='lower right')
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


if __name__ == '__main__':
    targeted_intervention_example()
    uniform_intervention_example()
    random_intervention_example()
    # absolute_time_intervention_example()
    rollout_intervention_examples()
