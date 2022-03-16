# EpIntervene
EpIntervene is a custom simulation framework in Python for Event-Driven simulation of epidemics on networks. 
What makes this framework unique is its functionality for introducing interventions in the middle of epidemic simulations, its customizable event classes, and its book-keeping of generations of infection.

EpIntervene is still a work in progress. The current version of the software is available for installation
on the PyPi index, available [here](https://pypi.org/project/epintervene/), available via a recommended `pip` install in an `anacondas` virtual environment. 

## Usage
For quick examples on how to use EpIntervene, we recommend walking through the examples in the `examples/sandbox.py` file [here](https://github.com/andrea-allen/epintervene/blob/main/examples/sandbox.py)
for standard demonstrations on SIR and SEIR epidemics. For the intervention features, you can walk through the examples
in `examples/interventions_sandbox.py` found [here.](https://github.com/andrea-allen/epintervene/blob/main/examples/interventions_sandbox.py)
For an interactive tool to let you explore what the results of using EpIntervene would look like,
check out my [Streamlit app](https://share.streamlit.io/andrea-allen/epintervene-interactive/main/app.py)
with source code found [here.](https://github.com/andrea-allen/epintervene-interactive)

This package and its current capabilities support SIR or SEIR model simulations, and purely SIR for the intervention-based simulations.
The two primary Simulation objects can be found in the classes `epintervene.simobjects.simulation.Simulation`
and `epintervene.simobjects.simulation.SimulationSEIR`. For the intervention feature, use any of the available Simulation classes in
`epintervene.simobjects.extended_simulation.py` once you have familiarized yourself with how to run a standard simulation.

Currently, calls to simulators `run_sim()` methods with the appropriate configurations and run time arguments
result in a single simulation. It is up to the user to aggregate results from an ensemble of runs on their own.

For sample usage, see the `examples` module. 

### Setting up a Simulation Object 
A single simulation object is accessed and created via `simobjects.simulation.Simulation`
and consists of an object encapsulating all information about the configuration, state, and results
for a single simulation. The base class, `Simulation`, will run an SIR disease model. 

#### Basic Configurations
To configure a `Simulation` object, you'll first need to define a network via a symmetric adjacency list.
The `network` module has a built-in function to assist with this and is recommended in order
to obtain an adjacency list in the right form. If supplying one's own adjacency list, make sure it is in the form
``` 
[line of node_label_index], node_label, neighbor, neighbor, ...
 [line 0]                   0           1         7
 [line 1]                   1           0         3
 [lines ...                 ...         ...       ...
 [line 3]                   3           1         22         26   ...
```

To use the built-in functions, supply a networkX object of a graph `G` and call 
``` adjlist = network.NetworkBuilder.create_adjacency_list(G)```
This will remove multi-edges, and make a symmetric adjacency list of the form displayed above. `G` must be
a NetWorkX Graph object. 

Then to create a simulation object, define
```your_simulation_object = simulation.Simulation(adj_list=adjlist, N=len(adjlist))```
where `N` is the number of nodes in your network and should be the length of the list.

Then add your `beta` and `gamma` rates for infection rate and recovery rate to your simulation via
```
        your_simulation_object.set_uniform_beta(beta=0.004)
        your_simulation_object.set_uniform_gamma(gamma=0.001)
```

#### Additional Configurations
It may be the case that your network model is partitioned into sub-populations, e.g. a stochastic block model network. Or, your network
is based on real data for which the labels indicate specific data on a node's sub-group within the network. 

For this reason, you can specify `membership_groups` and `node_memberships` when initializing a new `Simulation` object. 
Simply provide a list of ids or labels for each unique membership group, which can be any simple data type (string, its, floats, etc)
and then supply a vector (1-d numpy array) of the membership id for each node. The position in the vector will be read
as the node id, and the entry will indicate the membership group. For example, for a network with 10 nodes representing a zoo, with adjacency list `adjlist`, we would set up 

```zoo_sim = Simulation(adj_list=adjlist, N=10, membership_groups=['elephant', 'tiger', 'bird'], node_memberships=['tiger', 'tiger', 'bird', 'bird', 'bird', ..., 'elephant', 'tiger']```


#### Intervention Models and Configurations
This package offers a few types of simulations with customizable intervention regimes. To access
the intervention simulation objects, access via `simobjects.extended_simulation`.

`UniversalInterventionSim(Simulation)` will reduce `beta` for the entire network to the specified `beta_redux` at the specified 
intervention generation. How to configure:
```
adjlist = nb.create_adjacency_list(G)
sim = extended_simulation.UniversalInterventionSim(N=len(adjlist), adjlist=adjlist)
sim.set_uniform_beta(0.9)
sim.set_uniform_gamma(0.1)
sim.configure_intervention(intervention_gen=4, beta_redux=0.6)
```
This means that when the epidemic reaches generation `4`, `beta` the transmission rate will be reduced to `beta_redux`.

`RandomInterventionSim(Simulation)`
This simulation class lets the user configure an intervention at a given epidemic generation to model vaccinating a specified proportion
of the population to a transmissibility of zero, by selecting a random subset of the population. See the
`examples/interventions_sandbox` module for examples. 

Other interventions: 

`TargetedInterventionSim`, which works just as `RandomInterventionSim` but where the vaccination strategy targets the specified
proportion of the population in decreasing order of degree; highest degree nodes are vaccinated first. *Note: Under this
simulation framework, a desired percentage of nodes to vaccinate is specified, however the simulation will
fill that quota regardless of current vaccination or infection status of the node.*

`RandomRolloutSimulation` and `TargetedRolloutSimulation` allow the user to configure a phased rollout of either
random or targeted vaccinations. Specify the list of generations to intervene at, and the proportion of the population
that should be vaccinated at each generation.

`AbsoluteTimeNetworkSwitchSim` allows the user to switch the network on the same nodes at a specified time of intervention.
It is best to familiarize with oneself of the simulation on both networks first, to get a sense of the appropriate
time to intervene based on infection rates.

For all of the above, examples are provided in `examples/interventions_sandbox.py`.



#### SEIR Model
An SEIR model simulation object can be set up and configured much in the same way as the SIR one.
Currently, interventions are not yet supported on the SEIR framework, but this feature is coming soon.
The SEIR simulation framework does support membership groups.

The additional configurations needed from the user are additional rate matrices to describe
the node-pair rate for Exposed-Susceptible transmission, and a rate specifying the transition rate
for individual nodes from Exposed to Infected. 
Example:
```
seir_sim = SimulationSEIR(N=N, adjlist=adj_list)
seir_sim.set_uniform_gamma(gamma=0.00001)
seir_sim.set_uniform_beta(beta=0.5)
seir_sim.set_uniform_beta_es(beta_es=0.5)
seir_sim.set_uniform_gamma_ei(gamma_ei=1.0)
```

Similarly, one can specify `node_memberships=node_membership_vector, membership_groups=['tiger', 'bird', 'elephant']` with the SEIR model.

### Running a Simulation object
Once your Simulation object is set up, you are ready to run it. 

An important feature of EpIntervene is that even after running your simulation, the object will
be preserved with all of its state and attributes, and can be accessed repeatedly for different 
types of results. 
#### If you run a simulation again, you SHOULD create a NEW Simulation object in order to guarantee a clean state.

To run any of the simulations you've configured, call
```
your_sim_object.run_sim()
```
Calling `run_sim()` without changing any of its default arguments is the simplest way to run a single simulation
of your object. If you want to specify certain aspects of the simulation, there are a range of options you can choose from.

As of now, `run_sim(uniform_rate=True)` is the default. Non-uniform rates are no longer supported so there
is no need to specify a value for it.

If you want to track membership groups and you have configured them properly, you can run
```
your_sim_object.run_sim(with_memberships=True)
```

To stop the simulation when there are no more possible infection events, you can specify
```
your_sim_object.run_sim(wait_for_recovery=False)
```

Other optional arguments:
`kill_by=13` will stop the simulation after Generation 13 when there are no more active members of any generation from 0-13.

`p_zero=i` specifies the node label/index for patient zero. If not specified, a random node index will be chosen. To specify multiple
patient zeros, for user-defined node labels set
`p_zero=[i, j, k]` and to generate a specific amount (greater than 1) random patient zeros, specify
`p_zero=[None, None, None, None]` which will generate 4 random patient zero nodes. To specify a mix of random and non-random patient zeros,
specify `p_zero=[i, j, None, k, None]` in any order. If two copies of the same label are added, one copy will be replaced with another random patient zero.

`record_active_gen_sizes=True/False` a boolean, defaults to False. Will record the number and size of active generations over time.
An active node is defined as being infectious and still having susceptible neighbors. An active generation is defined as containing
one or more active nodes belonging to that epidemic generation.

### Obtaining results from a Simulation
#### Basic Results
The Epintervene framework is driven by an Event-Driven algorithm, in which
continuous time is tracked by drawing an exponential random variable for the waiting time until next event at each discrete time step
with rate parameter the sum of the rates of all potential individual events. As a result, the resulting
raw time series will not be normalized, making it hard to compare to an ensemble of results should the user
intend to average the results of an ensemble. Therefore, when calling the Simulation object for
continuous time results, results are binned into a time-resolution that is evenly split into bins.

To obtain a continuous-time time series of your infected and recovered nodes for an SIR simulation,
after you have run your simulation call
```
ts, infect_ts, recover_ts = sim.tabulate_continuous_time(time_buckets=1000)
```
which will return a time series time values (binned evenly between the start time and end time of the simulation at the resolution of the `time_buckets` optional parameter)
and the time series of infected and recovered nodes.
For the SEIR model, the same call returns four time series including the time values themselves, for example
```
ts, infect_ts, recover_ts, exposed_ts = seir_sim.tabulate_continuous_time(time_buckets=1000)
```

#### Tracking Generations of Infection 
EpIntervene also provides the user the option to obtain a generational time series from a
Simulation object. The form of this return value will be a single vector, where the index of the vector indicates the epidemic generation, in ascending order. 
Each entry corresponds to the number of cumulatively infected nodes
belonging to a generation that is less than or equal to the entry's index. For example, 
```
ts_by_gen = sim.tabulate_generation_results(max_gens=20)
```
might look like `[1, 2, 4, 8, 13, 18, 37]`, which indicates that 1 node was infected in generation 0,
2 nodes were infected in generation 0 and 1 combined, and 13 nodes were infected including generations 0,1, 2, 3, and 4. 
It is worth noting for the user that, to obtain the number of nodes in each unique generation,
one can call a `diff()` method on the result, where `ts_by_gen[i+1] - ts_by_gen[i]` will return the
number of nodes infected in generation `i+1`.

This is not a traditional time series in this form, it is a chronological list of cumulative infections
by the birth of each epidemic generation. It does not formally correspond to continuous time results.


#### Tracking results for population groups
To obtain distinct time series results for each membership group, if specified in the simulation,
the user may call
```
ts, membership_ts_infc = sim.tabulate_continuous_time_with_groups(time_buckets=1000)
```
which will return the time series values, and a dictionary of time series of infections per group.
An example of how to iterate through the results and plot each group is here:
```
for group in membership_ts_infc.keys():
    plt.plot(ts, membership_ts_infc[group], label=group)
```

For the SEIR model, there are 3 return values to include the exposed time series. Example:
```
ts, membership_ts_infc, membership_ts_exp = seir_sim.tabulate_continuous_time_with_groups(time_buckets=1000)
``` 


## Classes
See the docs directory for more detailed documentation on each class. Documentation will be updated
and is a living document. 

### Network Objects
#### NetworkBuilder
A static class to assist users in creating `NetworkX` graphs. The NetworkBuilder can return a network
built from a degree distribution, by calling `from_degree_distribution()` or from an adjacency matrix,
by calling `from_adjacency_matrix()`
#### Node
Internal object that represents a single node in the network, created and disposed during simulations.
#### Edge
Internal object that represents a pair of nodes, where the left node and right node are distinct.

## Examples
See module `examples.sandbox.py` for examples of basic SIR, SEIR, and group membership simulations.

See module `examples.interventions_sandbox.py` for examples of intervention-based simulations.

Results of an SIR simulation (single run) with membership groups:

![SEIR example](https://github.com/andrea-allen/epintervene/blob/main/docs/SEIR_group_example.png?raw=true)

![SEIR example](https://github.com/andrea-allen/epintervene/blob/main/docs/SEIR_total.png?raw=true)

![network](https://github.com/andrea-allen/epintervene/blob/main/docs/SEIR_network.png?raw=true)

![netswitch example](https://github.com/andrea-allen/epintervene/blob/main/docs/netswitch_samplefig.png?raw=true)

## Citing and using this code

If using all or any part of this code base in a published or non-published work, please cite as:

Cite this work as:

Andrea Allen. (2021). andrea-allen/epintervene: v1.0.3 EpIntervene Release (v1.0.3). Zenodo. https://doi.org/10.5281/zenodo.5514401