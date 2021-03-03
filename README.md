# EpIntervene
Custom simulation framework in Python for Event-Driven simulation of epidemics on networks. What makes this framework unique is its functionality for introducing interventions in the middle of epidemic simulations, its customizable event classes, and its book-keeping of generations of infection.

EpIntervene is still a work in progress. Please check back later.

## Usage
This package is meant for use for SIR or SEIR simulations. Currently, calls to simulators `run()` methods
result in a single simulation. It is up to the user to aggregate ensemble results on their own. 

### Setting up a Simulation Object 
A single simulation object is accessed and created via `simobjects.simulation.Simulation`
and consists of an object encapsulating all information about the configuration, state, and results
for a single simulation. The base class, `Simulation`, will run an SIR disease model. 

#### Basic Configurations
To create, give `Simulation(A)` an adjacency matrix in the form of a `numpy array` as the only required argument for construction.
You must also supply the simulation with a `Beta` matrix and `Gamma` vector that satisfies `Beta[i, j]` is the rate of infection
from node `i` to `j` (note that `Beta` need not be symmetric) and `Gamma[i]` is the recovery probability for node `i`.

Add these by calling the methods     `.add_infection_event_rates(Beta)` and
`.add_recover_event_rates(Gamma)` to your Simulation object.

#### Additional Configurations
It may be the case that your network model is partitioned into sub-populations, e.g. a stochastic block model network. Or, perhaps your network
is based on real data for which the labels indicate specific data on a node's sub-group within the network. 

For this reason, you can specify `membership_groups` and `node_memberships` when initializing a new `Simulation` object. 
Simply provide a list of ids or labels for each unique membership group, which can be any simple data type (string, its, floats, etc)
and then supply a vector (1-d numpy array) of the membership id for each node. The position in the vector will be read
as the node id, and the entry will indicate the membership group. For example, for a network with 10 nodes representing a zoo, with adjacency matrix `A`, we would set up 

```zoo_sim = Simulation(A, membership_groups=['elephant', 'tiger', 'bird'], node_memberships=['tiger', 'tiger', 'bird', 'bird', 'bird', ..., 'elephant', 'tiger']```


#### Intervention Models and Configurations
This package offers a few types of simulations with customizable intervention regimes. To access
the intervention simulation objects, access via `simobjects.extended_simulation`.

`RandomInterventionSim(Simulation)`
Random Intervention is an extension of a the base `SIR` simulation in which a random set of nodes according to a specified proportion
of the network have their transmission probabilities (their entries in the `Beta` matrix) reduced to a specified value.
To initialize, follow the same steps as outlined in the section for the base `Simulation` class.
The only additional configuration necessary is to configure the intervention. One needs to specify the generation of intervention,
which will be triggered when the first node corresponding to that generation becomes infected. 
User will also specify the `beta_redux` reduction in transmissibility for the affected nodes, and the proportion of the network
that should have the intervention applied. Example:

```
extended_simulation.RandomInterventionSim(A)
sim.configure_intervention(intervention_gen=4, beta_redux=0.0, proportion_reduced=0.2)
```

Other interventions: 

`UniversalInterventionSim(Simulation)` will reduce `Beta` for the entire network to the specified `beta_redux` at the specified 
intervention generation. 

`MultiInterventionSim(Simulation)` will implement a phased rollout of multiple interventions,
which should be specified via a list of intervention generations and beta reduction values, and proportion
of network to be affected. Currently, MultiInterventionSim only supports Random intervention. 
Example:
```
sim = extended_simulation.MultiInterventionSim(A)
sim.add_infection_event_rates(Beta)
sim.add_recover_event_rates(Gamma)
sim.configure_intervention(intervention_gen_list=[3, 4, 5], beta_redux_list=[0.0, 0.0, 0.0], proportion_reduced_list=[0.10, 0.15, 0.15])
```


#### SEIR Model
An SEIR model simulation object can be set up and configured much in the same way as the SIR one.
Currently, interventions are not yet supported on the SEIR framework, but this feature is coming soon.

The additional configurations needed from the user are additional rate matrices to describe
the node-pair rate for Exposed-Susceptible transmission, and a network-length vector specifying the transition rate
for individual nodes from Exposed to Infected. 
Example:
```
seir_sim = SimulationSEIR(Simulation)
seir_sim.add_infection_event_rates(Beta_IS)
seir_sim.add_exposed_event_rates(Beta_ES)
seir_sim.add_recover_event_rates(Gamma)
seir_sim.add_exposed_infected_event_rates(Theta_EI)
```

### Running a Simulation object
Once your Simulation object is set up, you are ready to run it. 

An important feature of EpIntervene is that even after running your simulation, the object will
be preserved with all of its state and attributes, and can be accessed repeatedly for different 
types of results. 
#### If you run a simulation again, you SHOULD create a NEW Simulation objectin order to guarantee a clean state.

To run any of the simulations you've configured, call
```
your_sim_object.run_sim()
```
If you want to track membership groups and you have configured them, you can run
```
your_sim_object.run_sim(with_memberships=True)
```

### Obtaining results from a Simulation
#### Basic Results
Intrinsic to the EpIntervene framework is the Event-Driven simulation framework, in which
continuous time is tracked by drawing an exponential random variable for the time until next event at each discrete time step
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
An Edge also can double as an Event, equipped with an Event Rate (infection rate from left node to right susceptible node)


## Examples
Coming soon.

![alt text](https://github.com/andrea-allen/epintervene/blob/main/docs/sample_img.png?raw=true)
