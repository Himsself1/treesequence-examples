'''
Intermediate model of 4 populations.
'''
from pathlib import Path
import msprime
import tskit
import random
import time

def simulation_loop( pop1_size, pop2_size, pop3_size, pop4_size,
                     pop34_size, pop234_size, pop1234_size,
                     split_time_34, split_time_234, split_time_1234,
                     migration_proportion, seed ):
    '''
    Loop that performs simulations.
    One can also use the 'num_replicates' argument but I find it
    easier to do it in a loop that saves the output everytime.
    '''
    ## Building demographic scenario
    ## In msprime you need to treat internal branches as ancestral populations.
    ## So we need to specify
    demography = msprime.Demography()
    demography.add_population( name="anc_1234",
                               initial_size = pop1234_size)
    demography.add_population( name="anc_234",
                               initial_size = pop234_size)
    demography.add_population( name="anc_34",
                               initial_size = pop34_size)
    demography.add_population(name="Pop1",
                              initial_size = pop1_size,
                              default_sampling_time = 0 )
    demography.add_population(name="Pop2",
                              initial_size = pop2_size,
                              default_sampling_time = 0 )
    demography.add_population(name="Pop3",
                              initial_size = pop3_size,
                              default_sampling_time = 0 )
    demography.add_population(name="Pop4",
                              initial_size = pop4_size,
                              default_sampling_time = 0 )
    
    demography.add_population_split(time = split_time_34,
                                    derived = ["Pop3", "Pop4"],
                                    ancestral = "anc34" )
    demography.add_population_split(time = split_time_234,
                                    derived = ["Pop2", "anc_34"],
                                    ancestral = "anc_234" )
    demography.add_population_split(time = split_time_1234,
                                    derived = ["Pop1", "anc_234"],
                                    ancestral = "anc_1234" )
    ## Command to add migration. Note that migrations is backwards in time
    ## that means that source-destination is decided as you move backwards
    ## in terms of traversing a tree.
    demography.add_migration_rate_change(
        time = 0,
        source = "Pop4",
        dest = "Pop1",
        rate = migration_proportion
    )
    ## ALWAYS RUN "sort_events" starting the simulation!
    demography.sort_events()
    ts = msprime.sim_ancestry(
        demography = demography,
        samples = {"Pop1" : 5, "Pop2" : 5,
                   "Pop3" : 5, "Pop4" : 5},  ## Specify number of samples from each population. There are smarter ways to sample; see manual on sampling.
        recombination_rate = 1e-8,
        sequence_length = 1e+7,
        random_seed = seed
    ) ## Command that starts simulation. The output is a tree-sequence
    return ts

## Run a loop that saves treesequences in files
reps = 10 ## Number of repetition
for i in range(reps):
    ## 1st line in for loop is modern sizes
    ## 2nd line is ancestral sizes
    ## 3rd line is split times
    temp_ts = simulation_loop( 1000, 1000, 1000, 1000,
                               100, 100, 100,
                               1000, 2000, 3000,
                               0.05, 50 )
    temp_ts.dump( "path/to/save/location/name".join(i) )
