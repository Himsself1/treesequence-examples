'''
Simple model of 2 populations
'''
from pathlib import Path
import msprime
import tskit
import random
import time

def simulation_loop( pop1_size, pop2_size, ancestral_size, split_time, seed ):
    '''
    Loop that performs simulations.
    One can also use the 'num_replicates' argument but I find it
    easier to do it in a loop that saves the output everytime.
    '''
    ### Building demographic scenario
    demography = msprime.Demography()
    demography.add_population( name="ANCESTRAL",
                               initial_size = ancestral_size)
    demography.add_population(name="Pop1",
                              initial_size = pop1_size,
                              default_sampling_time = 0 )
    demography.add_population(name="Pop2",
                              initial_size = pop2_size,
                              default_sampling_time = 0 )
    demography.add_population_split(time = split_time,
                                    derived = ["Pop1", "Pop2"],
                                    ancestral = "ANCESTRAL" )
    ## ALWAYS RUN "sort_events" starting the simulation!
    demography.sort_events()
    ts = msprime.sim_ancestry(
        demography = demography,
        samples = {"Pop1" : 5, "Pop2" : 5},  ## Specify number of samples from each population
        recombination_rate = 1e-8,
        sequence_length = 1e+7,
        random_seed = seed
    ) ## Command that starts simulation. The output is a tree-sequence
    return ts

## Run a loop that saves treesequences in files
reps = 10 ## Number of repetition
for i in range(reps):
    temp_ts = simulation_loop( 1000, 1000, 1000, 1000, 50)
    temp_ts.dump( "path/to/save/location/name".join(i) )
    
