
description = '''

============//============ cs2Rotamer ============//============

Side chain chemical shifts to populations

The aim of this algorithm is to determine the rotamer populations of 
isoleucine from chemical shifts. 

The experimentally observed chemical shifts are listed in the 
file which by defult is called 'Measured_shifts.txt'.
The Calculated chemical shifts are given in XXX.

This program was written by Lucas Siemons. 
If you have any question or issues please contact me at 

Lucas.Siemons@googlemail.com

============//============             ============//============

'''

import os
import random as rn
import numpy as np

from cs2rotamer import fiIeIO as IO
from cs2rotamer import chemicalShifts as cs
import cs2rotamer.distanceFunctions as distFuncs


def First_Gen(Sim_number_iteration, params):
    
    '''
    This function creates the first generation. 
    It is based on the python random number generator random.random()
    
    Parameters
    ==========
    Sim_number_iteration : int
        simulaition iteration
    params : Parameters class 

    Returns
    =======
    individuals : dict 
        contains the individuals
    '''
    
    WorkingRandomList = []
    Random_Sum = 0.0
    RandomFloat = 0.0
    NormalisedList = []
    
    individuals = {}
    
    for indiv in range(params.pop_size):

        NormalisedList = []
        WorkingRandomList = []
        
        #This is a sum of the random numbers generated and is used to 
        #normalise the individual. 
        Random_Sum = 0
        for i in range(params.number_of_states):
            #here the first generation is made randomly 
            chance = rn.randint(0,100)
            if chance <= params.zeros:
        
                RandomFloat = 0.0000000001
            else:
                RandomFloat = rn.random()
            WorkingRandomList.append(RandomFloat )
            Random_Sum = Random_Sum + RandomFloat 
        
        #note here we normalise the vector so all the elements sum to 1 but this is the only 
        #time it is done. Later the random numbers are taken and they are normalised in
        #the fitness function when we compute the shifts in Computed_shifts
        for i in range(params.number_of_states):
            normalised_state_pop =  np.divide(WorkingRandomList[i],Random_Sum)
            NormalisedList.append(normalised_state_pop)
        
        individuals[indiv] = NormalisedList
    return individuals
    
def select_parents(daughter_generation, probabilities, params):
    

    '''This function aims to produce a list of 'parent couples' which are 
    then crossed over to make each individual
    
    Parameters
    ==========
    daughter_generation : dict 
        contains the next generation 
    
    probabilities : dict
        contains the 1/distance 

    params : Parameters class 

    Returns
    =======
    parent_list : list 
        list of parents
    '''
    
    sum = 0
    accumulated_reciprocals = []
    parent_list = []
    
    #this is used to do the weighted PDF 
    for i in range(params.pop_size):
    
        sum = sum + probabilities[i]
        accumulated_reciprocals.append(sum)
    
    for indiv in range(params.pop_size):
        parent_selector = rn.uniform(0.0, accumulated_reciprocals[params.pop_size-1])
        counting_along_for_parent = -1
        
        
        #so here we select a parent. The probability is represented by the size of the gaps
        for i in range(len(accumulated_reciprocals)):
            if parent_selector <= float(accumulated_reciprocals[i]):
                counting_along_for_parent = counting_along_for_parent + 1
        parent_1 = counting_along_for_parent
        
        #temporary - moved by the while statement ... just saves lines 
        parent_2 = counting_along_for_parent
        
        #here we make sure the two parents are not the same
        while parent_2 == parent_1:
            parent_selector = rn.uniform(0.0, accumulated_reciprocals[params.pop_size-1])
            counting_along_for_parent = -1
            for i in range(len(accumulated_reciprocals)):
                if parent_selector <= float(accumulated_reciprocals[i]):
                    counting_along_for_parent = counting_along_for_parent + 1
            parent_2 = counting_along_for_parent

        parent_list.append([parent_1, parent_2])
    return parent_list
                
def mate_and_mutate(parent_list, parent_generation, MeanResults, Measured_shifts, params):
    ''' 
    The aim of this function is to take the list of 'parents couples' and produce the 
    new individual from each couple

    Parameters
    =========
    parent_list : list  
    
    parent_generation : dict 
        the previous generation 

    MeanResults : dict 
        contains the theoretical chemical shifts

    Measured_shifts : dict 
        contains the experimental shifts

    params : Parameters class 

    Returns
    =======
    new_generation : dict

    '''
    new_generation = {}
    individual_counter = 0
    number_of_states = params.number_of_states


    #here we use a single break point which is randomly chosen
    for parent in parent_list:

        current_individual = []
        break_point = rn.randint(0,number_of_states-1)
        for i in range(number_of_states):
            if i <= break_point: 
                current_individual.append(parent_generation[parent[0]][i])
            else: 
                current_individual.append(parent_generation[parent[1]][i])
            
        sum = 0
        current_individual_2 = []
        
        #here is where the mutation occurs
        mutation_deicder = rn.uniform(0,99)
        if mutation_deicder <= params.Mutation_rate:
            
            state_decider = rn.randint(0, number_of_states-1)
            current_individual[state_decider] = rn.uniform(0,1)
        
        #so in order to make this GA converge I had to include this section
        #here if the new individual is less fit than one of the parents we take the 
        #most fit parent
        family = {}
        family[0] = parent_generation[parent[0]]
        family[1] = parent_generation[parent[1]]
        family[2] = current_individual

        family_shifts = cs.Computed_shifts(family,  MeanResults,params)     
        family_distances = distFuncs.Distance(family_shifts, Measured_shifts,params)
        
        Shortest_distance_in_family = 1000
        family_member = 4
        for item in family_distances:
            if family_distances[item] <= Shortest_distance_in_family: 
                Shortest_distance_in_family = family_distances[item]
                family_member = item
                
        #new_generation[individual_counter] = current_individual_2
        new_generation[individual_counter] = family[family_member]
        individual_counter =  individual_counter + 1
        
    return new_generation
    
def Simulation_end_conditions(gen_counter, distances, params):
    '''
    Determines the end of the simulation 
    
    Parameters
    =========
    gen_counter : int
        counts the generation we are on 
    
    distances : dict
        dictionary of the distances 
    
    params : Parameters class 
    
    Returns
    =======
    output : int
    '''
    
    #find the smallest distance
    smallest = 100
    for item in distances:
        if distances[item] <= smallest:
            smallest = distances[item]

    #give a maximum length to the number of generations 
    if gen_counter == params.Max_gen_number:
        output = 1  
        
    #the cut off for when things are 'equal'
    elif smallest <= params.Close_enough_to_measured:
        output = 1
    else: 
        output = 0
    
    #need to add a section which detects when it converges to some value 
    #which is not close enough to the observed values but there is little difference between the
    #generations
    return output


def algo(Mean_CS, Measured_CS, params): 
    '''
    This is the main algorithm 

    Parameters
    ==========
    Mean_CS : dict 
        calculated chemical shifts 
    Measured_CS : dict 
        measured chemical shift 
    params : Parameters class 

    '''

    print 'Max gen number: ' , params.Max_gen_number
    print 'population size: ', params.pop_size

    try:
        os.mkdir( params.write_out_location )
    except OSError:
        print 'directory already exists'
    
    if params.verbose == '1':
        print '''#====================================================================
#starting ...
#===================================================================='''

    #smallest distance log
    if params.writeDistance != None:
        smallest_dist_in_generations = []

    print 'using %s to drive the fitness function' %(params.disMeasure.__name__)
    initial_pop = {}             #starting population
    
    for i in range(params.number_of_simulations):
     
        #the first gen is produced randomly
        Generation_1 = First_Gen(i, params) 
        
        IO.Write_out(1,Generation_1, i, 'Generations', params)
    
        #now we work out the shifts 
        indiv_shifts = cs.Computed_shifts(Generation_1,  Mean_CS,params)
    
        #these are the distances bewteen the generated shifts and the observed shifts 
        

        distances = distFuncs.Distance(indiv_shifts, Measured_CS,params)
        IO.Write_out(1, distances, i, 'Distances', params)
    
        #this is used to make the probabilities 
        recpirical = distFuncs.reciprocal(distances)    
    
        #a list of the 'parent couples'
        parents = select_parents(Generation_1, recpirical, params)
    
        #this line is just for keeping the names the same 
        Parent_Generation = Generation_1.copy()



        generation_Counter = 1
    
        end_of_sim = 0
        while end_of_sim != 1:       
        
            generation_Counter = generation_Counter + 1
            
            if params.verbose == '1':
                div_by_50 =  generation_Counter % 10
                if div_by_50 == 0: 
                    print generation_Counter
        
            #Here we make the new generation
            new_gen = {}
            new_gen = mate_and_mutate(parents, Parent_Generation, Mean_CS, Measured_CS, params) 
        
        
            IO.Write_out(generation_Counter, new_gen, i, 'Generations', params)
        
            #work out the shifts and the distances
        
            indiv_shifts = {}
            indiv_shifts = cs.Computed_shifts(new_gen,  Mean_CS,params)
        
            distances = {}
            distances =  distFuncs.Distance(indiv_shifts, Measured_CS,params)
            IO.Write_out(generation_Counter, distances, i, 'Distances', params)
    
            #do the reciprocal
            recpirical =  distFuncs.reciprocal(distances)    
     
            #a list of the 'parent couples'
            parents = {}
            parents = select_parents(Parent_Generation, recpirical, params)
    
            Parent_Generation = {}
            Parent_Generation = new_gen.copy()
            
            if params.writeDistance != None:
                s_indiv, smallest_dist =  distFuncs.get_smallest_distance(distances)
                smallest_dist_in_generations.append((generation_Counter,smallest_dist))

            end_of_sim = Simulation_end_conditions(generation_Counter, distances,params)
    
        #so here we wite out the best individual in the last generation to solutions
        #I guess we could also look for the best solution in all the generations but this might late longer

        best_individual, dist =  distFuncs.get_smallest_distance(distances)
        
        IO.Wite_solutions(best_individual, Parent_Generation, Mean_CS,Measured_CS,params)
        
        if params.writeDistance != None:
            IO.write_smallest_distances(smallest_dist_in_generations, params)


    #print 'coffee time is over!' 
    #print '==================='
    #print 'It took', time.time()-start, 'seconds.'
    #print '==================='
