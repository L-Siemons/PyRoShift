
'''
This module is for functions which act to manipulate or calculate chemical 
shifts.
'''

import numpy as np

def Computed_shifts(individuals, rotamer_shifts,params): 

    '''
    
    This function converts the populations to chemical shifts

    Parameters
    ==========
    individuals : dictionary
        The key is the identifier for each member of the population
        The entry is the pseudo population distribution
    rotamer_shifts : dictionary
        Contains the chemical shifts used by the algorithm 
    
    Returns
    =======
    Shifts : dictionary
        key - individual id 
        entry - list of chemical shifts
    '''

    Shifts = {} #shifts given by each individual
    
    for key in individuals:
        
        #just in case mostly appeared in an attempt to debug         
        side_chain_shifts = [0.] * params.number_of_states
        side_chain_shifts_counter = 0
        
        #so here i now need to normalsie the input vector so the sum of the elements is 1      
        element_sum = sum(individuals[key])
        current_individual_norm = []       

        for i in range(params.number_of_states):
            working_element = np.divide(individuals[key][i], element_sum)
            current_individual_norm.append(working_element)
        
        #calculate the chemical shifts as a weighted average
        for sidechain_Carbon in params.Side_chain_carbons:
            totalShift = 0
            
            for i in range(params.number_of_states):
                current_state_shift = 0.0
                pop = float(current_individual_norm[i])
                state_shift = float(rotamer_shifts[sidechain_Carbon, params.usedCC[i][0], params.usedCC[i][1]])
                current_state_shift = pop * state_shift
                totalShift = totalShift + current_state_shift


            side_chain_shifts[side_chain_shifts_counter] = totalShift
            side_chain_shifts_counter = side_chain_shifts_counter + 1

        Shifts[key] = side_chain_shifts

    return Shifts