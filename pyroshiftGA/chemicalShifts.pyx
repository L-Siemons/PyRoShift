
'''
This module is for functions which act to manipulate or calculate chemical
shifts.
'''

def Computed_shifts(individuals, rotamer_shifts, params):

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


    cdef int side_chain_shifts_counter
    cdef float element_sum
    cdef int i
    cdef float working_element
    cdef str sidechain_Carbon
    cdef float current_state_shift
    cdef float state_shift
    cdef int key

    Shifts = [] #shifts given by each individual

    for key, entry in enumerate(individuals):

        #just in case mostly appeared in an attempt to debug
        side_chain_shifts = [0.] * params.number_of_states
        side_chain_shifts_counter = 0

        #so here i now need to normalsie the input vector so the sum of the elements is 1
        element_sum = sum(entry)
        current_individual_norm = []

        for i in range(params.number_of_states):
            working_element = entry[i]/element_sum
            current_individual_norm.append(working_element)

        #calculate the chemical shifts as a weighted average
        for sidechain_Carbon in params.Side_chain_carbons:
            totalShift = 0

            for i in range(params.number_of_states):
                current_state_shift = 0.0
                pop = float(current_individual_norm[i])

                Rotmaerkey = tuple([sidechain_Carbon] + params.usedCC[i])
                state_shift = float(rotamer_shifts[Rotmaerkey])

                current_state_shift = pop * state_shift
                totalShift = totalShift + current_state_shift

            side_chain_shifts[side_chain_shifts_counter] = totalShift
            side_chain_shifts_counter = side_chain_shifts_counter + 1

        Shifts.append(side_chain_shifts)
    return Shifts
