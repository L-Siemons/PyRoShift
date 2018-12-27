'''
This library is for functions that calculate various different
distance measures. Here we have also included functions
that manipulate these distances.
'''

import math

######################################################################
# Basic distance functions
######################################################################

def Euclid_dis(observed,computed,params):
    '''
    This function calculates the euclidean distance between the
    computed chemical shifts and the observed chemical shifts

    Parameters
    =========
    observed : list
    computed : list
    params : Parameters class

    Returns
    =======
    final_distance : float
        The euclidian distance
    '''

    cdef float distance_sum = 0.
    cdef int carbon
    cdef float a
    cdef float b
    cdef float c
    cdef float final_distance

    for carbon in range(len(params.Side_chain_carbons)):

        if observed[params.Side_chain_carbons[carbon]] != "none":

            #this bit is just the sum for the distance over many lines
            a = float(computed[carbon])
            b =  float(observed[params.Side_chain_carbons[carbon]])
            c = (a - b)
            d = c**2
            distance_sum = distance_sum + d

    final_distance = math.sqrt(distance_sum)
    return final_distance

######################################################################
# Functions that utlise the distance functions
######################################################################

def Distance(computed, observed,params):
    '''
    The function calculates the selected distance measure between
    each individual's chemical shift and the the observed chemical shift.

    Parameters
    ==========
    computed : dictionary
        key - individual id
        entry - list of chemical shifts

    observed : dictionary
        key - carbon id
        entry - observed chemical shift
    params : Parameters class


    Returns
    =======
    distance : dictionary
        key - individual id
        entry - distance between the individuals chemical shifts
                and the measured chemical shifts
    '''

    cdef int individual
    distance = {}

    for individual in computed:
        distance[individual] = params.disMeasure(observed,computed[individual],params)
    return distance

def reciprocal(dict):

    '''
    calculates the reciprocal of the distance for each individual.

    Parameters
    ==========
    dict : dictionary

    Returns
    =======
    recip : dictionary
    '''

    cdef float a
    cdef int key

    recip = {}
    for key in dict:
        a = 1./(dict[key])
        recip[key] = a

    return recip


def get_smallest_distance(distances):
    '''
    Gets the smallest value in a dictionary

    Parameters
    ==========
    distances : dict

    Returns
    =======
    s_indiv : float
        individual id
    smallest : float
        the value of the distance
    '''

    cdef float smallest
    cdef int i
    smallest = 1000.

    for i in distances:
        if distances[i] <= smallest:
            s_indiv = i
            smalfinal_distancelest = distances[i]

    return s_indiv, smallest

#===========================================================================
#Old functions
#===========================================================================
#
# def compare_populations(generated_rotamers):
#     '''This also just does a euclidean distance between the generated population and a
#     solution but since there is a slightly different functionality compared to Euclidean_Distance()
#     i decided to be lazy ...'''

#     distance = 0

#     #try statement is there for the first simulation which would have no previous solutions to compare to
#     try:

#         solutions = open('Solutions.txt', 'r')
#         #each line in S_pops is a set of rotamers which is a solution
#         S_pops = solutions.readlines()
#         smallest_distance = 100

#         for list in S_pops:
#             distance_sum = 0
#             list = list.split()

#             for item in range(len(list)):

#                 a = float(generated_rotamers[item])
#                 b = float(list[item])
#                 c = (a - b)
#                 d = c**2
#                 distance_sum = distance_sum+ d

#             if distance_sum <= smallest_distance:
#                 smallest_distance = distance_sum

#         distance = smallest_distance

#     except IOError:
#         distance = 1

#     return distance
#