'''
This module contains functions that read and write files
'''

import os
import sys

import numpy as np
from pyroshiftGA import chemicalShifts as cs
import pyroshiftGA.distanceFunctions as distFuncs


#===============================================================
# Read in file
#===============================================================

def ReadIn_measured_C_shifts(input_file):
    '''
    Reads in the measured/observed chemical shifts

    Parameters
    =========
    input_file : str
        Name of the file containing the shifts

    Returns
    =======
    ReadIn_shifts : dict
        Key - atom id
        entry is the chemical shift

    '''

    ReadIn_shifts = {}
    inputshifts = open(input_file, 'r')
    for line in inputshifts.readlines():
        split = line.split()
        ReadIn_shifts[split[0]] = split[1]
        if split[0] == 'SSE':
            SSE = split[1]

    return ReadIn_shifts, SSE

def setup_Fitted_DFT(dft_CS_file, sec_struct, params, read_measured_cs='yes'):

    '''
    This function reads in the appropriate calculated chemical shifts
    Given the secondary structure.

    It does this be reading all the theoretical chemical shifts and then
    selecting the right set based on the secondary structure given in the
    file with the measured chemical shifts.

    These can be either R, A or B for random coil, alpha helix or beta sheet.

    Parameters
    =========
    dft_CS_file : str
        file with the theoretical chemical shifts
    sec_struct :
        secondry structure
    params : Parameters class
    read_measured_cs : str
        select if to read in the measured chemical shift or not

    Return
    ======
    selected_Mean_CS : dict
        contains the selected theoretical chemical shifts
    '''

    alphaDFT = {}
    betaDFT = {}
    rcDFT = {}

    usedCC = params.usedCC
    dft = open(dft_CS_file, 'r')
    #print 'states used:   ',  usedCC

    for line in dft.readlines():
        split = line.split()
        angles = list(split[1])
        for indx, atom in enumerate(params.Side_chain_carbons):
            keys = tuple([atom] + angles)

            if split[0] == 'alpha':

                if angles in usedCC:
                    alphaDFT[keys] = float(split[2+indx])

            if split[0] == 'beta':
                if angles in usedCC:
                    betaDFT[keys] = float(split[2+indx])

            if split[0] == 'rCoil':
                if angles in usedCC:
                    rcDFT[keys] = float(split[2+indx])


    dft.close()

    if sec_struct == 'R':
        selected_Mean_CS = rcDFT
    elif sec_struct == 'A':
        selected_Mean_CS = alphaDFT
    elif sec_struct == 'B':
        selected_Mean_CS = betaDFT
    else:
        print '''ERROR: one of the options selected was not a recognized secondary structure type.
        Pease check the set_up() function'''
        sys.exit()

    return selected_Mean_CS


#===============================================================
# write out file
#===============================================================

def Write_out(gen_number, populaton_dict, sim_number, dict_type, params):

    '''
    This function writes out the populations of each generation
    Parameters
    =========
    gen_number : int
        generation
    populaton_dict : dict
        key - individual id
        value - the rotamer distribution
    sim_number : int
        the simulation we are on, so far we only tested running a single one
    dict_type : str
        label whether generations or distances are written out
    params : Parameters class
    '''

    out_direc = params.write_out_location

    if params.write_gen == 0:
        return None

    try:
        os.mkdir(out_dir+str(dict_type))
    except OSError:
        pass

    try:
        os.mkdir(out_dir+str(dict_type)+'/Simulation_'+str(sim_number)+'/')
    except OSError:
        pass

    dir_name = out_dir+str(dict_type)+'/Simulation_'+str(sim_number)+'/gen_'+str(gen_number)+'.txt'
    file = open(dir_name, 'w')

    for i in range(params.pop_size):
        file.write(str(i)+'  ')
        try:
            for index in range(len(populaton_dict[i])):
                file.write(str(populaton_dict[i][index])+ '  ')
            file.write('\n')
        except TypeError:
            #this just makes the function more flexible,
            #as some dicts only have single values in them
            file.write(str(populaton_dict[i])+ '\n')
    file.close()

def Wite_solutions(best_indiv_id, final_gen, Mean_CS, Measured_CS, params):
    '''
    This function writes out the solution from the GA

    Parameters
    =========
    best_indiv_id : int
        the id of the best individual
    final_gen : dict
        the final generation of the GA
    Mean_CS : dict
        the chemical shifts of the states
    Measured_CS : dict
        the measured chemical shifts
    params : Parameters class
    '''

    #here we write out the populations
    out_name = params.write_out_location+params.sol_name
    print 'outname:  ', out_name
    S_file = open(out_name, 'a')

    total = sum(final_gen[best_indiv_id])

    for i in range(len(final_gen[best_indiv_id])):

        S_file.write(str(np.divide(final_gen[best_indiv_id][i], total)) + '   ')
    S_file.write('\n')
    S_file.close()

    #the file for the best individual to be written to
    best_indiv_name = out_name.split('.')[0] + '_bestIndiv.txt'

    #this is just some reformatting so we can use Computed_shifts()
    dummyDict = {}
    dummyDict['dummy'] = final_gen[best_indiv_id]

    #calculate the chemical shift
    closest_shifts = cs.Computed_shifts(dummyDict, Mean_CS, params)

    print 'File with best individual: ', best_indiv_name
    print 'Calculated populations:'
    for i, j in zip(params.usedCC, dummyDict['dummy']):
        name = ''.join(i)
        print '%4s: %0.3f' % (name, j)

    print 'Best chemical shifts:'
    for i in closest_shifts['dummy']:
        print i

    bestIndiv = open(best_indiv_name, 'w')

    #write out the chemical shift
    for atom, chemicalShift in zip(params.Side_chain_carbons, closest_shifts['dummy']):
        bestIndiv.write('%s   %s\n'% (atom, chemicalShift))


    #calculate the distance
    distance = distFuncs.Distance(closest_shifts, Measured_CS, params)

    bestIndiv.write('Distance: %f' % (distance['dummy']))
    bestIndiv.close()

def write_smallest_distances(smallest_dist_in_generations, params):
    '''
    This function writes the smallest distances to a file specified in
    params.writeDisrance

    Parameters
    =========
    smallest_dist_in_generations : list
        [generation , distance]
    params : Parameters class
    '''

    f = open(params.writeDistance, 'w')
    for i in smallest_dist_in_generations:
        f.write('%i  %0.5f\n'%(i[0], i[1]))
    f.close()
