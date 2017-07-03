#! /usr/bin/python
# written by Lucas Siemons

import os
import sys
import glob
import re
import random as rn 
import numpy as np 
import time
import pickle as pic

#for a time stamp
start = time.time()

#this meas we will be writing everything out!
#verbose_text = 0

#write out each generation 1 = yes, 0 = now
#write_gen = 0

#chance of zero when making the initial individual state being assigned a zero populaiton
#zero_chance = 5

#some variables
#states used: this is the order of the staets and the  number used in the algorithm  
#usedCC = [[-175,167],[-167,65],[-64,-169],[-54, 91],[-76, 68],[-58,-60],[60,174],[55,84],[-64,169]]
#usedCC = [[-175,167],[-167,65],[-76, 68],[-58,-60],[60,174],[55,84],[-64,169]]

#usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]

#for later use these are the backbone angles
#usedPP = np.array([[-64, -45], [-71, -18],[-112, 8],[-105, -45],[-130, 157],[-126, 129],[-115, 126],[-100, 126],[-72, 132]])

#these are the carbon identifiers
#cList = ['1','2','4','5','6','8','10','12','24']
#Side_chain_carbons = ['4','5', '6','12','24']

#pop_size_global = 20  #size of initial population

#i think the filtering in the cross over stage give a preference to higher mutation rates
#but ... again maybe its just an idea, not tested because when we filter them out a bad mutation will be disregarded
#so we only take good mutations ...?
#reddit seems to suggest that 1-5% is best but 10 also seems good aswel as much higher  
#when i raised to to 80% it took 2/3 of the generations (~200)
#when it was 95 it took 170 gen


#Mutation_rate_global = 60.0
#number_of_simulations_global = 1
simulation_Number = 0
simulation_Path ='Generations/Simulation_'
#number_of_states = len(usedCC)



#these are used to determine the end conditions
#Max_gen_number_global = 20
Close_enough_to_measured = 0.000000001

#====================================================================
# Class to read in the input file
#====================================================================

#this class reads in the input file and bundles the 
#parameters for the GA.

class ReadInputFile:
    def __init__(self, inputFile):
        f = open(inputFile, 'r')
        
        for line in f.readlines():
            
            if line[0] != ';':
                split = line.split()    
                
                #these are the input feilds 
                
                #distance measure
                if split[0] == 'dist_func':
                    if split[1] == 'Euclid_dis':
                        self.disMeasure = Euclid_dis
                    else:
                        print 'The distance measure does noe exist'
                        sys.exit()
                
                #population size
                if split[0] == 'pop_size':
                    self.popSize = int(split[1])
                
                if split[0] == 'Mutation_rate':
                    self.Mutation_rate = int(split[1])
                
                if split[0] == 'verbose':
                    self.verbose = int(split[1])
                
                if split[0] == 'Max_gen_number':
                    self.Max_gen_number = int(split[1])
                
                if split[0] == 'write_out_location':
                    self.write_out_location = split[1]
                
                if split[0] == 'residue':
                    if split[1] == 'ILE':
                        self.usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]
                        self.usedCC = [['t','t'], ['t','p'],['m','m'],['p','t'],['m','t']]
                        self.cList = ['1','2','4','5','6','8','10','12','24']
                        self.Side_chain_carbons = ['4','5', '6','12','24']
                
                if split[0] == 'CalcRes':
                    self.calcRes = split[1]
                    
                if split[0] == 'ExpShifts':
                    self.expShifts = split[1]
                
                if split[0] == 'zero_chance':
                    self.zeros = int(split[1])

                if split[0] == 'write_gen':
                    self.write_gen = int(split[1])
                
                if split[0] == 'sol_name':
                    self.sol_name = split[1]

        f.close()

#====================================================================
#Defining functions
#====================================================================

#This reads in the measured shifts
def ReadIn_measured_C_shifts(input_file):
    #note in the original configuration the inputfile is 'Measured_shifts.txt'
    #these are the target / observed chemical shiftes we want to determine the rotamer populations for
	ReadIn_shifts = {}
	inputshifts = open(input_file, 'r')
	for line in inputshifts.readlines():
	    split = line.split()
	    ReadIn_shifts[split[0]] = split[1]
	return ReadIn_shifts

#this function creates a dictionary of the first generation

def read_in_mean_results(res_file):
    #gonna take my mean results form the results.dat file
    #this should be in the same directory as the current script
    res = open(res_file,'r')
    file_checkpoint = 0
    MeanResults = {}
    for line in res.readlines():
        split = line.split()
        try:
            if file_checkpoint == 1:
                MeanResults[split[0], split[1],split[2]] = split[3]
            if split[0] == 'Mean':
                file_checkpoint = 1
        except IndexError:
            pass 
    res.close()
    return MeanResults

def read_in_sec_struc_res(res_file):
    ''' this function is here to read in the chemical shifts according to the secondry structure
    where as before we had only a random coil '''

    Results = {}
    res = open(res_file, 'r')
    
    #here we read in the results from the .dat file
    check = 0
    for line in res.readlines():
        split = line.split()
        
        try:
            if split[0] == '#Mean':
                check = 1
            #this check is so we ony read the first part of the data file
            if check != 1:
                if line[0] != '#':
                    Results[split[0], split[1], split[2],split[3],split[4]] = float(split[5])
        except IndexError:
            pass
    bSheetData = {}
    aHelixData = {}
    #mean dictionaries for the secondry structures
    MeanBSheetData = {}
    MeanAHelixData = {}
    
    for key in Results:

        # these just define the redion of the aHelix in the Ramachadran plot
        if float(key[1]) >= float(-140) and float(key[1]) <= float(-30):
            if float(key[2]) >= float(-91) and float(key[2]) <= float(38):
                try:

                    aHelixData[key[0],key[3],key[4]].append(Results[key])
                except KeyError:
                    aHelixData[key[0],key[3],key[4]] = []
                    aHelixData[key[0],key[3],key[4]].append(Results[key])
                #print 'added something to a list'
            
        # these just define the redion of the bsheet in the Ramachadran plot
        if float(key[1]) >= float(-170) and float(key[1]) <= float(-42):
            if float(key[2]) >= float(65) and float(key[2]) <= float(180):
                try:
                    bSheetData[key[0],key[3],key[4]].append(Results[key])
                except KeyError:
                    bSheetData[key[0],key[3],key[4]] = []
                    bSheetData[key[0],key[3],key[4]].append(Results[key])
                #print 'added something to b list'

    #making the means results for 

    for key in bSheetData:
        meanBCs = np.divide(sum(bSheetData[key]),len(bSheetData[key]))
        MeanBSheetData[key] = meanBCs
        
        #these lines just gives the standard deviation should you need it 
        #meanStd = np.std(bSheetData[key])
        #MeanBSheetData[key].append(meanStd)

    for key in bSheetData:
        meanACs = np.divide(sum(aHelixData[key]),len(aHelixData[key]))
        MeanAHelixData[key] = meanACs    
    res.close()
    

    return MeanAHelixData, MeanBSheetData


def First_Gen(Sim_number_iteration, pop_size):
    
    ''' this function aims to produce the first generation, these are produced by the random number generator
    rn.random()'''
    
    WorkingRandomList = []
    Random_Sum = 0.0
    RandomFloat = 0.0
    NormalisedList = []
    
    individuals = {}
    
    for indiv in range(pop_size):

        NormalisedList = []
        WorkingRandomList = []
        
        Random_Sum = 0
        for i in range(number_of_states):
            #here the first generation is made randomly 
            chance = rn.randint(0,100)
            if chance <= zero_chance:
        
                RandomFloat = 0.0000000001
            else:
                RandomFloat = rn.random()
            WorkingRandomList.append(RandomFloat )
            Random_Sum = Random_Sum + RandomFloat 
        
        #note here we normalise the vector so all the elements sum to 1 but this is the only 
        #time it is done. Later the random numbers are taken and they are normalised in
        #the fitness function when we compute the shifts in Computed_shifts
        for i in range(number_of_states):
            normalised_state_pop =  np.divide(WorkingRandomList[i],Random_Sum)
            NormalisedList.append(normalised_state_pop)
        
        individuals[indiv] = NormalisedList
    return individuals

def Computed_shifts(individuals, rotamer_shifts): #not right yet

    '''
    The aim of this function is to take all the individuals, ie the populations 
    of each state and work out the chemical shifts for the 5 carbon atoms on the ILE 
    side chain. This is stored in the Dictionary Shifts, where the key is the individual identifier 
    and the entry is the list of the chemical shifts. These are in the same order as the 
    list of carbons in Side_chain_carbons
    
    individuals is a dict where the key identifies the individual and the entry is the list of values which are
    used to give the populations
    
    rotamer_shifts is a dict containing the chemical shifts from the DFT, currently it uses mean results 
    which were printed out at the end of my summer project'''

    Shifts = {} #shifts given by each individual, the are calculated using the 'MeanResults'
    
    for key in individuals:
        
        #just in case mostly appeared in an attempt to debug 
        side_chain_shifts = [0,0,0,0,0]
        side_chain_shifts_counter = 0
        totalShift = 0
        
        #so here i now need to normalsie the input vector so the sum of the elements is 1      
        current_individual = []
        current_individual = individuals[key]
        element_sum = 0.0
        for i in range(len(current_individual)):
            element_sum = element_sum + current_individual[i]
        
        #here we normalise the individual because in generations after the first corss over 
        #does not require the elements to sum to 1
        current_individual_norm = []
        for i in range(len(current_individual)):
            working_element = np.divide(current_individual[i], element_sum)
            current_individual_norm.append(working_element)
                    
        for sidechain_Carbon in Side_chain_carbons:
            totalShift = 0
            
            for i in range(len(usedCC)):
                current_state_shift = 0.0
                #print rotamer_shifts[sidechain_Carbon, str(usedCC[i][0]), str(usedCC[i][1])]
                
                current_state_shift =  float(current_individual_norm[i]) * float(rotamer_shifts[sidechain_Carbon, str(usedCC[i][0]), str(usedCC[i][1])])
                totalShift = totalShift + current_state_shift


            side_chain_shifts[side_chain_shifts_counter] = totalShift
            side_chain_shifts_counter = side_chain_shifts_counter + 1
        

        
        Shifts[key] = side_chain_shifts

    return Shifts
      
def Distance(computed, observed,dist_func):
    '''here we produce a set of Euclidian distances between the set of shifts for each population member and 
    the observed chemical shift this is written to a dictionary
    
    this forms the basis of the PDF function used to select parents
    
    computed is a dictionary with the computed shifts, the key is the individual it refers to and the 
    entry is the list of chemical shifts in the same order as sidechain carbons
    
    observed is the observed / measured chemical shifts
    
    NOTE - should check this part... when i get a real example

    dist_func - this allows you to choose the distance function that compares the 
                individual's CS with the observed one


    '''
    
    distance = {}

    for individual in computed: 
        args=(individual,observed,computed)
        distance[individual] = dist_func(*args)
    
    return distance

def Euclid_dis(individual,observed,computed):
    distance_sum = 0
    for carbon in range(len(Side_chain_carbons)):
        if observed[Side_chain_carbons[carbon]] != "none":
                
            #this bit is just the sum for the distance over many lines
            a = float(computed[individual][carbon])
            b =  float(observed[Side_chain_carbons[carbon]])
            c = (a - b)
            d = c**2
            distance_sum = distance_sum+ d

    final_distance = np.sqrt(distance_sum)
    return final_distance

def Chi_squared(individual,observed,computed):
    distance_sum = 0
    for carbon in range(len(Side_chain_carbons)):
        if observed[Side_chain_carbons[carbon]] != "none":
                
            #this bit is just the sum for the distance over many lines
            a = float(computed[individual][carbon])
            b =  float(observed[Side_chain_carbons[carbon]])
            c = (a - b)
            d = c**2
            distance_sum = distance_sum+ d

    final_distance = distance_sum
    return final_distance

def compare_populations(generated_rotamers):
    '''This also just does a euclidian distance between the generated population and a 
    solution but since there is a slightly different functionality compared to Euclidean_Distance()
    i decided to be lazy ...'''
    
    distance = 0 
    
    #try statement is there for the first simulation which would have no previous solutions to compare to
    try:
        
        solutions = open('Solutions.txt', 'r')
        #each line in S_pops is a set of rotamers which is a solution 
        S_pops = solutions.readlines()
        smallest_distance = 100

        for list in S_pops:
            distance_sum = 0
            list = list.split()
            
            for item in range(len(list)):

                a = float(generated_rotamers[item])
                b = float(list[item])
                c = (a - b)
                d = c**2
                distance_sum = distance_sum+ d        
            
            if distance_sum <= smallest_distance:
                smallest_distance = distance_sum
        
        distance = smallest_distance
        
    except IOError: 
        distance = 1

    return distance
    
def reciprocal(dict, previous_gen):
    
    '''I do this so that later it can be used to make large values less likely to be selected
    than small distances, ie selecting for vectors closer to the terget "observed value" 
    
    here the key of the dict also refers to the individual
    
    also here i need to add the term which adds a bit to compare it to previous solutions'''
    
    
    Recip = {}
    #print dict
    for key in dict: 
        #a = np.divide(1,np.exp((   dict[key])))

        a = np.divide(1,(dict[key]))
        
        #the next two lines allow us to select for values which are not close 
        #to previously found solutions ... need to check but could give a better
        #sampling of the solution space 

        b = compare_populations(previous_gen[key])
        c = a*b
        Recip[key] = c
    return Recip

def select_parents(daughter_generation, probabilities, pop_size):
    
    '''This function aims to produce a list of 'parent couples' which are 
    then crossed over to make each individual
    
    daughter gen is a dictionary like individuals
    
    the probabilities are the reciprocals ... need to write a bit up on the cross over 
    will probably do this in a Tex document '''
    
    sum = 0
    accumulated_reciprocals = []
    parent_list = []
    
    #this is used to do the weighted PDF 
    for i in range(pop_size):
    
        sum = sum + probabilities[i]
        accumulated_reciprocals.append(sum)
    
    for indiv in range(pop_size):
        parent_selector = rn.uniform(0.0, accumulated_reciprocals[pop_size-1])
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
            parent_selector = rn.uniform(0.0, accumulated_reciprocals[pop_size-1])
            counting_along_for_parent = -1
            for i in range(len(accumulated_reciprocals)):
                if parent_selector <= float(accumulated_reciprocals[i]):
                    counting_along_for_parent = counting_along_for_parent + 1
            parent_2 = counting_along_for_parent

        parent_list.append([parent_1, parent_2])
    return parent_list
                
def mate_and_mutate(parent_list, daughter_generation, MeanResults, Measured_shifts, Mutation_rate,func):
    
    ''' the aim of this function is to take the list of 'parents couples' and produce the 
    new individual from each couple'''
    
    #NOTE 'daughter_generation' really ought to be called parent gen

    new_generation = {}
    individual_counter = 0
    
    #here we use a single break point which is randomly chosen
    for parent in parent_list:

        current_individual = []
        break_point = rn.randint(0,number_of_states-1)
        for i in range(number_of_states):
            if i <= break_point: 
                current_individual.append(daughter_generation[parent[0]][i])
            else: 
                current_individual.append(daughter_generation[parent[1]][i])
            
        sum = 0
        current_individual_2 = []
        
        #here is where the mutation occurs
        mutation_deicder = rn.uniform(0,99)
        if mutation_deicder <= Mutation_rate:
            
            state_decider = rn.randint(0, number_of_states-1)
            current_individual[state_decider] = rn.uniform(0,1)
        
        #so in order to make this GA converge I had to include this section
        #here if the new individual is less fit than one of the parents we take the 
        #most fit parent
        family = {}
        family[0] = daughter_generation[parent[0]]
        family[1] = daughter_generation[parent[1]]
        family[2] = current_individual

        family_shifts = Computed_shifts(family,  MeanResults)     
        family_distances = Distance(family_shifts, Measured_shifts,func)
        
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

def Write_out(gen_number, populaton_dict, sim_number, dict_type, write_out, pop_size, out_dir):
    
    if write_out == 0:
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
    file = open(dir_name,'w')
    
    for i in range(pop_size):
        file.write(str(i)+'  ')
        try: 
            for index in range(len(populaton_dict[i])):
                file.write(str(populaton_dict[i][index])+ '  ')
            file.write('\n')
        except TypeError:
            #this just makes the function more flexible, as some dicts only have single vlaues in them
            file.write(str(populaton_dict[i])+ '\n')
    file.close() 
        
def Simulation_end_conditions(gen_counter, max_generations, Close_enough_to_measured, distances):
    

    
    #find the smallest distance
    smallest = 100
    for item in distances:
        if distances[item] <= smallest:
            smallest = distances[item]

    #give a maximum length to the number of generations 
    if gen_counter == max_generations:
        output = 1  
        
    #the cut off for when things are 'equal'
    elif smallest <= Close_enough_to_measured:
        output = 1
    else: 
        output = 0
    
    #need to add a section which detects when it converges to some value 
    #which is not close enough to the observed values but there is little difference between the
    #generations
    return output

def Wite_solutions(best_indiv_id, final_gen, out_dir, sol_name, Mean_CS,func,Measured_CS): 
    
    #here we write out the populations
    out_name = out_dir+sol_name
    print 'outname:  ', out_name
    S_file = open(out_name, 'a')
    
    total = sum(final_gen[best_indiv_id])

    for i in range(len(final_gen[best_indiv_id])):

        S_file.write( str(np.divide(final_gen[best_indiv_id][i], total) ) + '   ' )
    S_file.write('\n')
    S_file.close()
    
    #here we write out the final chemical shifts we obtainand the final distance between these and the 
    #input chemical shift
    
    #the file for the best individual to be written to
    best_indiv_name = out_name + '_bestIndiv.txt'
    
    #this is just some reformatting so we can use Computed_shifts()
    dummyDict = {}
    dummyDict['dummy'] = final_gen[best_indiv_id]

    #calculate the chemical shift
    closest_shifts = Computed_shifts(dummyDict, Mean_CS)
    
    print 'best f name: ', best_indiv_name
    print 'cs' , closest_shifts
    bestIndiv = open(best_indiv_name, 'w')
    
    #write out the chemical shift
    print closest_shifts
    for atom, cs in zip(Side_chain_carbons, closest_shifts['dummy']):
        bestIndiv.write('%s   %s\n'% (atom, cs))


    #calculate the distance
    distance = Distance(closest_shifts, Measured_CS,func)
    
    bestIndiv.write('Distance: %f' % (distance['dummy']))
    bestIndiv.close()


   

#====================================================================
#Doing stuff
#====================================================================

#setting up the folder structure here ------ 
#Writing out the first population ----------------------

def setup_Fitted_DFT(Mean_CS_file,  Measured_CS_file, sec_struct,params, read_measured_cs='yes'):

    '''
    This setup file used the chemical shifts derived from the fitted DFT minima 
    these are currently held in CS_ile_DFT.dat

    the secondary structure can be A B or R for alpha beta or random coil 

    '''

    alphaDFT = {}
    betaDFT  = {}
    rcDFT    = {}

    usedCC = params.usedCC
    dft = open(Mean_CS_file, 'r')
    #print 'states used:   ',  usedCC
    
    for line in dft.readlines():
        split = line.split()
        #print line
        if split[0] == 'alpha':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                alphaDFT['4', angles[0], angles[1]] = float(split[2])
                alphaDFT['5', angles[0], angles[1]] = float(split[3])
                alphaDFT['6', angles[0], angles[1]] = float(split[4])
                alphaDFT['12', angles[0], angles[1]] = float(split[5])
                alphaDFT['24', angles[0], angles[1]] = float(split[6])

        if split[0] == 'beta':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                betaDFT['4', angles[0], angles[1]] = float(split[2])
                betaDFT['5', angles[0], angles[1]] = float(split[3])
                betaDFT['6', angles[0], angles[1]] = float(split[4])
                betaDFT['12', angles[0], angles[1]] = float(split[5])
                betaDFT['24', angles[0], angles[1]] = float(split[6])

        if split[0] == 'rCoil':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                rcDFT['4', angles[0], angles[1]] = float(split[2])
                rcDFT['5', angles[0], angles[1]] = float(split[3])
                rcDFT['6', angles[0], angles[1]] = float(split[4])
                rcDFT['12', angles[0], angles[1]] = float(split[5])
                rcDFT['24', angles[0], angles[1]] = float(split[6])

    dft.close()

    if sec_struct == 'R':
        selected_Mean_CS = rcDFT
    elif sec_struct == 'A':
        selected_Mean_CS = alphaDFT
    elif sec_struct == 'B':
        selected_Mean_CS = betaDFT
    else:
        print '''ERROR: one of the options selected was not a recognised secondry structure type. 
        Pease check the set_up() function'''
        sys.exit()
    #read in the exprimental results
    if read_measured_cs == 'yes':
        Measured_CS_1 = ReadIn_measured_C_shifts(Measured_CS_file)
    else:
        print 'attention not readin measured CS'
        Measured_CS_1 = ''
    dft.close()

    #print 'alpha ------------------'
    #print alphaDFT
    #print '------------------'

    #print 'beta ------------------'
    #print betaDFT
    #print '------------------'

    #print rcDFT
    #print '------------------'


    return selected_Mean_CS, Measured_CS_1

def set_up(Mean_CS_file,  Measured_CS_file, sec_struct):
    ''' Note that sec_struct can be set to A, B or R for alpha helix 
        beta sheet of random coil / unknown '''
    #this part checks that secondary  structure. -- could be done better because it reads Results.dat twice
    R_coil_Mean_CS = read_in_mean_results(Mean_CS_file)
    a_helix_CS, b_sheet_CS = read_in_sec_struc_res(Mean_CS_file)
    Measured_CS_1 = ReadIn_measured_C_shifts(Measured_CS_file)
    
    if sec_struct == 'R':
        selected_Mean_CS = R_coil_Mean_CS
    elif sec_struct == 'A':
        selected_Mean_CS = a_helix_CS
    elif sec_struct == 'B':
        selected_Mean_CS = b_sheet_CS
    else:
        print '''ERROR: one of the options selected was not a recognised secondry structure type. 
        Pease check the set_up() function'''
        sys.exit()
    return selected_Mean_CS, Measured_CS_1

def algo(Mean_CS, Measured_CS, params): 
    
    #setting up input params
    global verbose_text
    verbose_text = params.verbose
   
    global verbose
    verbose = params.verbose
    
    global write_gen 
    write_gen = params.write_gen
    
    global zero_chance
    zero_chance = params.zeros
    
    global usedCC
    usedCC = params.usedCC

    global cList
    cList = params.cList
    
    global Side_chain_carbons
    Side_chain_carbons = params.Side_chain_carbons

    global pop_size
    pop_size = params.popSize
    
    global Mutation_rate
    Mutation_rate = params.Mutation_rate
    
    global number_of_states
    number_of_states = len(usedCC)
    
    global out_direc
    out_direc = params.write_out_location
    
    global func
    func = params.disMeasure
    
    global number_of_simulations
    number_of_simulations = 1
    
    global Max_gen_number
    Max_gen_number = params.Max_gen_number
    
    global sol_name
    sol_name = params.sol_name

    print 'Max gen number: ' , Max_gen_number
    print 'population size: ', pop_size

    try:
        os.mkdir( out_direc )
    except OSError:
        print 'directory already exists'
    
    if verbose == '1':
        print '''#====================================================================
#starting ...
#===================================================================='''

    
    print 'using %s to drive the fitness function' %(func.__name__)
    initial_pop = {}             #starting population
    
    for i in range(number_of_simulations):
     
        #the first gen is produced randomly
        Generation_1 = First_Gen(i, pop_size) 
        
        Write_out(1,Generation_1, i, 'Generations', write_gen, pop_size, out_direc)
    
        #now we work out the shifts 
        indiv_shifts = Computed_shifts(Generation_1,  Mean_CS)
    
        #these are the distances bewteen the generated shifts and the observed shifts 
        distances = Distance(indiv_shifts, Measured_CS,func)
        Write_out(1, distances, i, 'Distances', write_gen, pop_size, out_direc)
    
        #this is used to make the probabilities 
        recpirical = reciprocal(distances, Generation_1)    
    
        #a list of the 'parent couples'
        parents = select_parents(Generation_1, recpirical, pop_size)
    
        #this line is just for keeping the names the same 
        Parent_Generation = Generation_1.copy()

        generation_Counter = 1
    
        end_of_sim = 0
        while end_of_sim != 1:       
        
            generation_Counter = generation_Counter + 1
            
            if verbose == '1':
                div_by_50 =  generation_Counter % 10
                if div_by_50 == 0: 
                    print generation_Counter
        
            #Here we make the new generation
            new_gen = {}
            new_gen = mate_and_mutate(parents, Parent_Generation, Mean_CS, Measured_CS, Mutation_rate,func) 
        
        
            Write_out(generation_Counter, new_gen, i, 'Generations', write_gen, pop_size, out_direc)
        
            #work out the shifts and the distances
        
            indiv_shifts = {}
            indiv_shifts = Computed_shifts(new_gen,  Mean_CS)
        
            distances = {}
            distances = Distance(indiv_shifts, Measured_CS,func)
            Write_out(generation_Counter, distances, i, 'Distances', write_gen, pop_size, out_direc)
    
            #do the reciprocal
            recpirical = reciprocal(distances,new_gen)    
     
            #a list of the 'parent couples'
            parents = {}
            parents = select_parents(Parent_Generation, recpirical, pop_size)
    
            Parent_Generation = {}
            Parent_Generation = new_gen.copy()
        
            end_of_sim = Simulation_end_conditions(generation_Counter, Max_gen_number, Close_enough_to_measured, distances)
    
        #so here we wite out the best individual in the last generation to solutions
        #I guess we could also look for the best solution in all the generations but this might late longer

        shortest_distance = 100
        best_individual = -1
        for item in distances: 
            if distances[item] <= shortest_distance:
                shortest_distance = distances[item]
                best_individual = item
        
        Wite_solutions(best_individual, Parent_Generation, out_direc, sol_name, Mean_CS,func,Measured_CS)
    
    #print 'coffee time is over!' 
    #print '==================='
    #print 'It took', time.time()-start, 'seconds.'
    #print '==================='


if __name__ == '__main__':
    
    params = ReadInputFile('InputFile.txt')
    calRes = params.calcRes
    expshifts = params.expShifts
    
    dist_func = Euclid_dis
    #print 'using '+ str(dist_func.__name__)
    #average_CS, experimental_CS = set_up(calRes,  expshifts,'R')
    #algo(Mean_CS,   Measured_CS,     write_gen, pop_size,        Mutation_rate,        number_of_simulations,        verbose,      Max_gen_number, )
    average_CS, experimental_CS = setup_Fitted_DFT(calRes,  expshifts,'R', params)
    algo(average_CS, experimental_CS, params)





    
