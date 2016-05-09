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
sys.path.append("..") 
import Master_2 as M


#for a time stamp
start = time.time()

verbose = 0

#write out each generation 1 = yes, 0 = now
write_gen = 0

#some variables

#states used: this is the order of the staets and the  number used in the algorithm  
#usedCC = [[-175,167],[-167,65],[-64,-169],[-54, 91],[-76, 68],[-58,-60],[60,174],[55,84],[-64,169]]
#usedCC = [[-175,167],[-167,65],[-76, 68],[-58,-60],[60,174],[55,84],[-64,169]]

usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]

#for later use these are the backbone angles
usedPP = np.array([[-64, -45], [-71, -18],[-112, 8],[-105, -45],[-130, 157],[-126, 129],[-115, 126],[-100, 126],[-72, 132]])

#these are the carbon identifiers
cList = ['1','2','4','5','6','8','10','12','24']
Side_chain_carbons = ['4','5', '6','12','24']

pop_size_global = 250  #size of initial population

Mutation_rate_global = 60.0
number_of_simulations_global = 1
simulation_Number = 0
simulation_Path ='Generations/Simulation_'
number_of_states = len(usedCC)



#these are used to determine the end conditions
Max_gen_number_global = 100
Close_enough_to_measured = 0.000001


grid_search_cuttoff = 2
data_point_number = 1000
os.system('mkdir grid_search_results')
os.system('mkdir synth_cs')
os.system('mkdir algo_solutions_all')
os.system('mkdir algo_solutions_cg2_cd')

synth_pop = {}
synth_pop_CS = {}
grid_search_res = {}

#functions

def get_synth_CS(rotamer_pops, sideChain_atoms, mean_CS, id):

    cs_input = open('synth_shifts_'+str(id)+'.txt', 'w')
    for atom in sideChain_atoms:
        CS_sum = 0.0
        for i2 in range(len(rotamer_pops)):
            CS_sum = CS_sum + rotamer_pops[i2] * float(mean_CS[atom, str(usedCC[i2][0]), str(usedCC[i2][1])])
        
        #now CS sum is the total chemical shift for that atom.
        
        cs_input.write(str(atom)+'  '+str(CS_sum)+'\n')
    cs_input.close()

  
for i, element_1 in enumerate(usedCC):
    print i 
    empty = []
    for indx, element in enumerate(usedCC):    
        empty.append(0.0)
    
    empty[i] = 1.0
    print empty
    indiv = empty
    synth_pop[i] = indiv
    
    input_name = 'synth_shifts_'+str(i)+'.txt'
    
    
    #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
    #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
    #do the set up here
    average_CS, experimental_CS = M.set_up('../Results.dat',  '../Measured_shifts.txt', 'R')

    #make the input file
    print 'run: ' + str(i)
    get_synth_CS(indiv, Side_chain_carbons, average_CS, i)    
    
    average_CS, experimental_CS = M.set_up('../Results.dat',  'synth_shifts_'+str(i)+'.txt', 'R')


    
    #start the algorithm 
    M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, number_of_simulations_global, verbose, Max_gen_number_global)
    
    os.system('mv Solutions.txt algo_solutions_all/Solutions_'+str(i)+'.txt ')
    os.system('mv synth_shifts_'+str(i)+'.txt synth_cs/synth_shifts_'+str(i)+'.txt' )

    #now we do a run where the only input is the cg2 and cd
    for cs in experimental_CS:
        if cs != str(12):
            if cs != str(24):
                experimental_CS[cs] ='none'

    #start the algorithm 
    M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, number_of_simulations_global, verbose, Max_gen_number_global)    
    os.system('mv Solutions.txt algo_solutions_cg2_cd/Solutions_'+str(i)+'.txt ')


pic.dump( synth_pop, open( "synth_pop.p", "wb" ) )
pic.dump( synth_pop_CS, open( "synth_cs.p", "wb" ) )
print 'It took', time.time()-start, 'seconds.'
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







