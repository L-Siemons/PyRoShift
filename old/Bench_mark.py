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
import Master_2 as M
import multiprocessing
import itertools


#for a time stamp
start = time.time()

#this meas we will be writing everything out!
verbose_text = 0

#write out each generation 1 = yes, 0 = now
write_gen = 0

#chance of zero when making the initial individual state being assigned a zero populaiton
zero_chance = 5

#this is the number of times you would like to run the simulation
data_point_number = 1000


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

pop_size_global = 20  #size of initial population

#i think the filtering in the cross over stage give a preference to higher mutation rates
#but ... again maybe its just an idea, not tested because when we filter them out a bad mutation will be disregarded
#so we only take good mutations ...?
#reddit seems to suggest that 1-5% is best but 10 also seems good aswel as much higher  
#when i raised to to 80% it took 2/3 of the generations (~200)
#when it was 95 it took 170 gen


Mutation_rate_global = 60.0
number_of_simulations_global = 1
simulation_Number = 0
simulation_Path ='Generations/Simulation_'
number_of_states = len(usedCC)

#these are used to determine the end conditions
Max_gen_number_global = 100
Close_enough_to_measured = 0.000001

#========================================================================================
#functions used in the benchmarking runs
#========================================================================================
def get_synth_CS(rotamer_pops, sideChain_atoms, mean_CS, id, dir):

    cs_input = open(dir+'synth_shifts_'+str(id)+'.txt', 'w')
    for atom in sideChain_atoms:
        CS_sum = 0.0
        for i2 in range(len(rotamer_pops)):
            CS_sum = CS_sum + rotamer_pops[i2] * float(mean_CS[atom, str(usedCC[i2][0]), str(usedCC[i2][1])])
        
        #now CS sum is the total chemical shift for that atom.
        
        cs_input.write(str(atom)+'  '+str(CS_sum)+'\n')
    cs_input.close()

def use_cd_cg2_only(experi_CS):
    #now we do a run where the only input is the cg2 and cd
    for cs in experi_CS:
        if cs != str(12):
            if cs != str(24):
                experi_CS[cs] ='none'
    return experi_CS

def get_off_CS_dicts(CS_file):
    ''' this function gets chemical shifts based on the backbone CS distribution
    so now the value of the CS is taken from a normal distribution '''
    
    # key = [atom , x1, x2]
    res_cs = {}
    res_std =  {}
    
    check_1 = 0
    check_2 = 0
    f = open(CS_file, 'r')  
    for i in f:
        split = i.split()
        try:
            if split[0][0] != '#':
                
                if check_1 == 1:
                    
                    res_cs[split[0], split[1], split[2]] = float(split[3])
                    res_std[split[0], split[1], split[2]] = float(split[4])
                
                                    
                if split[0] == 'Mean':
                    check_1 = 1
                            
        
        except IndexError:
            pass
        
    return res_cs, res_std
        
def get_one_off_cs(cs, std):
    
    off_cs = {}
    for i in cs:
        off_cs[i] = np.random.normal(cs[i], std[i], 1)[0]
    
    
    return off_cs
        
#========================================================================================
#bench marking runs
#========================================================================================

def pure_states(dir): 
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    try:
        os.mkdir(dir)
    except OSError:
        pass

    for i, element_1 in enumerate(usedCC):
    
        print i 
        empty = []
        for indx, element in enumerate(usedCC):    
            empty.append(0.0)
    
        empty[i] = 1.0
        indiv = empty
        synth_pop[i] = indiv
    
        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        #make the input file
        print 'run: ' + str(i)
        
        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
        average_CS, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')

        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')
       
       
        experimental_CS = use_cd_cg2_only(experimental_CS)
        
        
        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
        pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
        pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    return ''

def two_states(dir):
    try:
        os.mkdir(dir)
    except OSError:
        pass
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    #number of states you want to have a population of 
    zero_states = 3
    for i in range(data_point_number):
    
        #this section makes the individual.
        
        indiv = []
        zeros = rn.sample(range(0,len(usedCC)), zero_states)
        for i3 in range(len(usedCC)):
            if i3 in zeros:
                indiv.append(0.0)
            else:
                indiv.append(rn.randint(0, 1000))
        indiv_len = float(sum(indiv))
        for i4 in range(len(indiv)):
            indiv[i4] = float(indiv[i4])/indiv_len

        synth_pop[i] = indiv

        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
    
        average_CS, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')


    
        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')


        experimental_CS = use_cd_cg2_only(experimental_CS)

        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)

        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
    pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
    pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    
def all_mixed_states(dir):
    try:
        os.mkdir(dir)
    except OSError:
        pass
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    #number of states you want to have a population of 
    zero_states = 0
    for i in range(data_point_number):
    
        #this section makes the individual.
        
        indiv = []
        zeros = rn.sample(range(0,len(usedCC)), zero_states)
        for i3 in range(len(usedCC)):
            if i3 in zeros:
                indiv.append(0.0)
            else:
                indiv.append(rn.randint(0, 1000))
        indiv_len = float(sum(indiv))
        for i4 in range(len(indiv)):
            indiv[i4] = float(indiv[i4])/indiv_len

        synth_pop[i] = indiv

        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
    
        average_CS, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')


    
        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')


        experimental_CS = use_cd_cg2_only(experimental_CS)

        #start the algorithm 
        M.algo(average_CS, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)

        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
    pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
    pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    
def off_cs_two_states(dir):
    
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    try:
        os.mkdir(dir)
    except OSError:
        pass
    
  #make the grid of individuals:
    
    diff = 0.1
    list = [[1.0, 0.0]]
    
    sum = 1
    while sum > 0.0:
        
        sum =round( sum - diff, 3)
        
        list.append([ round(sum, 3), round(1-sum, 3)  ])


    for i in list:
        for a in range(len(usedCC)-2):
            i.append(0.0)
        
    #get all the permutations of this list
   
    list_2 = []
    for i in list:
        entry =  [x for x in itertools.permutations(i)]
        list_2 = list_2 + entry
    
    reduced_list = set(list_2)
    number_of_repeats = data_point_number / len(reduced_list)
    
    all_indivs = []
    
    for i in range(number_of_repeats):
        for i2 in reduced_list:
            all_indivs.append(i2)
    
    #normalise i and make it the individual    
    norm_list = []
    for i in all_indivs:
        sum = 0.0
        for i2 in i:
            sum =  sum +i2
        
        entry = []
        for i2 in i:

            entry.append(float(i2)/sum)
        
        norm_list.append(entry)
    
    cs_mean, cs_std = get_off_CS_dicts('Results.dat')
    for i, element_1 in enumerate(norm_list):

        indiv = norm_list[i]
        synth_pop[i] = indiv
    
        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        average_CS = get_one_off_cs(cs_mean, cs_std)

        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
    
        average_CS_2, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')
        
        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')


        experimental_CS = use_cd_cg2_only(experimental_CS)

        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)

        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
    pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
    pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    return ''

def off_cs_pure_states(dir):
    
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    try:
        os.mkdir(dir)
    except OSError:
        pass
    
  #make the grid of individuals:
    
    diff = 0.1
    list = [[1.0]]
    for i in list:
        for a in range(len(usedCC)-1):
            i.append(0.0)
        
    #get all the permutations of this list   
    list_2 = []
    for i in list:
        entry =  [x for x in itertools.permutations(i)]
        list_2 = list_2 + entry
    
    reduced_list = set(list_2)
    number_of_repeats = data_point_number / len(reduced_list)
    
    all_indivs = []
    
    for i in range(number_of_repeats):
        for i2 in reduced_list:
            all_indivs.append(i2)
    
    #normalise i and make it the individual    
    norm_list = []
    for i in all_indivs:
        sum = 0.0
        for i2 in i:
            sum =  sum +i2
        
        entry = []
        for i2 in i:

            entry.append(float(i2)/sum)
        
        norm_list.append(entry)
    
    cs_mean, cs_std = get_off_CS_dicts('Results.dat')
    for i, element_1 in enumerate(norm_list):

        indiv = norm_list[i]
        synth_pop[i] = indiv
    
        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        average_CS = get_one_off_cs(cs_mean, cs_std)

        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
    
        average_CS_2, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')
        
        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')


        experimental_CS = use_cd_cg2_only(experimental_CS)

        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)

        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
    pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
    pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    return ''

def off_cs_all_states(dir):
    
    synth_pop = {}
    synth_pop_CS = {}
    grid_search_res = {}
    #number of states you want to have a population of 
    zero_states = 0
    try:
        os.mkdir(dir)
    except OSError:
        pass
    
    cs_mean, cs_std = get_off_CS_dicts('Results.dat')
    for i in range(data_point_number):
        
            
        #this section makes the individual.
        
        indiv = []
        zeros = rn.sample(range(0,len(usedCC)), zero_states)
        for i3 in range(len(usedCC)):
            if i3 in zeros:
                indiv.append(0.0)
            else:
                indiv.append(rn.randint(0, 1000))
        indiv_len = float(sum(indiv))
        for i4 in range(len(indiv)):
            indiv[i4] = float(indiv[i4])/indiv_len


        synth_pop[i] = indiv
    
        input_name = dir+'synth_shifts_'+str(i)+'.txt'
    
    
        #ok here I have a dirty fix because I need to use the mean_CS file before I can read In 
        #the shifts so to avoid re-writing the function I call it twice and redefine the input the second time ... sorry
        #do the set up here
        average_CS, experimental_CS = M.set_up('Results.dat',  'Measured_shifts.txt', 'R')

        average_CS = get_one_off_cs(cs_mean, cs_std)

        get_synth_CS(indiv, Side_chain_carbons, average_CS, i, dir)    
    
        average_CS_2, experimental_CS = M.set_up('Results.dat',  dir+'synth_shifts_'+str(i)+'.txt', 'R')
        
        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)
        os.system('mv '+dir+'Solutions.txt  '+dir+'sol_'+str(i)+'_all.txt')


        experimental_CS = use_cd_cg2_only(experimental_CS)

        #start the algorithm 
        M.algo(average_CS_2, experimental_CS, write_gen, pop_size_global, Mutation_rate_global, 
              number_of_simulations_global, verbose_text, Max_gen_number_global, dir)

        os.system('mv '+dir + 'Solutions.txt  '+dir+'sol_'+str(i)+'_methly.txt')
    pic.dump( synth_pop, open( dir+"synth_pop.p", "wb" ) )
    pic.dump( synth_pop_CS, open( dir+"synth_cs.p", "wb" ) )
    return ''
    
#========================================================================================
#executing the runs
#========================================================================================


#define the number of cores
pool = multiprocessing.Pool( 2 )


#define the jobs
jobs = []

jobs.append( pool.apply_async( pure_states, args=(['mean_cs_pure_states/']) ))
jobs.append( pool.apply_async( two_states, args=(['mean_cs_two_states/']) ))
jobs.append( pool.apply_async( all_mixed_states, args=(['mean_cs_all_states/']) ))

jobs.append( pool.apply_async( off_cs_pure_states, args=(['off_cs_pure_states/']) ))
jobs.append( pool.apply_async( off_cs_two_states, args=(['off_cs_two_states/']) ))
jobs.append( pool.apply_async( off_cs_all_states, args=(['off_cs_all_states/']) ))


#do the jobs
for i in jobs:
    i.get()
