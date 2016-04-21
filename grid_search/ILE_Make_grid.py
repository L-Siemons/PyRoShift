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
import datetime
sys.path.append("..") 
import Master_2 as M

#for a time stamp
start = time.time()

usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]

#for later use these are the backbone angles
usedPP = np.array([[-64, -45], [-71, -18],[-112, 8],[-105, -45],[-130, 157],[-126, 129],[-115, 126],[-100, 126],[-72, 132]])

#these are the carbon identifiers
cList = ['1','2','4','5','6','8','10','12','24']
Side_chain_carbons = ['4','5', '6','12','24']

#====================================================================
#Defining functions
#====================================================================

def set_up(Mean_CS_file,  Measured_CS_file):
    Mean_CS_1 = read_in_mean_results(Mean_CS_file)
    Measured_CS_1 = ReadIn_measured_C_shifts(Measured_CS_file)
    
    return Mean_CS_1, Measured_CS_1

if __name__ == '__main__':
    ts = datetime.datetime.now()
    print '=============================================================================='
    
    grid_pop = {}
    
    #very important to check these you you might waste lots of time!!
    step_size = 0.025
    sum_1 = []
    grid = np.empty([41,41,41,41,41])
    
    
    counter = 0
    print 'starting the huge loop ...'
    for idx in np.ndindex(grid.shape):
        if round(float(sum(idx)) * step_size, 5) == 1.0:
            counter = counter + 1
            grid_pop[counter] = list(idx)
    print 'finished the huge loop ...'
    print 'now saving...' 
        
    pic.dump( grid_pop, open( "grid_pop.p", "wb" ) )
    
    print 'reasing in mean CS'
    mean_cs = M.read_in_mean_results('../Results.dat')
    print mean_cs
    
    print 'calculating grid CS'

    grid_CS = M.Computed_shifts(grid_pop, mean_cs)

    print 'saving grid CS'   
    pic.dump( grid_CS, open( "grid_cs.p", "wb" ) )
    
    tf = datetime.datetime.now()
    te = tf - ts
    print 'whole script took:' + str(te)
    
    
    
    
    
    
    
    