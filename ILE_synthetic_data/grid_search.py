#! /usr/bin/python
# written by Lucas Siemons

import os
import sys
import re
import random as rn 
import numpy as np 
import time
import pickle as pic
sys.path.append("..") 
import Master_2 as M
import matplotlib.pyplot as plt
import argparse
from scipy.stats import norm
import datetime

from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


#for a time stamp
start = time.time()
ts = datetime.datetime.now()

usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]
cut_off = 2.0
#these are the carbon identifiers
cList = ['1','2','4','5','6','8','10','12','24']
Side_chain_carbons = ['4','5', '6','12','24']


grid_search_cuttoff = 1
data_point_number = 1000

#functions


if __name__ == '__main__':
    
    input_pops = pic.load(open('synth_pop.p', 'r'))
    average_CS, experimental_CS = M.set_up('../Results.dat',  '../Measured_shifts.txt', 'A')
    input_CS = M.Computed_shifts(input_pops, average_CS)
    grid_cs = pic.load(open('../grid_search/grid_cs.p','r'))
    grid_pop =  pic.load(open('../grid_search/grid_pop.p','r'))
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(input_pops)).start()
    
    best_grid_pop = {}
    within_cut_off = {}
    gird_dist = {}
    best_dist = {}
    
    for indiv in input_pops:
        pbar.update(indiv+1)
        target = input_CS[indiv]
        closest_pop = []
        smallest_dist = 100.0
        for populations in grid_cs:
            
            dist = np.linalg.norm(np.array(input_CS[indiv])-np.array(grid_cs[populations]))**2
            if smallest_dist >= dist:
                smallest_dist = dist
                closest_pop = grid_pop[populations]
            if dist <= cut_off:
                try:
                    within_cut_off[indiv].append(grid_pop[populations])
                except KeyError:
                    within_cut_off[indiv] = []
                    within_cut_off[indiv].append(grid_pop[populations])
        best_dist[indiv] = smallest_dist
        best_grid_pop[indiv] = closest_pop
    pbar.finish()
    tf = datetime.datetime.now()
    print tf - ts
    pic.dump( best_grid_pop, open('best_grid_pops.p','w'))
    pic.dump( within_cut_off, open('within_cut_off_grid.p','w'))
    pic.dump( best_dist, open('best_dist_grid.p','w'))

    
    
    
    
    
    
    
    
    
    
    
    







