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
import matplotlib.pyplot as plt
import argparse
from scipy.stats import norm
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import metrics
from pandas.tools.plotting import parallel_coordinates
import pandas as pd
from sklearn.preprocessing import StandardScaler



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


synth_pop = {}
synth_pop_CS = {}
grid_search_res = {}

d = '''This script is here to do some basic comparison'''
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
parser.add_argument("-all", action='store_true', help="shows the RMSD comparing the input all output population districution as a whole")
parser.add_argument("-states", action='store_true', help="shows the RMSD of each state seperately")
parser.add_argument("-minima", action='store_true', help="shows the RMSD of each state seperately")

args = parser.parse_args()

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

  
def get_algo_res(path):
    synth_pop = {}
    for i in range(data_point_number):
        input_name = path+str(i)+'.txt'
        input = open(input_name, 'r')

        for pop in input.readlines():
            synth_pop[i] = pop.split()
    return synth_pop

def RMSD_all(a,b):
    sum = 0
    for indx, element in enumerate(a):
        dist = np.linalg.norm(float(a[indx])-float(b[indx]))**2
        sum = sum + dist
    rmsd = np.sqrt( np.divide(sum, float(len(a))) )
    return rmsd
    
def fit_hist(data):
    # Fit a normal distribution to the data:
    mu, std = norm.fit(data)

    # Plot the histogram.
    plt.hist(data, bins=25, normed=0, alpha=0.6, color='g')
    plt.hist(data, bins=25, normed=0, alpha=0.6, color='g')
    plt.xlabel('distance')
    plt.ylabel('frequency')
    '''# Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)'''
    plt.show()

def get_minima(raw_data):
    
    data  = np.array(raw_data)
    

    bandwidth = estimate_bandwidth(data, quantile=0.1, n_samples=20000)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(data)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    
    print("number of estimated clusters : %d" % n_clusters_)
    if n_clusters_ != 1:
        for i in range(n_clusters_):
            print 'cluster__ '+ str(i)
            T=X[ms.labels_==i].indices
            for ind in T:
                print terms[ind]
            print
    
    return n_clusters_
    
if __name__ == '__main__':
    
    grid_dist = pic.load(open('best_dist_grid.p','r'))
    best_grid_pop =  pic.load(open('best_grid_pops.p','r'))
    synth = pic.load(open("synth_pop.p", 'r'))
    
    
    if args.all:
        all_rmsd = []
        X = []
        counter = 0
        for indiv in synth:
            a1 = best_grid_pop[indiv]
            a2 = [x * 0.025 for x in a1] 
            
            all_rmsd.append(RMSD_all(synth[indiv], a2))
            X.append(counter)
            counter = counter + 1

        plt.xlabel('individual')
        plt.ylabel('overall RMSD')
        plt.scatter(X, all_rmsd)
        plt.show()
        fit_hist(all_rmsd)

    if args.states:
        
        total_data = []
        
        for indx, state in enumerate(usedCC):
            working = []
            
            
            
            for indiv in synth:
                working.append(abs(float(synth[indiv][indx]) - (float(best_grid_pop[indiv][indx] * 0.025) )))         
            total_data.append(working)
        
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

        # generate some random test data
        all_data = total_data

        # plot violin plot
        axes[0].violinplot(all_data,
                           showmeans=False,
                           showmedians=True)
        axes[0].set_title('violin plot')

        # plot box plot
        axes[1].boxplot(all_data, showfliers=False)
        axes[1].set_title('box plot')

        #axes[1].set_ylim([0.0,0.5])
        # adding horizontal grid lines
        for ax in axes:
            ax.yaxis.grid(True)
            ax.set_xticks([y+1 for y in range(len(all_data))])
            ax.set_xlabel('xlabel')
            ax.set_ylabel('average dist')

        # add x-tick labels
        plt.setp(axes, xticks=[y+1 for y in range(len(all_data))],
                 xticklabels = [str(x) for x in usedCC])
        plt.show()

    
    if args.minima:
        col_names = []
        for i in usedCC:
            col_names.append(str(i))
        grid_minima = pic.load(open("within_cut_off_grid.p", 'r'))
        for i in  grid_minima:
            a = get_minima(grid_minima[i])
            print a
            if a != 1:
            
                df = pd.DataFrame(grid_minima[i], index=range(len(grid_minima[i])), columns=col_names)
                entry_id = []
                counter = 0
                
                print 'making the entry id list'
                for entries in grid_minima[i]:
                    counter = counter + 1
                    entry_id.append(str(counter))
                df['Name'] = entry_id
                print df
            
                fig, ax = plt.subplots(1)
            
                #ax.plot(xdata, ydata)
                parallel_coordinates(df, 'Name')
                ax.legend().set_visible(False)
                fd.plot(color='red')
                print 'did I need to show?'
                plt.show()
                #parallel_coordinates(grid_minima[i], 'Name')
          
    
    
    
    
    
    