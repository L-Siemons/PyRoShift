#! /usr/bin/python
# written by Lucas Siemons

import os
import sys
import glob
import re
import random as rn 
import numpy as np 
import pickle as pic
import matplotlib.pyplot as plt

b = int(sys.argv[1])

distances = {}


for i in range(b):
    
    if i % 5 == 0:
        sum = 0.0
        a = 1000.0
        file = open('Distances/Simulation_0/gen_'+str(i+1)+'.txt','r')
        counter = 0
        for line in file.readlines():
            counter = counter +1
            split = line.split()
        
            sum = sum + float(split[1])
        
            if float(split[1]) <= a:
                a = split[1]
            
        average = np.divide(sum, float(counter))    
        distances[i+1] = a
        file.close()
    
    
        plt.scatter(i+1, average, c='r')
    
for item in distances:
    plt.scatter(item, distances[item], c='b')
plt.ylabel('distance measure')
plt.xlabel('generations')
plt.show()