#! /usr/bin/python
# written by Lucas Siemons

import numpy as np

#this fiel contains some python functions which may be useful for analysis.

def get_chi_2(J_exp, J_trans, J_g):
    
    '''
    
    this is taken from flemming's paper on ILE 
    3_J_exp is the measured 3J(Ca,Cd) scalar coupling 
    J_trans and J_g are the theoretical ones for the trans and g- states
    
    based on Flemming's ILE paper
    
    '''
    
    top = J_exp - J_trans
    bottom = J_g - J_trans 
    pg_m  = np.divide(top, bottom)
    p_trans = 1- pg_m
    return p_trans, pg_m