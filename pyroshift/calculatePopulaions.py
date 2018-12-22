'''
This is a module for the functions that deal determining populations
from chemical shifts
'''

import numpy as np
from scipy.optimize import least_squares

import pyroshift
from .fileIO import Input

class isoleucine(Input):

    '''
    This is the class for the user and calculates the rotamer distributions
    using a linear approach.

    Attributes:
    ===========

    Methods:
    ========

    '''

    def __init__(self, file,shift_matrix='default', ref_opt_file='default'):

        '''
        This collects all the required inputs.
        '''

        #read everything in!
        Input.__init__(self,file, shift_matrix=shift_matrix, ref_opt_file=ref_opt_file)
        self.atoms = ['ca', 'cb', 'cg1', 'cg2', 'cd1']

    def eq(self,state_matrix, pops_vector):
        '''
        The equation for the linear system:

        shift_matrix * populations = observed shifts

        Args:
        =====
        state_matrix : np.array
            contains the chemical shifts for each state

        Returns:
        ========
        out : np.array
            chemical shifts corresponding to that
            population distribution


        '''
        return np.dot(state_matrix.T, pops_vector)

    def resid(self, pops_vector,cs_vector, matrix,):
        '''
        Calculate the residual for the chi squared. Here
        we also add a restraint so the populations sum to 1.

        Args:
        =====
        cs_vector : np.array
            observed chemical shifts
        pops_vector : np.array
            population distribution

        Returns:
        ========
        total : np.array
            this has the residuals and the restraint
        '''

        cs_model = self.eq(matrix, pops_vector)
        resdidual =  cs_model - cs_vector
        restraint = 1e4*(np.sum(pops_vector) - 1. )
        total = np.append(resdidual, restraint)
        return total

    def calc_pops(self, state_matrix, cs_vector):
        '''
        Least squares fitting to determine the populations

        Args:
        =====
        state_matrix : np.array
            matrix containing the chemical shifts for each state
        cs_vector : np.array
            observed chemical shift vector

        Returns:
        ========
        result.x : np.array
            populations from the fit.s
        '''
        lb=np.zeros(5)
        ub=lb+1

        start_pops = np.array([0.2,0.2,0.2,0.2,0.2])
        args = [cs_vector, state_matrix]

        result = least_squares( self.resid, start_pops, bounds=(lb,ub), args=args )
        return result.x

    def calc_pop_for_all_opts(self, state_matrix, cs_vector, sse):
        '''
        Here we calculate the populations using all the shielding tensor optimizations
        the final reported values are the mean and std

        Args:
        =====
        state_matrix : np.array
            matrix containing the chemical shifts for each state
        cs_vector : np.array
            observed chemical shift vector
        sse : str
            the backbone conformation

        Returns:
        ========
        pops_dict : dictionary
            keys are the states and entries are a tuple (value, std)
        '''

        pops_list = []
        for opt in self.opts:
            current_shift_matrix = state_matrix + opt
            pops = self.calc_pops(current_shift_matrix, cs_vector)
            pops_list.append(pops)

        pops_list = np.array(pops_list).T
        pops_dict = {}

        for i,j in zip(pops_list, self.state_order[sse]):
            val = np.mean(i)
            std = np.std(i)
            pops_dict[j] = (val, std)

        return pops_dict

    def calc_pops_for_all_residues(self, ):
        '''
        calculate populations for all residues
        '''

        print 'calculating the populations ...'
        pops = {}
        for residue in self.shifts:

            res_sse = self.sse[residue]
            shift_matrix = self.shift_matrix[res_sse]
            shifts = []

            for i in self.atoms:
                shifts.append(self.shifts[residue][i])
            shifts = np.array(shifts)
            pops[residue] = self.calc_pop_for_all_opts(shift_matrix, shifts, res_sse)

        self.populations = pops
        return pops