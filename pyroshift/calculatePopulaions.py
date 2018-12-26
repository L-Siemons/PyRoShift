'''
This is a module for the functions that deal determining populations
from chemical shifts
'''

import numpy as np
from scipy.optimize import least_squares
from .fileIO import Input, Output

class Isoleucine(Input, Output):

    '''
    This is the class for the user and calculates the rotamer distributions
    using a linear approach.

    Note that in order for this method to be effective all 5 chemical shifts must be
    provided and that should be referenced to DSS (not the TMS).

    This class initializes two classes to deal with the inputs and outputs for this method.
    The methods and attributes from these are documented in their respective class
    in pyroshift/FileIo.py



    Attributes:
    ===========
    atoms : list
        - contains all the names for the side-chain

    populations : nested dictionary
    - key is the residue
    - second key is the population
    - entry is a tuple containing floats :
        - (population, error)
    '''

    def __init__(self, file, shift_matrix='default', ref_opt_file='default'):

        '''
        This collects all the required inputs and initializes the classes
        Input and Output.

        The attributes these two classes create are documents under their respective class
        '''

        #read everything in!
        Input.__init__(self, file, shift_matrix=shift_matrix, ref_opt_file=ref_opt_file)
        Output.__init__(self,)
        self.atoms = ['ca', 'cb', 'cg1', 'cg2', 'cd1']

    def eq(self, state_matrix, pops_vector):
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

    def resid(self, pops_vector, cs_vector, matrix):
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
        resdidual = cs_model - cs_vector
        restraint = 1e4*(np.sum(pops_vector) - 1.)
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
        lb = np.zeros(5)
        ub = lb+1

        start_pops = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
        args = [cs_vector, state_matrix]

        result = least_squares(self.resid, start_pops, bounds=(lb, ub), args=args)
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

        for i, j in zip(pops_list, self.state_order[sse]):
            val = np.mean(i)
            std = np.std(i)
            pops_dict[j] = (val, std)

        return pops_dict

    def calc_pops_for_all_residues(self, ):
        '''
        Calculate populations for all residues that are specified in the the input file.
        We do this by solving the linear system

        C_states * P = C_exp

        Under the constraint sum(P_i) = 1. where;
        - C_states is a matrix containing the chemical shifts for each state
        - P is a vector containing the populations
        - C_exp is a vector with the experimental populations

        This method defined self.populations

        Returns:
        =======
        pops : nested dictionary
            - key is the residue
            - second key is the population
            - entry is a tuple containing floats :
                - (population, error)
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
