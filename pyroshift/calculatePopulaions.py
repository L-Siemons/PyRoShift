'''
This is a module for the functions that deal determining populations
from chemical shifts
'''

import numpy as np
from scipy.optimize import least_squares
from .fileIO import Input, Output
import random as rn

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

    def __init__(self, file, shift_matrix='default', ref_opt_file='default',force_constant='default'):

        '''
        This collects all the required inputs and initializes the classes
        Input and Output.

        The attributes these two classes create are documents under their respective class
        '''

        #read everything in! these are inuput and output functions
        Input.__init__(self, file, shift_matrix=shift_matrix, ref_opt_file=ref_opt_file)
        Output.__init__(self,)

        #define a few things
        self.atoms = ['ca', 'cb', 'cg1', 'cg2', 'cd1']
        if force_constant == 'default':
            self.force_constant = np.array([2./2.68,2./2.00,2./1.70,2./1.34,2./1.66])
        else: 
            self.force_constant = force_constant

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
        resdidual = (cs_model - cs_vector)*self.force_constant
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
            populations from the fits
        shifts : np array
        	contains the calculated chemical shifts
        '''
        lb = np.zeros(state_matrix.shape[0])
        ub = lb+1

        start_pops = [1./state_matrix.shape[0] for _ in range(state_matrix.shape[0])]
        start_pops = np.array(start_pops)
        args = [cs_vector, state_matrix]

        result = least_squares(self.resid, start_pops, bounds=(lb, ub), args=args)
        shifts = np.dot(result.x, state_matrix)

        return result.x, shifts

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
        shifts_dict : dictionary
        '''

        pops_list = []
        shift_list = []
        for opt in self.opts:
            current_shift_matrix = state_matrix + opt
            pops, shifts = self.calc_pops(current_shift_matrix, cs_vector)
            pops_list.append(pops)
            shift_list.append(shifts)

        pops_list = np.array(pops_list).T
        pops_dict = {}

        shift_list = np.array(shift_list).T
        shifts_dict = {}

        for i, j in zip(pops_list, self.state_order):
            val = np.mean(i)
            std = np.std(i)
            pops_dict[j] = (val, std)

        for i, j in zip(shift_list, self.atoms):
            val = np.mean(i)
            std = np.std(i)
            shifts_dict[j] = (val, std)

        return pops_dict, shifts_dict

    def calc_pops_for_all_residues(self, ):
        '''
        Calculate populations for all residues that are specified in the the input file.
        We do this by solving the linear system

        C_states * P = C_exp

        Under the constraint sum(P_i) = 1. where;
        - C_states is a matrix containing the chemical shifts for each state
        - P is a vector containing the populations
        - C_exp is a vector with the experimental populations

        This method defines self.populations and self.calc_shifts

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
        calc_shifts = {}
        for residue in self.shifts:

            res_sse = self.sse[residue]
            shift_matrix = self.shift_matrix['alpha']*res_sse + self.shift_matrix['beta']*(1.-res_sse)
            shifts = []

            for i in self.atoms:
                shifts.append(self.shifts[residue][i])
            shifts = np.array(shifts)
            current_pops, current_shifts = self.calc_pop_for_all_opts(shift_matrix, shifts, res_sse)

            pops[residue] = current_pops
            calc_shifts[residue] = current_shifts

        self.populations = pops
        self.calc_shifts = calc_shifts

        return pops

    def combine_along_chi_angles(self, mc_loop=1000):

        '''
        THis function projects the populations along each of the chi angles.
        This is done using a Monte Carlo approach. Note that at the end of this
        we normalise the populations so that they still sum to 1.

        Args:
        =====
        mc_loop : int
            number of cycles in the Monte Carlo loop, default is 1000

        Returns:
        =======
        combined : dict
            this is a nested dictionary where the keys are; residue, chi angle, state.
            the entry is tuple (value, entry)

        '''

        combined = {}
        for res in self.populations:
            states =  self.populations[res].keys()
            combined[res] = {}

            for current_angle, _ in enumerate(states[0]):
                single_angle_states = list(set([a[current_angle] for a in states]))
                chi_angle_name = str(current_angle+1)

                # learning point; Never use fromkeys() with mutables as default values!
                # as it behaves really weirdly; 1hr was spent right here!
                all_values = dict.fromkeys(single_angle_states, ())

                for _ in range(mc_loop):

                    current_values = dict.fromkeys(single_angle_states, 0.)

                    for state in  self.populations[res]:
                        current_state = state[current_angle]
                        mu = self.populations[res][state][0]
                        sig = self.populations[res][state][1]

                        #this is the mc bit
                        mc_value = rn.gauss(mu,sig)
                        if mc_value < 0.:
                            mc_value = 0.
                        elif mc_value > 1.:
                            mc_value = 1.

                        current_values[current_state] += mc_value

                    for i in all_values:
                        all_values[i] += (current_values[i],)

                total = 0.

                for i in all_values:
                    val = np.mean(all_values[i])
                    err = np.std(all_values[i])
                    total += val
                    all_values[i] = (val, err)

                #normalize the values and the error
                for i in all_values:
                    all_values[i] = [a/total for a in all_values[i]]

                combined[res][chi_angle_name] = all_values

        self.population_projections = combined
        return combined