'''
This is a module for the functions that deal with reading and writing
files
'''

import re
import sys
import os.path
import numpy as np

def strip_comments(string):
    '''
    remove # comments from strings

    Args:
    =====
        string : str
    '''

    string = str(string)
    return re.sub(r'(?m)^ *#.*\n?', '', string)

class Input(file):
    '''
    This reads in the chemical shifts
    and stores them

    Args:
    =====
        file : str
            - file containing chemical shifts

        ref_opt_file : str
            - file containing the optimizations for the shielding tensors

        shift_matrix : str
            - path to the file with the chemical shift matrices, unless you
              are super sure 'default' is highly recommended

        ref_opt_file : str
            - path to the file with the shielding tensor tensor optimizations, unless you
              are super sure 'default' is highly recommended


    Attributes:
    ===========
        sse : dictionary
            keys are the residues. Entries are the backbone conformations
            So far the following are implemented
                - Alpha helix (A)
                - Beta sheet (B)
                - Random Coil (R)

        shifts : dictionary
            keys are the residues and the entries are nested dictionaries
            where the keys are the atom name (lower case) and the entires are the
            chemical shifts.

        opts : list of numpy arrays
            This list contains the optimizations for the reference shielding
            tensors

        shifts_matrix : dictionary
            The key is the backbone conformation and the entries are the matrices
            containing the chemical shifts for each state

        state_order : list
            a list of the states in the shifts_matrix and their order
    '''

    def __init__(self, shift_file, shift_matrix='default', ref_opt_file='default'):

        '''
        Read in the input file and set up the correct attributes

        Args:
        =====
        shift_file : str
            - file containing the experimental observed 13C chemical shifts

        shift_matrix : str
            - file containing the chemical shifts for each state

        ref_opt_file : str
            - file containing the optimizations for the chemical shift matrix
        '''

        sse_map = {}
        sse_map['a'] = 'alpha'
        sse_map['b'] = 'beta'
        sse_map['r'] = 'coil'

        self.sse = {}
        self.shifts = {}
        self.opts = []
        self.shift_matrix = {}
        self.state_order = {}

        f = open(shift_file)
        for line in f.readlines():

            #get rod of empty lines
            line = strip_comments(line)
            s = line.split()
            if s != []:
                res = s[0]
                atom = s[1].lower()
                if atom == 'ca':
                    current_sse = s[3].lower()
                    self.sse[res] = sse_map[current_sse]

                if res not in self.shifts:
                    self.shifts[res] = {}

                try:
                    self.shifts[res][atom] = float(s[2])
                except ValueError:

                    sys.tracebacklimit = 0
                    str1 = 'The input file (%s) contains a non-float\n' %(file)
                    str2 = 'in the Third column. Please correct this to proceed!'
                    raise Exception(str1+str2)
        f.close()

        resources = os.path.join(os.path.split(__file__)[0], "resources/")

        if shift_matrix == 'default':
            shift_matrix = resources + 'CS_DFT_surf_average_aceIleNmePdb_final.cs'

        matrix_file = open(shift_matrix)
        for line in matrix_file:
            s = line.split()
            state = s[1]
            sse = s[0]

            try:
                if state not in self.state_order[sse]:
                    self.state_order[sse].append(state)
            except KeyError:
                self.state_order[sse] = [state]

            if sse not in self.shift_matrix:
                self.shift_matrix[sse] = []

            self.shift_matrix[sse].append([float(a) for a in s[2:]])

        for backbone in self.shift_matrix:
            self.shift_matrix[backbone] = np.array(self.shift_matrix[backbone])

        if ref_opt_file == 'default':
            ref_opt_file = resources +  'compiled_ref_opts.dat'

        ref_opts = open(ref_opt_file)
        for line in ref_opts:
            s = line.split()
            self.opts.append(np.array([float(a) for a in s[2:]]))
        ref_opts.close()

class Output():
    '''
    This class contains functions for writing out the results from the
    chemical shift calculations.

    In general it is intended to be used as part of calculatePopulaions.isoleucine
    and uses it's attribute isoleucine.populations

    There are two output modes: print to screen and write to file

    Attributes:
    ===========
    output_lines : list
        - list of lines containing the populations and the relative error
    '''

    def __init__(self):
        '''
        Initialize the class
        '''

    def generate_lines(self,):
        '''
        This function compiles the results from the calculation
        into a set of lines that can be written to a file or printed
        '''

        states = list(set(self.state_order['alpha']))
        top_line = '# residue'  + '  ' + '  '.join(states)
        lines = [top_line]

        for res in self.populations:
            top_line = '#'
            sse = self.sse[res]
            top_line = '#'  + '  ' + '  '.join(self.state_order[sse])
            total_list = []

            for pop in states:
                val =  self.populations[res][pop][0]
                err =  self.populations[res][pop][1]
                total = '%0.3f+/-%0.3f' % (val, err)
                total_list.append(total)

            lines =lines + [res +  '    '+ '    '.join(total_list)]

        self.output_lines = lines

    def print_lines(self):
        '''
        prints the calculation results to stdout (terminal screen)
        '''

        #generate the lines if need be
        if hasattr(self, 'output_lines') == False:
            self.generate_lines()


        for i in self.output_lines:
            print i

    def write_lines(self, file):
        '''
        writes the calculation results to a file

        Args:
        =====
        file : str
            - name of output file
        '''

        #generate the lines if need be
        if hasattr(self, 'output_lines') == False:
            self.generate_lines()

        out = open(file, 'w')
        for i in self.output_lines:
            out.write(i)
            out.write('\n')
        out.close()

    def write_latex_pops(self, file):
        '''
        writes the populations into a latex table.
        
        Args:
        =====
        file : str
            - name of output file
        '''

        map_ = {}
        map_['m'] = 'g_{m}'
        map_['t'] = 't'
        map_['p'] = 'g_{p}'



        latex_header = '\\begin{table}[th]\n\\centering\n\\begin{tabular}{|l |l |l |l |l |l |l |l  |l |l| }\n\\hline\n'
        out = open(file,'w')
        out.write(latex_header)

        for i in self.output_lines:
            i = i.lower().replace('ile', '')
            j = i.replace('#', '')
            k = ' & '.join(j.replace('+/-','\\pm' ).split()) + '\\\\\n'
            out.write(k)
        latex_footer = '\\end{tabular}\n\\caption{}\\label{tab:RCpops}\n\\end{table}\n'
        out.write(latex_footer)
        out.close()
