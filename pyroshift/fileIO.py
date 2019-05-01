'''
This is a module for the functions that deal with reading and writing
files
'''

import re
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from math import pi
import re
from . import util
from . import exceptions

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

        self.sse = {}
        self.shifts = {}
        self.opts = []
        self.shift_matrix = {}
        state_order = {}

        alpha_list = ('a', 'alpha', 'h', 'helix')
        beta_list = ('b', 'beta', 'sheet', 'strand')
        coil_list = ('c', 'coil', 'r')

        f = open(shift_file)
        for line in f.readlines():

            #get rod of empty lines
            line = strip_comments(line)
            s = line.split()
            if s != []:
                res = s[0]
                atom = s[1].lower()

                # Here we read in the secondary structure. 
                if atom == 'ca':
                    current_sse = s[3].lower()
                    if current_sse in alpha_list:
                        self.sse[res] = 1. 
                    elif current_sse in beta_list:
                        self.sse[res] = 0.
                    elif current_sse in coil_list:
                        self.sse[res] = 0.5

                    # Do a few checks
                    else:
                        check = True
                        try:
                            
                            self.sse[res] = float(current_sse)
                            if 1. <= self.sse[res]:
                                check = False
                            if 0. >= self.sse[res]:
                                check = False

                        except ValueError:
                            check = False

                        if check == False:
                            print res, current_sse
                            raise exceptions.Secondry_Structure_Error(res, current_sse)

                    #self.sse[res] = sse_map[current_sse]

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

        if ref_opt_file == 'default':
            if shift_matrix == 'default':
                ref_opt_file = resources +  'ref_opts_5_states.dat'

            elif shift_matrix == 'default_4_states':
                ref_opt_file = resources +  'ref_opts_4_states.dat'

        ref_opts = open(ref_opt_file)
        for line in ref_opts:
            s = line.split()
            self.opts.append(np.array([float(a) for a in s[2:]]))
        ref_opts.close()

        if shift_matrix == 'default':
            shift_matrix = resources + 'CS_DFT_ile_5_states.cs'

        if shift_matrix == 'default_4_states':
            shift_matrix = resources + 'CS_DFT_ile_4_states.cs'

        matrix_file = open(shift_matrix)
        for line in matrix_file:
            s = line.split()
            state = s[1]
            sse = s[0]

            try:
                if state not in state_order[sse]:
                    state_order[sse].append(state)
            except KeyError:
                state_order[sse] = [state]

            if sse not in self.shift_matrix:
                self.shift_matrix[sse] = []

            self.shift_matrix[sse].append([float(a) for a in s[2:]])
        matrix_file.close()

        for backbone in self.shift_matrix:
            self.shift_matrix[backbone] = np.array(self.shift_matrix[backbone])
        
        #check all the state orders are the same
        #if this is not true we cannot add and subtract the matrices

        orders = [str(i[1]) for i in state_order.items()]
        identity_check = len(set(orders))
        if identity_check != 1:
            raise exceptions.ShiftMatrixStatesOrder(shift_matrix)
        else:
            self.state_order = state_order.items()[0][1]


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

    def generate_lines(self,order=False,strip=False):
        '''
        This function compiles the results from the calculation
        into a set of lines that can be written to a file or printed
        '''

        states = list(set(self.state_order))
        top_line = '# residue'  + '  ' + '  '.join(states)
        lines = [top_line]

        if order == False:
            keys = self.populations
        else:
            keys = [(''.join([i for i in a if not i.isdigit()]), int(re.sub("\D", "", a))) for a in self.populations]
            keys = sorted(keys, key=lambda x: x[1])
            keys = [i+str(k) for i,k in keys]

        for res in keys:
            top_line = '#'
            sse = self.sse[res]
            top_line = '#'  + '  ' + '  '.join(self.state_order)
            total_list = []

            for pop in states:
                val =  self.populations[res][pop][0]
                err =  self.populations[res][pop][1]
                total = '%0.3f+/-%0.3f' % (val, err)
                total_list.append(total)
            
            if strip == True: 
                res = re.sub("\D", "", res)
            lines =lines + [res +  '    '+ '    '.join(total_list)]

        self.output_lines = lines

    def print_lines(self, order=False, strip=False):
        '''
        prints the calculation results to stdout (terminal screen)
        '''

        #generate the lines if need be
        if hasattr(self, 'output_lines') == False:
            self.generate_lines(order=order, strip=strip)


        for i in self.output_lines:
            print i

    def write_lines(self, file, order=False, strip=False):
        '''
        writes the calculation results to a file

        Args:
        =====
        file : str
            - name of output file
        '''

        #generate the lines if need be
        if hasattr(self, 'output_lines') == False:
            self.generate_lines(order=order, strip=strip)

        out = open(file, 'w')
        for i in self.output_lines:
            out.write(i)
            out.write('\n')
        out.close()

    def write_latex_pops(self, file, name_tag=None,):
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

        if name_tag == None:
            tag = '\\caption{}'
        else: 
            tag = '\\caption{Chemical shift populations for %s.}' %(name_tag)



        latex_header = '\\begin{table}[th]\n%s\\label{}\n\\centering\n\\begin{tabular}{|l |l |l |l |l |l |l |l  |l |l| }\n\\hline\n' %(tag)
        out = open(file,'w')
        out.write(latex_header)

        for i in self.output_lines:
            i = i.lower().replace('ile', '')

            if '#' in i:
                j = i.replace('#', '')
                for entry in map_:
                    j = j.replace(entry, '$%s$'%(map_[entry]))
                j = j.replace('$$', '/')
            else:
                j = i 

            k = ' & '.join(j.replace('+/-','$\\pm$' ).split()) + '\\\\\n'
            out.write(k)
        latex_footer = '\\end{tabular}\n\\end{table}\n'
        out.write(latex_footer)
        out.close()

    def write_shift_tex_table(self, out, name_tag=None, strip=False,order=False):
        '''

        '''

        atom_map = {}
        atom_map['ca'] = 'C$_{\\alpha}$'
        atom_map['cb'] = 'C$_{\\beta}$'
        atom_map['cg1'] = 'C$_{\\gamma 1}$'
        atom_map['cg2'] = 'C$_{\\gamma 2}$'
        atom_map['cd1'] = 'C$_{\\delta}$'

        if name_tag == None:
            tag = '\\caption{}'
        else: 
            tag = '\\caption{Experimental and calculated chemical shifts for %s.}' %(name_tag)


        
        header = '\\begin{table}[th]\n%s\\label{}\n\\centering\n\\begin{tabular}{|l |l |l |l |l |}\n\\hline\n' %(tag)
        out_file = open(out, 'w')
        out_file.write(header)
        out_file.write('residue & atom & $\\delta_{calc}$ & $\\delta_{exp}$ & $|\\delta_{calc}-\\delta_{exp}|$ \\\\ \n')

        if order == False:
            keys = self.calc_shifts

        else:
            keys = [(''.join([i for i in a if not i.isdigit()]), int(re.sub("\D", "", a))) for a in self.calc_shifts]
            keys = sorted(keys, key=lambda x: x[1])
            keys = [i+str(k) for i,k in keys]


        for res in keys:
            for atom in self.atoms:
                calc = self.calc_shifts[res][atom][0]
                exp = self.shifts[res][atom]
                diff = abs(exp-calc)
                if strip == True:
                    res_entry = re.sub("\D", "", res)
                else:
                    res_entry = res
                line =  '%s & %s & %0.2f & %0.2f & %0.2f \\\\ \n' % (res_entry, atom_map[atom], calc, exp, diff)
                out_file.write(line)
        
        end =  '\\end{tabular}\n\\end{table}\n'
        out_file.write(end)

    def plot_radar(self, tex=False, top8000=False, size=(10,8), label_size=(12, 10),save=False, set_level=False):
        '''
        plot radar plots for each of the populations

        Args:
        =====

        tex : bool
            - use tex to make the axis, note this will require tex
        top8000 : bool or str 
            - plot the top8000 rotamer distribution. False or the 3 letter code for the residue.
        size : tuple 
            - size of the figure 
        label_size : tuple
            - text sizes for the axis
        save : bool 
            - save the figures instead of displaying them 

        '''

        matplotlib.rc('xtick', labelsize=label_size[0])
        matplotlib.rc('ytick', labelsize=label_size[1])
        matplotlib.rcParams['axes.linewidth'] = 2
    
        if tex == True:
            
            conv_dict = {}
            conv_dict['t'] = 't' 
            conv_dict['m'] = 'g_{m}' 
            conv_dict['p'] = 'g_{p}' 


        for res in self.populations:
            
            categories = self.state_order
            pops = [self.populations[res][a][0] for a in categories]
            N = len(categories)

            if top8000 != False:
                pdb_pops = util.read_pdb_distribution(top8000)
                pdb_pops = [pdb_pops[a] for a in categories]
            else:
                pdb_pops = []


            if tex == True:
                categories = ['/'.join([conv_dict[b] for b in list(a)]) for a in categories]
                categories = ['$%s$'%(a) for a in categories]

            max_pop = max(pops+pdb_pops)

            # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
            angles = [(n / float(N) *2* pi)for n in range(N)]

            angles += angles[:1]

            # Initialise the spider plot
            plt.figure(figsize=(util.cm2inch(size[0]), util.cm2inch(size[1])))
            plt.rc('text', usetex=True)
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams['text.latex.preamble'] = [
                            r'\usepackage{tgheros}',    # helvetica font
                            r'\usepackage{sansmath}',   # math-font matching  helvetica
                            r'\sansmath'                # actually tell tex to use it!
                            r'\usepackage{siunitx}',    # micro symbols
                            r'\sisetup{detect-all}',    # force siunitx to use the fonts
                        ]  


            ax = plt.subplot(111, polar=True)
            ax.yaxis.grid(linewidth=2)
            ax.xaxis.grid(linewidth=2)
            # Draw one axe per variable + add labels labels yet
            plt.xticks(angles[:-1], categories, color='k')

            # Draw ylabels
            ax.set_rlabel_position(angles[2]*(180./pi))
            deg_angles = [a*(180./pi) for a in angles[:-1]]
            ax.set_thetagrids(deg_angles, frac=1.24)
            

            if 0.75 < max_pop or set_level == 1.: 
                plt.yticks([0.25,0.5,0.75], ['0.25','0.5','0.75'], color="k")
                plt.ylim(0,1)
            elif 0.5 < max_pop <= 0.75 or set_level == 0.75:
                plt.yticks([0.25,0.5], ['0.25','0.5'], color="k")
                plt.ylim(0,0.75)          
            elif 0 < max_pop <= 0.5 or set_level == 0.5:
                plt.yticks([0.1, 0.2, 0.3, 0.4], ['0.1', '0.2', '0.3', '0.4'], color="k")
                plt.ylim(0,0.5)
            else:
                print 'you should not have got here ...'

            #plot rc 
            if top8000 != False:
                pdb_values = pdb_pops + [pdb_pops[0]]
                ax.plot(angles, pdb_values, linewidth=3, linestyle='solid',c='k')
                ax.fill(angles, pdb_values, 'k', alpha=0.1)            

            values = pops + [pops[0]]
            ax.plot(angles, values, linewidth=3, linestyle='solid',c='b')
            ax.fill(angles, values, 'b', alpha=0.1)

            if save == False:
                plt.show()
            else:
                name = res+'.pdf'
                plt.savefig(name,bbox_inches='tight')