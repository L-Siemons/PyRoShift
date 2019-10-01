'''
Some utility functions
'''

import os.path


def cm2inch(value):
    '''
    convert from cm to inches
    '''
    return value / 2.54


def read_pdb_distribution(residue):
    '''
    read the pdb distributions
    '''

    if residue.lower() == 'ile':
        path = os.path.join(
            os.path.split(__file__)[0],
            "resources/ile_top8000_distribution.dat")
        file_ = open(path)
        data = file_.readlines()
        file_.close()
        pops = {a.split()[0]: float(a.split()[1]) for a in data}

    return pops
