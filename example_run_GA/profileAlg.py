#! /usr/bin/python
# written by Lucas Siemons

import cs2rotamer as c
import argparse
import time
from cs2rotamer.GeneticAlgorithm import *
import cProfile

if __name__ == '__main__':
    
    #====================================================================
    # Arguments
    #====================================================================

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=" ")
    parser.add_argument("-f", type=str, default='InputFile.txt',help="configuration file containing the algorithm configuration")    
    args = parser.parse_args()
    config_file = args.f

    #====================================================================
    # Run the program
    #====================================================================

    i = 1
    params = c.ReadInputFile(config_file)
    calRes = params.calcRes
    expshifts = params.expShifts
    average_CS, experimental_CS = c.setup_Fitted_DFT(calRes,  expshifts,'R', params)
    
    Mean_CS = average_CS
    Measured_CS = experimental_CS

    #the first gen is produced randomly
    Generation_1 = First_Gen(i, params) 


    IO.Write_out(1,Generation_1, i, 'Generations', params)

    #now we work out the shifts 
    cProfile.run('cs.Computed_shifts(Generation_1,  Mean_CS,params)')

    
