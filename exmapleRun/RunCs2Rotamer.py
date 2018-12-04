#! /usr/bin/python
# written by Lucas Siemons

import cs2rotamer as c
import argparse
import time


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

    params = c.ReadInputFile(config_file)
    calRes = params.calcRes
    expshifts = params.expShifts
    experimental_CS, SSE = c.ReadIn_measured_C_shifts(params.expShifts)
    average_CS = c.setup_Fitted_DFT(calRes,SSE, params)

    start = time.time()
    c.algo(average_CS, experimental_CS, params)
    end = time.time()
    print 'time taken: ', end - start



