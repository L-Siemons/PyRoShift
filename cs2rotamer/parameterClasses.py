import os
from cs2rotamer import distanceFunctions as dist


class ReadInputFile:

    '''
    This class contains a bunch of parameters that are read from the 'Input file'
    This is often refered to as the parameters class

    Parameters
    ==========
    InputFile : str
        name of the input file

    '''

    def __init__(self, inputFile):
        
        '''
        This function reads the input file and based on the configuration 
        sets the parameters for the algorithm 

        The parameters that are set are: 

        self.disMeasure : func
            The distance measure used to drive the algoruthm
            currenlty can only be Euclid_dis

        self.pop_size : int
            This is the population size of the GA

        self.Mutation_rate : int
            mutation rate for the GA
        
        self.verbose : int
            verbosity 
        
        self.Max_gen_number : int
            maximum generations that the GA runs for 
        
        self.write_out_location : str
            where the GA writes out the results 
        
        self.usedCC : list
            the states the GA is set to determine 
        
        self.cList : list
            ids for all the carbon atoms 
        
        self.Side_chain_carbons : list
            the used side chain carbon ids 
        
        self.CalcRes : str
            File with the Calculated DFT results
        
        self.ExpShifts :str 
            file with the experimental/observed chemical shifts
        
        self.zero_chance : int
            frequency of 0s in the first generation
        
        self.write_gen : int
            write out the generaitons 
        
        self.sol_name : str
            name of the solution file
        '''

        f = open(inputFile, 'r')
        
        #set a defult value here so it doesn't need to be in the input file 
        self.writeDistance = None
        
        for line in f.readlines():
            
            if line[0] != ';':
                split = line.split()    
                
                #these are the input feilds 
                
                #distance measure
                if split[0] == 'dist_func':
                    if split[1] == 'Euclid_dis':
                        self.disMeasure = dist.Euclid_dis
                    else:
                        print 'The distance measure does not exist'
                        sys.exit()
                
                #population size
                if split[0] == 'pop_size':
                    self.pop_size = int(split[1])
                
                if split[0] == 'Mutation_rate':
                    self.Mutation_rate = int(split[1])
                
                if split[0] == 'verbose':
                    self.verbose = int(split[1])
                
                if split[0] == 'Max_gen_number':
                    self.Max_gen_number = int(split[1])
                
                if split[0] == 'write_out_location':
                    self.write_out_location = split[1]
                
                if split[0] == 'residue':
                    if split[1] == 'ILE':
                        #self.usedCC = [[-175,167],[-167,65],[-58,-60],[60,174], [-64,169]]
                        self.usedCC = [['t','t'], ['t','p'],['m','m'],['p','t'],['m','t']]
                        self.cList = ['1','2','4','5','6','8','10','12','24']
                        self.Side_chain_carbons = ['4','5', '6','12','24']
                        self.number_of_states = len(self.usedCC)
                
                if split[0] == 'CalcRes':

                    if  split[1].lower() == 'default':
                        self.calcRes =os.path.join(os.path.dirname(__file__),'resources/PdbFittedChemicalShifts.dat') 
                    else:    
                        self.calcRes = split[1]
                    
                if split[0] == 'ExpShifts':
                    self.expShifts = split[1]
                
                if split[0] == 'zero_chance':
                    self.zeros = int(split[1])

                if split[0] == 'write_gen':
                    self.write_gen = int(split[1])
                
                if split[0] == 'sol_name':
                    self.sol_name = split[1]

                if split[0] == 'write_distance':
                    if split[1] == 'None':
                        self.writeDistance = None
                    else:
                        self.writeDistance = split[1]

        #this is curretly a defult value
        self.Close_enough_to_measured = 0.000000001
        self.number_of_simulations = 1

        f.close()