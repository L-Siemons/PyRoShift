'''
This module contains functions that read and write files
'''
import numpy as np
from cs2rotamer import chemicalShifts as cs
import cs2rotamer.distanceFunctions as distFuncs

#===============================================================
# Read in file
#===============================================================

def ReadIn_measured_C_shifts(input_file):
    '''
    Reads in the measured/observed chemical shifts 

    Parameters
    =========
    input_file : str
        Name of the file containing the shifts
    
    Returns
    ======= 
    ReadIn_shifts : dict 
        Key - atom id 
        entry is the chemical shift

    '''

    ReadIn_shifts = {}
    inputshifts = open(input_file, 'r')
    for line in inputshifts.readlines():
        split = line.split()
        ReadIn_shifts[split[0]] = split[1]

    return ReadIn_shifts

def setup_Fitted_DFT(Mean_CS_file,  Measured_CS_file, sec_struct,params, read_measured_cs='yes'):

    '''
    This function reads in the appropriate calculated chemical shifts 
    Given the secondry structure.

    It does this be reading all the theoretical chemical shifts and then 
    selecting the right set based on the secondry structure given in the
    file with the measured chemical shifts. 

    These can be either R, A or B for random coil, alpha helix or beta sheet.
    
    Parameters
    =========
    Mean_CS_file : str
        file with the theoretical chemical shifts
    Measured_CS_file : str
        file with the measured chemical shifts 
    sec_struct :
        secondry structure 
    params : Parameters class 
    read_measured_cs : str
        select if to read in the measured chemical shift or not 

    Return
    ======
    selected_Mean_CS : dict
        contains the selected theoretical chemical shifts 
    Measured_CS_1 : dict
        contains the read chemical shifts
    '''

    alphaDFT = {}
    betaDFT  = {}
    rcDFT    = {}

    usedCC = params.usedCC
    dft = open(Mean_CS_file, 'r')
    #print 'states used:   ',  usedCC
    
    for line in dft.readlines():
        split = line.split()
        #print line
        if split[0] == 'alpha':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                alphaDFT['4', angles[0], angles[1]] = float(split[2])
                alphaDFT['5', angles[0], angles[1]] = float(split[3])
                alphaDFT['6', angles[0], angles[1]] = float(split[4])
                alphaDFT['12', angles[0], angles[1]] = float(split[5])
                alphaDFT['24', angles[0], angles[1]] = float(split[6])

        if split[0] == 'beta':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                betaDFT['4', angles[0], angles[1]] = float(split[2])
                betaDFT['5', angles[0], angles[1]] = float(split[3])
                betaDFT['6', angles[0], angles[1]] = float(split[4])
                betaDFT['12', angles[0], angles[1]] = float(split[5])
                betaDFT['24', angles[0], angles[1]] = float(split[6])

        if split[0] == 'rCoil':
            angles = [split[1][0], split[1][1]]
            if angles in usedCC:
                rcDFT['4', angles[0], angles[1]] = float(split[2])
                rcDFT['5', angles[0], angles[1]] = float(split[3])
                rcDFT['6', angles[0], angles[1]] = float(split[4])
                rcDFT['12', angles[0], angles[1]] = float(split[5])
                rcDFT['24', angles[0], angles[1]] = float(split[6])

    dft.close()

    if sec_struct == 'R':
        selected_Mean_CS = rcDFT
    elif sec_struct == 'A':
        selected_Mean_CS = alphaDFT
    elif sec_struct == 'B':
        selected_Mean_CS = betaDFT
    else:
        print '''ERROR: one of the options selected was not a recognised secondry structure type. 
        Pease check the set_up() function'''
        sys.exit()
    #read in the exprimental results
    if read_measured_cs == 'yes':
        Measured_CS_1 = ReadIn_measured_C_shifts(Measured_CS_file)
    else:
        print 'attention not readin measured CS'
        Measured_CS_1 = ''
    dft.close()

    return selected_Mean_CS, Measured_CS_1

#===============================================================
# write out file
#===============================================================

def Write_out(gen_number, populaton_dict, sim_number, dict_type, params):
    
    '''
    This function writes out the populations of each generation
    Parameters
    =========
    gen_number : int 
        generation
    populaton_dict : dict 
        key - individual id
        value - the rotamer distribution 
    sim_number : int 
        the simulation we are on, so far we only tested running a single one
    dict_type : str 
        label wheather generations or distances are written out
    params : Parameters class 
    '''

    out_direc = params.write_out_location

    if params.write_gen == 0:
        return None
        
    try:
        os.mkdir(out_dir+str(dict_type))
    except OSError:
        pass
    
    try:
        os.mkdir(out_dir+str(dict_type)+'/Simulation_'+str(sim_number)+'/')
    except OSError:
        pass
        
    dir_name = out_dir+str(dict_type)+'/Simulation_'+str(sim_number)+'/gen_'+str(gen_number)+'.txt'
    file = open(dir_name,'w')
    
    for i in range(params.pop_size):
        file.write(str(i)+'  ')
        try: 
            for index in range(len(populaton_dict[i])):
                file.write(str(populaton_dict[i][index])+ '  ')
            file.write('\n')
        except TypeError:
            #this just makes the function more flexible, as some dicts only have single vlaues in them
            file.write(str(populaton_dict[i])+ '\n')
    file.close() 

def Wite_solutions(best_indiv_id, final_gen, Mean_CS,Measured_CS,params): 
    '''
    This function writes out the solution from the GA

    Parameters
    =========  
    best_indiv_id : int 
        the id of the best individual
    final_gen : dict 
        the final generation of the GA 
    Mean_CS : dict
        the chemical shifts of the states
    Measured_CS : dict 
        the measured chemical shifts
    params : Parameters class     
    '''

    #here we write out the populations
    out_name = params.write_out_location+params.sol_name
    print 'outname:  ', out_name
    S_file = open(out_name, 'a')
    
    total = sum(final_gen[best_indiv_id])

    for i in range(len(final_gen[best_indiv_id])):

        S_file.write( str(np.divide(final_gen[best_indiv_id][i], total) ) + '   ' )
    S_file.write('\n')
    S_file.close()
    
    #the file for the best individual to be written to
    best_indiv_name = out_name + '_bestIndiv.txt'
    
    #this is just some reformatting so we can use Computed_shifts()
    dummyDict = {}
    dummyDict['dummy'] = final_gen[best_indiv_id]

    #calculate the chemical shift
    closest_shifts = cs.Computed_shifts(dummyDict, Mean_CS,params)
    
    print 'best f name: ', best_indiv_name
    print 'cs' , closest_shifts
    bestIndiv = open(best_indiv_name, 'w')
    
    #write out the chemical shift
    print closest_shifts
    for atom, chemicalShift in zip(params.Side_chain_carbons, closest_shifts['dummy']):
        bestIndiv.write('%s   %s\n'% (atom, chemicalShift))


    #calculate the distance
    distance = distFuncs.Distance(closest_shifts, Measured_CS,params)
    
    bestIndiv.write('Distance: %f' % (distance['dummy']))
    bestIndiv.close()

def write_smallest_distances(smallest_dist_in_generations, params):
    '''
    This function writes the smnallest distances to a file specified in 
    params.writeDisrance

    Parameters
    ========= 
    smallest_dist_in_generations : list 
        [generation , distance]
    params : Parameters class     
    '''

    f = open(params.writeDistance,'w')
    for i in smallest_dist_in_generations:
        f.write('%i  %0.5f\n'%(i[0], i[1]))
    f.close()

# def read_in_mean_results(res_file):
#     #gonna take my mean results form the results.dat file
#     #this should be in the same directory as the current script
#     res = open(res_file,'r')
#     file_checkpoint = 0
#     MeanResults = {}
#     for line in res.readlines():
#         split = line.split()
#         try:
#             if file_checkpoint == 1:
#                 MeanResults[split[0], split[1],split[2]] = split[3]
#             if split[0] == 'Mean':
#                 file_checkpoint = 1
#         except IndexError:
#             pass 
#     res.close()
#     return MeanResults
    
# def read_in_sec_struc_res(res_file):
#     ''' this function is here to read in the chemical shifts according to the secondry structure
#     where as before we had only a random coil '''

#     Results = {}
#     res = open(res_file, 'r')
    
#     #here we read in the results from the .dat file
#     check = 0
#     for line in res.readlines():
#         split = line.split()
        
#         try:
#             if split[0] == '#Mean':
#                 check = 1
#             #this check is so we ony read the first part of the data file
#             if check != 1:
#                 if line[0] != '#':
#                     Results[split[0], split[1], split[2],split[3],split[4]] = float(split[5])
#         except IndexError:
#             pass
#     bSheetData = {}
#     aHelixData = {}
#     #mean dictionaries for the secondry structures
#     MeanBSheetData = {}
#     MeanAHelixData = {}
    
#     for key in Results:

#         # these just define the redion of the aHelix in the Ramachadran plot
#         if float(key[1]) >= float(-140) and float(key[1]) <= float(-30):
#             if float(key[2]) >= float(-91) and float(key[2]) <= float(38):
#                 try:

#                     aHelixData[key[0],key[3],key[4]].append(Results[key])
#                 except KeyError:
#                     aHelixData[key[0],key[3],key[4]] = []
#                     aHelixData[key[0],key[3],key[4]].append(Results[key])
#                 #print 'added something to a list'
            
#         # these just define the redion of the bsheet in the Ramachadran plot
#         if float(key[1]) >= float(-170) and float(key[1]) <= float(-42):
#             if float(key[2]) >= float(65) and float(key[2]) <= float(180):
#                 try:
#                     bSheetData[key[0],key[3],key[4]].append(Results[key])
#                 except KeyError:
#                     bSheetData[key[0],key[3],key[4]] = []
#                     bSheetData[key[0],key[3],key[4]].append(Results[key])
#                 #print 'added something to b list'

#     #making the means results for 

#     for key in bSheetData:
#         meanBCs = np.divide(sum(bSheetData[key]),len(bSheetData[key]))
#         MeanBSheetData[key] = meanBCs
        
#         #these lines just gives the standard deviation should you need it 
#         #meanStd = np.std(bSheetData[key])
#         #MeanBSheetData[key].append(meanStd)

#     for key in bSheetData:
#         meanACs = np.divide(sum(aHelixData[key]),len(aHelixData[key]))
#         MeanAHelixData[key] = meanACs    
#     res.close()
    

#     return MeanAHelixData, MeanBSheetData

# def set_up(Mean_CS_file,  Measured_CS_file, sec_struct):
#     ''' Note that sec_struct can be set to A, B or R for alpha helix 
#         beta sheet of random coil / unknown '''
#     #this part checks that secondary  structure. -- could be done better because it reads Results.dat twice
#     R_coil_Mean_CS = read_in_mean_results(Mean_CS_file)
#     a_helix_CS, b_sheet_CS = read_in_sec_struc_res(Mean_CS_file)
#     Measured_CS_1 = ReadIn_measured_C_shifts(Measured_CS_file)
    
#     if sec_struct == 'R':
#         selected_Mean_CS = R_coil_Mean_CS
#     elif sec_struct == 'A':
#         selected_Mean_CS = a_helix_CS
#     elif sec_struct == 'B':
#         selected_Mean_CS = b_sheet_CS
#     else:
#         print '''ERROR: one of the options selected was not a recognised secondry structure type. 
#         Pease check the set_up() function'''
#         sys.exit()
#     return selected_Mean_CS, Measured_CS_1