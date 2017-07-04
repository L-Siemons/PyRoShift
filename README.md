                          Shifts2PopsPy


  What is it?
  -----------

  This program calculates rotamer populations from carbon 
  chemical shifts. This is implimented for the following residues: 

  - isoleucine

   If you use this program please cite the following: 
   XXX


  When to use
  -----------

  Our bench mark shows that the calculation of rotamer populations 
  is most accurate when all the carbon chemical shifts are provided 
  as this allows an acutate calculation for the following states: 

  t/t, t/p,m/m,p/t,m/t

  However we also noted that if the C delta is provided then we are able 
  to sum these states to obtain the rotamer populations for chi 2 
  independant of chi1. 

  Note that if the following

  C alpha 
  C beta 
  C gamma1 
  C gamma2 

  are not provided then the rotamer populations for chi1 are unreliable. 


  How to use
  -----------
  
  To use this program: 
  
  1)  Place the following into the working directory: 
          - Input.txt 
          - Master_2.py 
          - File containing the theoretical chemical shifts. 

  2)  Modify the file containing the experimentally observed 
      chemical shifts. This File is called 'Measured_shifts.txt'
      by defult. 

  3)  Configure the values in InputFile.txt. The defult values are
      the ones that were used during the bench mark. We suggest using these. 
      However if one sees that the distances between the best individual 
      and the experimental chemical shifts then it would be appropriate to 
      terminate the algorithm.

  4)  To run the algorithm simply execute the python script Master_2.py 
      in Mac/Linux enviroments this can be done as follows: 

      ./Master_2.py
      python Master_2.py



  Configurations
  -----------
  
  The File names Input.txt contians the configuration for the Genetic algorithm. 
  This section cantains a description of the configuration options in this file and
  the defult options which were used to bench mark the algorithm will be indictaed by *. 
  Changing these options is likely to affect the performance of the algorithm and so should
  be done with care. 

  Settings: 

  dist_func*             Euclid_dis      Distance function that drives the GA
  pop_size*              200             Total number of individuals in the generation
                                             - can be any integer above 0 

  Mutation_rate*         60              Mutation rate in the GA
                                             - Range 0-100

  verbose                1               0 write everything, 1 write less 
  
  Max_gen_number*        150             total number of generations 
                                             - can be any integer above 0
  
  residue                ILE             Which residue are we trying to predict the rotamer distribution
                                             - This takes the residue 3 letter code
  
  File IO:

  CalcRes                Results_DFTminima.dat       File with the calculated results
  
  ExpShifts              Measured_shifts.txt         File with the measured chemical shifts
  
  sol_name               Solutions.txt               file contaning the solution
  
  write_out_location     test/                       Where the results from the GA are written out 



  Authors
  -----------------

  This program was written by Lucas Siemons.

  This work was carried out in collaboration with 
  Dr Flemming Hansen 
  Ruth Dingle