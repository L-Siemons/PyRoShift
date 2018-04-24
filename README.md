

  ========== // ========== ========= ========== // ========== 
  ========== // ========== PyRoShift ========== // ========== 
  ========== // ========== ========= ========== // ========== 



  What is it?
  -----------

  This program calculates protein side chain rotamer populations from carbon 
  chemical shifts. This is implemented for the following residues: 

  - isoleucine

   If you use this program please cite the following: 
   XXX


  When to use
  -----------

  Our bench mark shows that the calculation of rotamer populations 
  is most accurate when all the carbon chemical shifts are provided 
  as this allows an actuate calculation for the following states: 

  t/t, t/p,m/m,p/t,m/t

  Note; 
  t = trans 
  m = gauche - 
  p = gauche + 

  However we also noted that if only the C delta is provided then we are able 
  to obtain the population of the m/m state. If not all the chemical shifts 
  are provided then then calculating the rotamer populations for chi1 is not advised.


  Installation
  ------------

  On Linux and Mac one can install the module as with any other module using 
  
  $ python setup.py build
  $ python setup.py install

  or 

  pip install .
  

  How to use
  -----------
  
  To use this program: 
  
  1)  An example of a working directory is shown in exmapleRun/
      This folder contains: 
      	- InputFile.txt
      		- contains the parameters
      	- Measured_shifts.txt
      		- lists the measured chemical shifts
      	- RunCs2Rotamer.py
      		- a script that executes the method

  2) To Run the program one simply executes RunCs2Rotamer.py
     For future use this can be moved to somewhere in the executable path if desired. 

  3) To use one's own data one can simply edit the Measured_shifts.txt file and the
     configurations in InputFile.txt. Note the that configurations in InputFile.txt
     are the ones used in our bench mark and the ones we suggest using.


  Configurations
  -----------
  The File names Input.txt contains the configuration for the Genetic algorithm. 
  This section contains a description of the configuration options in this file and
  the default options which were used to bench mark the algorithm will be indicated by * . 
  Changing these options is likely to affect the performance of the algorithm and so should
  be done with ca
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