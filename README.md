

  ========== // ========== ========= ========== // ==========
  ========== // ========== PyRoShift ========== // ==========
  ========== // ========== ========= ========== // ==========


  What is it?
  -----------

  Pyroshift (PYthon ROtamers from Chemical SHIFT) program
  calculates protein side chain rotamer populations from carbon
  chemical shifts. This is implemented for the following residues:

  - isoleucine

   If you use this program please cite the following:
   XXX


  When to use
  -----------

  This method can be used when the backbone conformation
  is known and the 5 side-chain chemical shifts can be provided.
  So far this method provides the rotamer distributions for the following
  five states in isoleucine.

  t/t, t/p,m/m,p/t,m/t

  Note;
  t = trans     (180 degrees)
  m = gauche -  (300 degrees)
  p = gauche +  ( 60 degrees)

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

  1)  An example of a working directory is shown in example_run/
        - Only two components are needed:
            1) The file with the chemical shifts
            2) the script that executes the functions in Pyroshift

  2) To Run the program one simply executes run.py
        - Note that one can edit the file names by changing the variables
            chemical_shifts_file
            output_file

  3) To use one's own data one can simply edit the file containing the chemical shifts
     the secondary structure only needs to be given in the row with the C alpha.
     The options for the secondary structure are alpha, beta and random coil (A,B and R)

  -----------------
  This program was written by Lucas Siemons.

  This work was carried out in collaboration with
  Dr Flemming Hansen