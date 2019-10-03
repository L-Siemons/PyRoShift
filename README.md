
# PyRoShift

PyRoShift (PYthon ROtamers from chemical SHIFT) is a module that
calculates protein side-chain rotamer populations from carbon
chemical shifts. Currently it is implemented for isoleucine.
If you want to add another residue please get in touch!

If you use this program please cite the following:
XXX

# What does it do?

Principally this method takes the chemical shifts for
Ca, Cb, Cg1, Cg2 and Cd1 and determines the population
for each of the four rotamers:

- t/t
- m/m
- p/t
- m/t

Where:
- t = trans     (180 degrees)
- m = gauche -  (300 degrees)
- p = gauche +  ( 60 degrees)

Note that unlike many other methods this approach considers each rotamer to
be defined by both chi angles! The first letter denote the state of chi 1 and the
second the state of chi 2.

# Installation

On Linux and Mac one can install this module as with any other module.
1. Download this repository from the github page.
2. Go to the directory where the setup.py is
3. Install as follows:
```
$ python setup.py build
$ python setup.py install
```
or
```
pip install .
```
Note that depending on how your system is set up
you might need to use sudo.

# How

For an example on how to use this module please see the example_run/ directory.
Here there should be an example called run.py and also a ipython notebook describing all the
steps. To run the example script simply use
```
python run.py
```

# Authors

This module is written and maintained by
Lucas Siemons

If you have any question please feel free to
email me at lucas.siemons@googlemail.com.

On a more personal note this is my first Python
module and if you have any critical advice on how to
improve this body of code (both functionally and stylistically)
please let me know!
