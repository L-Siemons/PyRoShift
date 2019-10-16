
PyRoShift
=========

PyRoShift (PYthon ROtamers from chemical SHIFT) is a module that
calculates protein side-chain rotamer populations from carbon
chemical shifts. Currently it is implemented for isoleucine.
If you want to add another residue please get in touch!

If you use this program please cite the following:

Determining isoleucine side-chain rotamer-sampling in proteins from 13C chemical shift.
Siemons L, Uluca-Yazgi B, Pritchard RB, McCarthy S, Heise H, Hansen F.
ChemCom 2019;

DOI:  10.1039/C9CC06496F

Site: https://pubs.rsc.org/en/Content/ArticleLanding/2019/CC/C9CC06496F#!divAbstract

What does it do?
----------------


Principally this method takes the isoleucine 13C chemical shifts
(Ca, Cb, Cg1, Cg2 and Cd1) and determines the population
for each of the four rotamers:

- t/t
- m/m
- p/t
- m/t

Where:
- t = trans     (180 degrees)
- m = gauche -  (300 degrees)
- p = gauche +  ( 60 degrees)

Unlike many other methods this approach considers each rotamer to
be defined by both chi angles! The first letter denotes the state of chi 1 and the
second the state of chi 2.


Installation
------------


On Linux and Mac one can install this module as with any other module. To do this follow these steps:

1. Download this repository from the github page. (https://github.com/L-Siemons/PyRoShift)
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


How to use?
-----------

For an example on how to use this module please see the example_run/ directory.
Here there should be an example called run.py and also a ipython notebook describing all the
steps. Sometimes the notebook doesn't render properly on github so I have also exported it as
an html file that can be opened in a browser.

To run the example script simply use
```
python run.py
```


Authors
-------

This module is written and maintained by
Lucas Siemons

If you have any question please feel free to
email me at lucas.siemons@googlemail.com.

On a more personal note this is my first Python
module and if you have any critical advice on how to
improve this body of code (both functionally and stylistically)
please let me know!
