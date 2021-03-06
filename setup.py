from setuptools import setup, find_packages

descrip='''
     ====// PyRoShift \\\\====

This package determines protein
side-chain rotamer populations from chemical
shifts.

Currently only isoleucine is implemented here.
If you wish to add more please get in touch!

If you use this work please cite our paper:

Determining isoleucine side-chain rotamer-sampling in 
proteins from 13C chemical shift.

Siemons L, Uluca-Yazgi B, Pritchard RB, McCarthy S, 
Heise H, Hansen F.  
ChemCom 2019; 

DOI:  10.1039/C9CC06496F

     ====// --------- \\\\====

Author:
Lucas Siemons

'''

setup(
    name='pyroshift',
    version='1.0',
    author='L. Siemons',
    author_email='lucas.siemons@googlemail.com',
    packages=find_packages(),
    license='LICENSE.txt',
    description=descrip,
    long_description=open('README.md').read(),
    install_requires=['numpy','scipy>=0.17.0', 'matplotlib'],
    package_data={'pyroshift': ['resources/*']}
)
