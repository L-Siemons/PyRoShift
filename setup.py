from distutils.core import setup
from setuptools import setup, find_packages

descrip='''

     ====// cs2rotamer \\\\====

This package looks to determine protien 
side-chain rotamer populations from chemical 
shifts. 

Currently only isoleucine is implimented here.
This package is principally maintained by 
Lucas Siemons 

Author:
Lucas Siemons 

'''

setup(
    name='cs2rotamer',
    version='1.0',
    author='L. Siemons',
    author_email='zcbtla0@ucl.ac.uk',
    packages=find_packages(),
    #scripts=[''],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description=descrip,
    long_description=open('README.md').read(),
    #install_requires=['os','sys','random','numpy','time',],
    package_data={'cs2rotamer': ['resources/*']}
)
