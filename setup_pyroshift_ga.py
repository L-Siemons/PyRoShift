from distutils.core import setup
from setuptools import setup, find_packages
from Cython.Build import cythonize

descrip='''

     ====// cs2rotamer \\\\====

This package looks to determine protein
side-chain rotamer populations from chemical
shifts.

Currently only isoleucine is implemented here.
This package is principally maintained by
Lucas Siemons

Author:
Lucas Siemons

'''

setup(
    name='pyroshiftGA',
    version='1.0',
    author='L. Siemons',
    author_email='zcbtla0@ucl.ac.uk',
    packages=find_packages(),
    #scripts=[''],
    #url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description=descrip,
    long_description=open('README.md').read(),
    install_requires=['cython','numpy',],
    package_data={'cs2rotamer': ['resources/*']},
    ext_modules=cythonize('pyroshiftGA/*.pyx')
)
