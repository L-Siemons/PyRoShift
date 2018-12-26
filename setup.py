from setuptools import setup, find_packages

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
    name='pyroshift',
    version='1.0',
    author='L. Siemons',
    author_email='zcbtla0@ucl.ac.uk',
    packages=find_packages(),
    license='LICENSE.txt',
    description=descrip,
    long_description=open('README.md').read(),
    install_requires=['os','sys','random','numpy','time','scipy'],
    package_data={'pyroshift': ['resources/*']}
)
