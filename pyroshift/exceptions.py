
class Secondry_Structure_Error(Exception):
   '''
   Raised when the Secondry structure is not correct
   '''
   def __init__(self, residue, value):

       messgae ='''
       ERROR: Secondary Structure Input is not parsed correctly. 
       Please make sure the value after the C alpha shift is 
       one of the following 

       - alpha : a, alpha, h, helix
       - beta : b, beta, sheet, strand
       - coil : c, coil, r, (assumed to be 50%% alpha and 50%% beta)
       - A number between 0 and 1. 1 = 100%% alpha helix, 0. = 100%% beta sheet

       The Value given for residue %s is %s ''' %(residue, value)

       Exception.__init__(self,messgae)