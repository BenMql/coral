import numpy as np
import matplotlib.pyplot as plt
from os import walk
from os import remove
from time import sleep


def stitch_one_profile(p2data_prefix, what_profile, Nfull, nMPI):
   """
   stitch_one_slice gathers the strips of slices output separately 
   by the different MPI processes, stitches them together correctly, and 
   saves the slice in a unique file.
      - p2data_prefix is a string that contains the path to the data,
        and the common prefix of each file, i.e., the strip saves by
        process 1234 is found in p2data_prefix+"core_1234.dat"
      - nMPI is the number of MPI processes
      - what_slice is a len=2 string: 'XY' or 'YZ' (XZ slices do not need stitching)
      - Nfull, is the "length" of the strip along the un-distributed dimension.
        For a 'XY' slice, this is NXAA.
        For a 'YZ' slice, this is NZAA.
   """
   full_profile = np.zeros([Nfull], dtype = np.float_)
   for icore in range (nMPI):
      full_file_name = p2data_prefix+"core"+str(icore).zfill(4)+".dat"
      full_profile += np.fromfile(full_file_name, dtype=np.float_)
      remove(full_file_name)
   full_profile.tofile(p2data_prefix+"full.dat")





#======================================================================================
#                                    . . . 
#======================================================================================

resolution = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]

profile_kind = 'Z'
N = resolution[3] * 3 / 2

path2data = "./"+profile_kind+"_Profiles/"
full_lof = []

for (dirpath, dirnames, filenames) in walk(path2data):
   full_lof.extend(filenames)
#filter out files that have already been stitched:
lof = [s for s in full_lof if not('full' in s)]
del full_lof
#sort by alphabetical order
lof.sort(key=str.lower)

if len(lof)==ncores:
  print ('I have found exactly ncores ='+str(ncores)+'files.')
else:
  print ('I have read ncores ='+str(ncores)+'.')
  print ('Here is the file ncores+1 files I found:')
  print ('..                                      ')
  print ('..                                      ')
  print(lof[ncores])
  print ('..                                      ')
  print ('..                                      ')
  print ('[file above should end with *core0000.dat]')
continue_bool = raw_input('Continue with ncores ='+str(ncores)+'? (yes/no)')

if (continue_bool=='yes'):
   # before starting, wait a bit, just to be sure that a file
   # is not being written by Fortran right now
   sleep(1)
   while lof!=[]:
      fname = lof[0]
      prefix = fname[:fname.find('core')]
      print(prefix)
      stitch_one_profile(path2data+prefix, profile_kind, N, ncores)
      lof = lof[ncores:]
else:
   print('You did not answer ''yes''. I am stopping now to avoid deleting anything')
               

               

