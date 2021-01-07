import numpy as np
import matplotlib.pyplot as plt
from os import walk
from os import remove
from time import sleep


def stitch_one_slice(p2data_prefix, what_slice, Nfull, nMPI):
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
   full_field = np.empty([0,Nfull], dtype = np.float_)
   for icore in range (nMPI):
      full_file_name = p2data_prefix+"core"+str(icore).zfill(4)+".dat"
      piece_of_field = np.fromfile(full_file_name, dtype=np.float_)
      if (what_slice == 'XY'):
         local_NY = piece_of_field.size/(Nfull+2)
         full_field = np.concatenate((full_field,np.reshape(piece_of_field,[local_NY,Nfull+2])[:,:-2] ), axis = 0)
      elif (what_slice == 'YZ'):
         local_NY = piece_of_field.size/Nfull
         full_field = np.concatenate((full_field,np.reshape(piece_of_field,[local_NY,Nfull]) ), axis = 0)
      else:
         print ("BAD what_slice ARGUMENT. CRASHING.")
         return
      remove(full_file_name)
   full_field.tofile(p2data_prefix+"full.dat")





#======================================================================================
#                                    . . . 
#======================================================================================

resolution = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]

slice_kind = 'XY'
N = resolution[1] * 3 / 2

path2data = "./"+slice_kind+"_Slices/"
full_lof = []

for (dirpath, dirnames, filenames) in walk(path2data):
   full_lof.extend(filenames)
#filter out files that have already been stitched:
lof = [s for s in full_lof if not('full' in s)]
del full_lof
#sort by alphabetical order
lof.sort(key=str.lower)

if (len(lof)==0):
   print('The list of file is empty. Skipping these slices.')
else:
   if (len(lof)==ncores):
      print('I have found exactly ncores = '+str(ncores)+' files')
      print ('Number of cores seemingly correct. Now proceeding to stitch slices')
      continue_bool=True
   else:
      print ('I have read ncores = '+str(ncores)+'.')
      print ('Here is the file ncores+1 files I found:')
      print(lof[ncores])
      if (lof[ncores][-12:] =='core0000.dat'):
         print ('Number of cores seemingly correct. Now proceeding to stitch slices')
         continue_bool=True
      else:
         print ('Something wrong with the number of cores. Stopping.')
         continue_bool=False
   
   if (continue_bool):
      # before starting, wait a bit, just to be sure that a file
      # is not being written by Fortran right now
      sleep(1)
      while lof!=[]:
         fname = lof[0]
         prefix = fname[:fname.find('core')]
         print(prefix)
         stitch_one_slice(path2data+prefix, slice_kind, N, ncores)
         lof = lof[ncores:]
                  

slice_kind = 'YZ'
N = resolution[3] * 3 / 2

path2data = "./"+slice_kind+"_Slices/"
full_lof = []

for (dirpath, dirnames, filenames) in walk(path2data):
   full_lof.extend(filenames)
#filter out files that have already been stitched:
lof = [s for s in full_lof if not('full' in s)]
del full_lof
#sort by alphabetical order
lof.sort(key=str.lower)

if (len(lof)==0):
   print('The list of file is empty. Skipping these slices.')
else:
   if (len(lof)==ncores):
      print('I have found exactly ncores='+str(ncores)+' files')
      print ('Number of cores seemingly correct. Now proceeding to stitch slices')
      continue_bool=True
   else:
      if (lof[ncores][-12:] =='core0000.dat'):
         print ('Number of cores seemingly correct. Now proceeding to stitch slices')
         continue_bool=True
      else:
         print ('Something wrong with the number of cores. Stopping.')
         continue_bool=False
   
   if (continue_bool):
      # before starting, wait a bit, just to be sure that a file
      # is not being written by Fortran right now
      sleep(1)
      while lof!=[]:
         fname = lof[0]
         prefix = fname[:fname.find('core')]
         print(prefix)
         stitch_one_slice(path2data+prefix, slice_kind, N, ncores)
         lof = lof[ncores:]
                              
