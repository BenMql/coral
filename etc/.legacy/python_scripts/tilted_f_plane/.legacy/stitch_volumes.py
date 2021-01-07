import numpy as np
import matplotlib.pyplot as plt
from os import walk
from os import remove
from time import sleep


def stitch_one_volume(p2data_prefix, NX, NZ, nMPI):
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
   #full_field = np.empty([0,NX, NZ], dtype = np.float_)
   list_of_fields = []
   list_of_file_names = []
   for icore in range (nMPI):
      full_file_name= p2data_prefix+"core"+str(icore).zfill(4)+".dat"
      piece_of_field = np.fromfile(full_file_name, dtype=np.float_)
      local_NY = piece_of_field.size/((NX+2)*NZ)
      list_of_fields.append(np.reshape(piece_of_field,[local_NY,NX+2,NZ])[:,:-2,:])
      list_of_file_names.append(full_file_name)
   full_field = np.concatenate(tuple(list_of_fields), axis = 0)
   full_field.tofile(p2data_prefix+"full.dat")
   for icore in range (nMPI):
      remove(list_of_file_names[icore])





#======================================================================================
#                                    . . . 
#======================================================================================



resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)

ncores = resolution[0]
NX = resolution[1]*3/2
NY = resolution[2]*3/2
NZ = resolution[3]*3/2

path2data = "./Volumes/"

# an empty list of files (l.o.f.) that contains everything
full_lof = []

for (dirpath, dirnames, filenames) in walk(path2data):
   full_lof.extend(filenames)

#filter out files that have already been stitched:
temp_lof = [s for s in full_lof if not('full' in s)]
lof      = [s for s in temp_lof if not('quart' in s)]
del full_lof, temp_lof

#sort by alphabetical order
lof.sort(key=str.lower)

print ('I have read ncores ='+str(ncores)+'.')
print ('The list of file has length :'+str(len(lof)))
if (len(lof)==0):
   print('The list of file is empty. Skipping these slices.')
   continue_bool = False
elif (len(lof)==ncores):
  print ('I have found exactly '+str(ncores)+' files.')
  continue_bool = True
else:
   print ('Here is the file ncores+1:')
   print ('~~                        ')
   print(lof[ncores])
   if (lof[ncores][-12:] =='core0000.dat'):
      print ('Number of cores seemingly correct. Now proceeding to stitch slices')
      continue_bool=True
   else:
      print ('Something wrong with the number of cores. Stopping.')
      continue_bool=False
   
print ('~~                        ')

if (continue_bool):

   # before starting, wait a bit, just to be sure that a file
   # is not being written by Fortran right now

   sleep(1)

   while lof!=[]:
      fname = lof[0]
      prefix = fname[:fname.find('core')]
      print(prefix)
      stitch_one_volume(path2data+prefix, NX, NZ, ncores)
      lof = lof[ncores:]
       
               
