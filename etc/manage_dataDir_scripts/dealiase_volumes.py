import numpy as np
from sys import argv
from scipy.fftpack import fft2, ifft2, dct, idct

def read_and_deAlias(path_to_vols, vol_name):

   ierror = 0
   possible_padding=[0,2]

   domain_decomp_infos = np.fromfile(path_to_vols+
                        '../Geometry/domDecmp.core0000', dtype=np.int32)
   NX   = domain_decomp_infos[0]
   NY   = domain_decomp_infos[1]
   NZ   = domain_decomp_infos[2]
   NXAA = domain_decomp_infos[3]
   NYAA = domain_decomp_infos[4]
   NZAA = domain_decomp_infos[5]
   readVec = np.fromfile(path_to_vols + vol_name, dtype=np.float_)
   curPhys = readVec.reshape(NYAA, NXAA, NZAA)
   volPhys = np.copy(curPhys)
   curPhys =  dct (curPhys, axis=2, type=2)
   curPhys[:,:,0]*=0.5
   curSpec = np.zeros(curPhys.shape, dtype=np.complex_)
   curSpec.real = curPhys
   curSpec =  fft2(curSpec, axes=(0,1))
   aux1 = np.copy  (curSpec[:,:,:NZ])
   aux2 = np.delete(aux1, np.arange(NXAA/3+1,2*NXAA/3+1) , axis = 1)
   del aux1
   aux1 = np.delete(aux2, np.arange(NYAA/3+1,(2*NYAA/3+1)) , axis = 0)
   deaSpec =  np.copy(aux1)
   del aux2, aux1
   deaSpec = ifft2(deaSpec, axes=(0,1))
   deaPhys = np.copy(deaSpec.real)
   deaPhys[:,:,0]*=2.0
   deaPhys = idct(deaPhys, axis=2)*4./9./2./NZAA

   deaPhys.tofile(path_to_vols+vol_name[:-8]+'twoThirds.dat')

   return ierror

def clean_this_dir(my_dir):
   from os import walk
   for root, dirs, files in walk(my_dir+"/Volumes/"):
      counter = 0
      files_with_full = [x for x in sorted(files, key=str.lower) if 'full' in x]
      for my_file in files_with_full:
         counter += 1
         print ('[ '+str(counter).zfill(3)+'/'+str(len(files_with_full)).zfill(3) +' ] :: '+ root+my_file)

         if (not(my_file[:-8]+'twoThirds.dat' in files)):
            dumm_integer = read_and_deAlias(root, my_file)


directory_to_clean = argv[1]

clean_this_dir(directory_to_clean)

