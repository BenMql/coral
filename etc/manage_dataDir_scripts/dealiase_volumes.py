import numpy as np
from sys import argv
from scipy.fftpack import fft2, ifft2, dct, idct

def read_and_deAlias(path_to_vols, vol_name):

   ierror = 0
   possible_padding=[0,2]

   NX = np.fromfile(path_to_vols+'../CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[1]
   NY = np.fromfile(path_to_vols+'../CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[2]
   NZ = np.fromfile(path_to_vols+'../CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[3]
   NXAA = 3*NX/2
   NYAA = 3*NY/2
   NZAA = 3*NZ/2
   #print ('read NX,   NY,   NZ=   '+str(NX)+', '+str(NY)+', '+str(NZ))
   #print ('read NXAA, NYAA, NZAA= '+str(NXAA)+', '+str(NYAA)+', '+str(NZAA))
   #print ('inferred size NXAA*NYAA*NZAA= '+str(NXAA*NYAA*NZAA))
   readVec = np.fromfile(path_to_vols + vol_name, dtype=np.float_)
   #print ('length of read data: '+str(readVec.size))
   for pad in possible_padding:
      if (readVec.size == (NXAA+pad)*NYAA*NZAA):
         curPhys = readVec.reshape(NYAA, NXAA+pad, NZAA)[:,:-pad,:]
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

