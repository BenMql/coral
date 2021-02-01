import numpy as np
from sys import argv, exit
try:
   from scipy.fft import dct, idct, rfft, irfft
except ImportError:
   import scipy as sp
   print ('The fft module was not found in scipy.')
   print ('This is possibly due to a scipy version older than 1.4.0')
   print ('For info, your scipy version is: '+sp.version.full_version)
   exit(1)
   


def read_and_resize_zeroModes(path_to_vols, vol_name, ntup, toDir='./'):
   '''
   Open and pad (or truncate) the Chebyshev coefficients of the mean fields
   stored in CheckPoints. In contrast to fluctations (kxky fileds), spectral 
   coefficients are stored directly. Hence, instead of NZAA (grid points),
   we store (NZ - nbc), where NZ = NZAA*3//2 and nbc is the number of 
   boundary conditions for a given field.
   '''
   ierror = 0
   domain_decomp_infos = np.fromfile(path_to_vols+
                        '../Geometry/domDecmp.core0000', dtype=np.int32)
   readVec = np.fromfile(path_to_vols + vol_name, dtype=np.float_)
   Nold= readVec.shape[0]
   Nnew = ntup[2]
   
   coscoefs = dct(readVec)
   if (Nnew>Nold):
      newcoefs = np.zeros((Nnew), dtype=np.float_)
      newcoefs[:Nold] = coscoefs
   else:
      newcoefs = coscoefs[:Nnew]

   aux2 = idct(newcoefs)*Nnew/Nold
   aux2.tofile(toDir+'/Restart/'+vol_name)
   return ierror

def read_and_resize(path_to_vols, vol_name, ntup, toDir='./'):
   '''
   Open and pad (or truncate) the Chebyshev-Fourier-Fourier coefficients 
   of the kxky fields stored in CheckPoints. As they are stored in physical
   space, the interpolation is done by cosine transform along z, followed by padding 
   or truncation, then inverse cosine transform. Along x and y, the grid is 
   evenly space, and we elect to use linear interpolation (which seems faster for
   large resolutions, according to a few tests).
   '''
   ierror = 0
   domain_decomp_infos = np.fromfile(path_to_vols+
                        '../Geometry/domDecmp.core0000', dtype=np.int32)
   NX   = ntup[0]
   NY   = ntup[1]
   NZ   = ntup[2]
   NXAA = domain_decomp_infos[3]
   NYAA = domain_decomp_infos[4]
   NZAA = domain_decomp_infos[5]
   print ('----------- :: Resizing (' + str(NXAA)+','+str(NYAA)+','+str(NZAA)+') into'+
                        ' (' + str(NX  )+','+str(NY  )+','+str(NZ  )+').')
   curPhys = np.fromfile(path_to_vols + vol_name, dtype=np.float_).reshape(NXAA, NYAA, NZAA)
   coscoefs = dct(curPhys, axis=-1)
   if (NZ>NZAA):
      newcoefs = np.zeros((NXAA,NYAA,NZ), dtype=np.float_)
      newcoefs[:,:,:NZAA] = coscoefs
   else:
      newcoefs = coscoefs[:,:,:NZ]

   aux2 = idct(newcoefs,axis=-1)*NZ/NZAA

   ## now we need to interpolate in the xy-plane...
   aux1 = np.zeros((NXAA, ntup[1],ntup[2]), dtype=np.float_)
   if (NYAA==ntup[1]):
     aux1=np.copy(aux2)
   elif (NYAA<ntup[1]):
     spectral_buffer = rfft(aux2, axis=1)
     padded_spectral_buffer = np.zeros((NXAA, ntup[1]//2+1, ntup[2]), dtype=np.complex_)
     padded_spectral_buffer[:,:(NYAA//2+1),:] = spectral_buffer
     aux1 = irfft(padded_spectral_buffer, axis=1)/NYAA*ntup[1]
     del padded_spectral_buffer, spectral_buffer
   else:
     spectral_buffer = rfft(aux2, axis=1)
     truncated_spectral_buffer = np.zeros((NXAA, ntup[1]//2+1, ntup[2]), dtype=np.complex_)
     truncated_spectral_buffer[:,:(ntup[1]//3),:] = spectral_buffer[:,:(ntup[1]//3),:]
     truncated_spectral_buffer.imag[:,-1,:]=0.
     truncated_spectral_buffer.imag[:, 0,:]=0.
     aux1 = irfft(truncated_spectral_buffer, axis=1)/NYAA*ntup[1]
     del truncated_spectral_buffer, spectral_buffer
   del aux2
   deaPhys = np.zeros((ntup[0], ntup[1],ntup[2]), dtype=np.float_)
   if (NXAA==ntup[0]):
     deaPhys=np.copy(aux1)
   elif (NXAA<ntup[0]):
     spectral_buffer = rfft(aux1, axis=0)
     padded_spectral_buffer = np.zeros((ntup[0]//2+1, ntup[1], ntup[2]), dtype=np.complex_)
     padded_spectral_buffer[:(NXAA//2+1),:,:] = spectral_buffer
     deaPhys = irfft(padded_spectral_buffer, axis=0)/NXAA*ntup[0]
     del padded_spectral_buffer, spectral_buffer
   else:
     spectral_buffer = rfft(aux1, axis=0)
     truncated_spectral_buffer = np.zeros((ntup[0]//2+1, ntup[1], ntup[2]), dtype=np.complex_)
     truncated_spectral_buffer[:(ntup[0]//3),:,:] = spectral_buffer[:(ntup[0]//3),:,:]
     truncated_spectral_buffer.imag[-1,:,:]=0.
     truncated_spectral_buffer.imag[ 0,:,:]=0.
     deaPhys = irfft(truncated_spectral_buffer, axis=0)/NXAA*ntup[0]
     del truncated_spectral_buffer, spectral_buffer
   del aux1
   deaPhys.tofile(toDir+'/Restart/'+vol_name)
   return ierror

def resize_checkpoints(fromDir='./', toDir = './', newResolution=(64,64,64)):
   from os import walk
   from time import time
   for root, dirs, files in walk(fromDir+"/CheckPoints/"):
      counter = 0
      files_with_kxky = [x for x in sorted(files, key=str.lower) if 'kxky' in x]
      t0 = time()
      for my_file in files_with_kxky:
         counter += 1
         print ('[ '+str(counter).zfill(3)+'/'+str(len(files_with_kxky)).zfill(3) +' ] :: '+ root+my_file)
         dumm_integer = read_and_resize(root, my_file, newResolution, toDir = toDir)
         print ('Elapsed time, interpolation: '+str(time()-t0)) 
         #dumm_integer = read_and_resize_fft(root, my_file, ntup)
         #print ('Elapsed time, fft          : '+str(time()-t0)) 
      files_with_zero = [x for x in sorted(files, key=str.lower) if 'zero' in x]
      for my_file in files_with_zero:
         counter += 1
         print ('[ '+str(counter).zfill(3)+'/'+str(len(files_with_zero)).zfill(3) +' ] :: '+ root+my_file)
         dumm_integer = read_and_resize_zeroModes(root, my_file, newResolution, toDir = toDir)



old_dir = argv[1]
new_dir = argv[2]

resize_checkpoints(fromDir = old_dir, toDir = new_dir, 
                   newResolution=(int(argv[3]), int(argv[4]), int(argv[5])))

