import numpy as np
from sys import argv, exit
try:
   from scipy.fft import dct, idct
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
   xold = np.linspace(0,1,num=NXAA, endpoint=False)
   yold = np.linspace(0,1,num=NYAA, endpoint=False)
   xnew = np.linspace(0,1,num=ntup[0], endpoint=False)
   ynew = np.linspace(0,1,num=ntup[1], endpoint=False)
   aux1 = np.zeros((NXAA, ntup[1],ntup[2]), dtype=np.float_)
   for ix in range(NXAA):
     for iz in range(ntup[2]):
        aux1[ix,:,iz] = np.interp(ynew, yold, aux2[ix,:,iz], period=1.)
   deaPhys = np.zeros((ntup[0], ntup[1],ntup[2]), dtype=np.float_)
   for iy in range(ntup[1]):
     for iz in range(ntup[2]):
        deaPhys[:,iy,iz] = np.interp(xnew, xold, aux1[:,iy,iz], period=1.)
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

