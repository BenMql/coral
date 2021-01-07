import numpy as np
import matplotlib.pyplot as plt
import read_TFP_input_file as rif

#___ Read the box dimension
Lx=rif.read_Lbox_X()
Ly=rif.read_Lbox_Y()
Lz=1.                     

#___ Read the resolution
resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)
NX = resolution[1]
NY = resolution[2]
NZ = resolution[3]
NXAA = resolution[1]*3/2
NYAA = resolution[2]*3/2
NZAA = resolution[3]*3/2

#___ Define grids accordingly
#... Equispaced in X and Y
xgridAA = np.linspace(0,Lx, num=NXAA, endpoint=False)
ygridAA = np.linspace(0,Ly, num=NYAA, endpoint=False)
xgrid = np.linspace(0,Lx, num=NX, endpoint=False)
ygrid = np.linspace(0,Ly, num=NY, endpoint=False)
#... Gauss-Chebyshev in Z (inner-points)
zgridAA = np.arange(NZAA)
zgridAA = np.cos((zgridAA*2+1.0)*np.pi/(2*NZAA))
zgridAA += 1.0
zgridAA *= 0.5
#... Gauss-Chebyshev in Z (inner-points)
zgrid = np.arange(NZ)
zgrid = np.cos((zgrid*2+1.0)*np.pi/(2*NZ))
zgrid += 1.0
zgrid *= 0.5


continue_bool = True

label1 = {'X': r'Y',
          'Y': r'X',
          'Z': r'Y'}             
label2 = {'X': r'depth z',
          'Y': r'depth z',
          'Z': r'X'}
grid1AA = {'X': ygridAA,
           'Y': xgridAA,
           'Z': ygridAA}
grid2AA = {'X': zgridAA,
           'Y': zgridAA,
           'Z': xgridAA}
grid1 = {'X': ygrid,
        'Y': xgrid,
        'Z': ygrid}
grid2 = {'X': zgrid,
         'Y': zgrid,
         'Z': xgrid}
twoThirds_flag = False

while continue_bool:
   #___ Inquire what should be plotted
   orthogonal_to = raw_input("Slice at a given {X, Y, or Z}? >")
   where_index = int(raw_input("Slice at index "+orthogonal_to+"= "))
   tstr = 't'+raw_input("Display what time? > ").zfill(8)
   vstr = 'var'+raw_input('Display what var?  > ').zfill(3)
   
   try:
      my_vol = np.fromfile('./Volumes/Vol_'+vstr+'_'+tstr+'_full.dat', 
                                     dtype=np.float_).reshape(NYAA,NXAA+2,NZAA)[:,:-2,:]
   except IOError:
      my_vol = np.fromfile('./Volumes/Vol_'+vstr+'_'+tstr+'_twoThirds.dat', 
                                     dtype=np.float_).reshape(NY,NX,NZ)
      twoThirds_flag = True
   if orthogonal_to=='X':
      my_slice = my_vol[:,where_index,:]
   elif orthogonal_to=='Y':
      my_slice = my_vol[where_index,:,:]
   elif orthogonal_to=='Z':
      my_slice = my_vol[:,:,where_index]
   plt.figure(figsize=(6,6)) 
   if (twoThirds_flag):
     plt.pcolormesh(grid1[orthogonal_to], grid2[orthogonal_to], my_slice.T, cmap='RdYlBu_r')
   else:
     plt.pcolormesh(grid1AA[orthogonal_to], grid2AA[orthogonal_to], my_slice.T, cmap='RdYlBu_r')
   plt.xlabel(label1[orthogonal_to])
   plt.ylabel(label2[orthogonal_to])
   plt.grid(True)
   plt.colorbar()
   plt.show()
