import numpy as np
import matplotlib.pyplot as plt
import read_TFP_input_file as rif

#___ Read the box dimension
Lx=rif.read_Lbox_X()
Lz=1.                     

#___ Read the resolution
resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)
NX = resolution[1]*3/2
NZ = resolution[3]*3/2

#___ Define grids accordingly
#... Equispaced in X
xgrid = np.linspace(0,Lx, num=NX, endpoint=False)
#... Gauss-Chebyshev in Z (inner-points)
zgrid = np.arange(NZ)
zgrid = np.cos((zgrid*2+1.0)*np.pi/(2*NZ))
zgrid += 1.0
zgrid *= 0.5

#___ Inquire what should be plotted
tstr = 't'+raw_input("Display what time? > ").zfill(8)
vstr = 'var'+raw_input('Display what var?  > ').zfill(3)
ystr = 'Y'+raw_input('Display what Y?    > ').zfill(4)

my_slice = np.fromfile('./XZ_Slices/'+vstr+'_'+ystr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NX+2,NZ)[:-2,:].T

slice_rms = np.sqrt(np.mean(my_slice**2))
print ('rms of the field sqrt('+vstr+'**2)'+str(slice_rms))
plt.figure(figsize=(6,6)) 
plt.pcolormesh(xgrid, zgrid, my_slice, cmap='RdYlBu_r')
plt.xlabel(r'azimuth $x$')
plt.ylabel(r'depth $z$')
plt.grid(True)
plt.colorbar()
plt.show()
