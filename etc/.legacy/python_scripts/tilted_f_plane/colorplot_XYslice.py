import numpy as np
import read_TFP_input_file as rif
import matplotlib.pyplot as plt
import numpy.fft as fft

#___ Read the box dimension
Lx=rif.read_Lbox_X()
Ly=rif.read_Lbox_X()

#___ Read the resolution
resolution = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
NX = resolution[1]*3/2
NY = resolution[2]*3/2

#___ Define grids accordingly
xgrid = np.linspace(0,Lx, num=NX, endpoint=False)
ygrid = np.linspace(0,Ly, num=NY, endpoint=False)

#___ Inquire what should be plotted
Zstr = 'Z'   + raw_input("What height?    > ").zfill(4)
tstr = 't'   + raw_input("What timestep?  > ").zfill(8)
vstr = 'var' + raw_input("What var.?      > ").zfill(3)

#___ Get the data                                      
path2file = './XY_Slices/'
my_slice = np.fromfile(path2file+vstr+'_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NY, NX+2)[:,:-2]

#___ Display some info                                 
print (vstr+" max is "+str(my_slice.max()))
print (vstr+" min is "+str(my_slice.min()))
print (vstr+" std is "+str(np.std(my_slice)))
print (vstr+" average is "+str(np.mean(my_slice)))

plt.figure()
plt.pcolormesh(xgrid, ygrid, my_slice)
plt.title(vstr)
plt.colorbar()
plt.show()


