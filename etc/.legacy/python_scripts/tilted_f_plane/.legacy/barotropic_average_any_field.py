import numpy as np
import numpy.fft as fft
import read_TFP_input_file as rif
import matplotlib.pyplot as plt

# =============================================================================================
# =============================================================================================
# ~~~~ Read the .sav binary file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)
NX = resolution[1]*3/2
NY = resolution[2]*3/2
NZ = resolution[3]*3/2

# ~~~~ Read the .in text file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iekman   = rif.read_iekman()
latitude = rif.read_latitude()
Lx_box   = rif.read_Lbox_X()
Ly_box   = rif.read_Lbox_Y()
# =============================================================================================
# =============================================================================================




tstr = 't'+raw_input("Display what time? > ").zfill(8)
vstr = 'var'+raw_input("Display what var.? > ").zfill(3)



grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5
grid_x = np.arange(NX)*Lx_box/NX/iekman**(1./3.)
grid_y = np.arange(NY)*Lx_box/NY/iekman**(1./3.)

grid_z_midpoints = grid_z[1:]+grid_z[:-1]
grid_z_midpoints*= 0.5
grid_z_midpoints = np.concatenate(([1],grid_z_midpoints,[0]))
dz = -np.diff(grid_z_midpoints)

GX, GY, GZ = np.meshgrid(grid_x, grid_y, grid_z)

T = np.fromfile('./Volumes/Vol_'+vstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NY,NX,NZ)*iekman**(-1.0/3.0)


T_spec = fft.rfft2(T,axes=(0,1))
kx = fft.rfftfreq(NX)*2.*np.pi/Lx_box*iekman**(1.0/3.0)*NX
ky = fft.fftfreq (NY)*2.*np.pi/Ly_box*iekman**(1.0/3.0)*NY
KX, KY, KZ = np.meshgrid(kx, ky, grid_z)
T_tilted_spec = T_spec*np.exp(+1.j*KY*np.tan(latitude/180.*np.pi)*KZ)
T_tilted_phys = fft.irfft2(T_tilted_spec,axes=(0,1))





# ########## FIG 1. ###############################
# ~~~~ plotting a YZ slice ~~~~
plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
plt.pcolormesh(grid_y, grid_z,T[:,0,:].T)
plt.plot(grid_y, grid_y,'--', color='k')
plt.plot(grid_y, 0.3+grid_y,'--', color='k')
plt.plot(grid_y, 0.6+grid_y,'--', color='k')
plt.ylim(0,1)
plt.xlabel('latitude y')
plt.ylabel('depth z')
plt.title('YZ slice')
plt.subplot(2,2,2)
plt.pcolormesh(grid_y, grid_z,T_tilted_phys[:,0,:].T)
plt.xlabel('latitude y')
plt.ylabel('depth z')
plt.title('YZ slice, un-tilted')
# ~~~~ Checking the XY slice (should be a mere translation) ~~~~
plt.subplot(2,2,3)
plt.pcolormesh(grid_y, grid_x,T[:,:,NZ/2].T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title('XY slice')
plt.subplot(2,2,4)
plt.pcolormesh(grid_y, grid_x,T_tilted_phys[:,:,NZ/2].T)
plt.xlabel('latitude y')
plt.ylabel(r'azimuth $x$')
plt.title('XY slice')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ########## FIG 2. ###############################
# ~~~~ plotting a YZ slice ~~~~
plt.figure(figsize=(6,6))
plt.pcolormesh(grid_y, grid_x,np.sum(T_tilted_phys*dz, axis=2).T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title('Barotropic component')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ########## SHOW #################################
plt.show()


