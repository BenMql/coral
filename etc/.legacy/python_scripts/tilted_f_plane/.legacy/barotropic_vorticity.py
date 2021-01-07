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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vstr = 'var001'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

F = np.fromfile('./Volumes/Vol_'+vstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NY,NX,NZ)*iekman**(-1.0/3.0)


F_spec = fft.rfft2(F,axes=(0,1))
kx = fft.rfftfreq(NX)*2.*np.pi/Lx_box*iekman**(1.0/3.0)*NX
ky = fft.fftfreq (NY)*2.*np.pi/Ly_box*iekman**(1.0/3.0)*NY
KX, KY, KZ = np.meshgrid(kx, ky, grid_z)
F_tilted_spec = F_spec*np.exp(+1.j*KY*np.tan(latitude/180.*np.pi)*KZ)
F_tilted_phys = fft.irfft2(F_tilted_spec,axes=(0,1))

u_baro = np.sum(F_tilted_phys*dz,axis=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vstr = 'var002'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

F = np.fromfile('./Volumes/Vol_'+vstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NY,NX,NZ)*iekman**(-1.0/3.0)


F_spec = fft.rfft2(F,axes=(0,1))
kx = fft.rfftfreq(NX)*2.*np.pi/Lx_box*iekman**(1.0/3.0)*NX
ky = fft.fftfreq (NY)*2.*np.pi/Ly_box*iekman**(1.0/3.0)*NY
KX, KY, KZ = np.meshgrid(kx, ky, grid_z)
F_tilted_spec = F_spec*np.exp(+1.j*KY*np.tan(latitude/180.*np.pi)*KZ)
F_tilted_phys = fft.irfft2(F_tilted_spec,axes=(0,1))

v_baro = np.sum(F_tilted_phys*dz,axis=2)


# ########## FIG 1. ###############################
# ~~~~ plotting a YZ slice ~~~~
plt.figure(figsize=(6,6))
plt.pcolormesh(grid_y, grid_x,u_baro.T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title(r'Barotropic component, $u$')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.figure(figsize=(6,6))
plt.pcolormesh(grid_y, grid_x,v_baro.T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title(r'Barotropic component, $v$')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ---{ Now computing the barotropic vorticity }---
del KX, KY, KZ
KX, KY = np.meshgrid(kx,ky)
iK2 = KX**2 + KY**2
iK2 = 1./iK2
iK2[0,0] = 0.
u_fourier = fft.rfft2(u_baro)
v_fourier = fft.rfft2(v_baro)
zeta_fourier = 1j*KX*v_fourier - 1j*KY*u_fourier
PSI_fourier  = iK2*zeta_fourier
zeta_baro = fft.irfft2(zeta_fourier)
PSI_baro = fft.irfft2(PSI_fourier)
# ---{ Now plotting the barotropic vorticity }----
plt.figure(figsize=(6,6))
plt.pcolormesh(grid_y, grid_x,zeta_baro.T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title(r'Barotropic vorticity, $\zeta=\partial_xv-\partial_yu$')
# ---{ Now plotting the barotropic vorticity }----
plt.figure(figsize=(6,6))
plt.pcolormesh(grid_y, grid_x,PSI_baro.T)
plt.xlabel('latitude y')
plt.ylabel('azimuth x')
plt.title(r'Barotropic Stream-function, $\Psi=\nabla_\perp^{-2}\zeta$')



# ########## SHOW #################################
plt.show()


