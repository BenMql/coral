import numpy as np
import matplotlib.pyplot as plt
import read_CLR_input_file as rif
import math
import matplotlib as mpl
mpl.rcParams.update({'font.size':30})

pi = math.pi

Nx =rif.read_NX()/2*3
Nz =rif.read_NZ()/2*3
Lz = rif.read_depth()
Lx= rif.read_Lbox_X()

grid_z = np.arange(Nz, dtype=np.float_)*2+1
grid_z = np.cos(grid_z*np.pi/2./Nz)
grid_z*= 0.5
grid_z-= 0.5
grid_z*= Lz

grid_x = np.linspace(0.,Lx, num=Nx)

Ek =rif.read_Ekman()
Ro =rif.read_Rossby()

Velx_t  = np.zeros(Nz)
Vely_t  = np.zeros(Nz)
Velz_t  = np.zeros(Nz)

thet = math.atan( ( ( (4.*pi)**2*Ek)**(2.) + (4.*pi)**2*Ek + 2  \
                )/( ( (4.*pi)**2*Ek)**(2.) - (4.*pi)**2*Ek +    \
                2) )

for jj in range(Nz):
	Velx_t[jj] = 4.*pi*Ek**(0.5)*( ( ( (4*pi)**2*Ek)**2 +    \
		4)/( ( (4*pi)**2*Ek)**2 + 1 ) )**(0.5)*np.cos(   \
		grid_z[jj]/(2.*Ek)**(0.5) - thet)*np.exp(grid_z[ \
		jj]/(2.*Ek)**(0.5) )*Ro - 1./(1. + ((4*pi)**2    \
		*Ek)**(2.0) )*Ro*np.exp(4.*pi*grid_z[jj])
	Vely_t[jj] = 4.*pi*Ek**(0.5)*( ( ( (4*pi)**2*Ek)**2 +    \
		4)/( ( (4*pi)**2*Ek)**2 + 1 ) )**(0.5)*np.sin(   \
		grid_z[jj]/(2.*Ek)**(0.5) - thet)*np.exp(grid_z[ \
		jj]/(2.*Ek)**(0.5) )*Ro + (4.*pi)**2*Ek/(1. +    \
		((4*pi)**2*Ek)**(2.0) )*Ro*np.exp(4.*pi*grid_z[  \
		jj])
	Velz_t[jj] = 0.0

t_str = raw_input('Display what time? >').zfill(8)
#v_str = raw_input('Display what var.? >').zfill(3)
#t_str = str(3000).zfill(8)

ux = np.fromfile('./XZ_Slices/var001_t'+t_str+'_full.dat', dtype= \
		np.float_).reshape(Nz, Nx)
uy = np.fromfile('./XZ_Slices/var002_t'+t_str+'_full.dat', dtype= \
		np.float_).reshape(Nz, Nx)
uz = np.fromfile('./XZ_Slices/var003_t'+t_str+'_full.dat', dtype= \
		np.float_).reshape(Nz, Nx)

fig1=plt.figure()
axe1=fig1.add_subplot(111)

fig2=plt.figure()
axe2=fig2.add_subplot(111)

fig3=plt.figure()
axe3=fig3.add_subplot(111)

axe1.plot(np.squeeze(np.mean(ux, axis = 1)),grid_z,color='b', \
        label = 'V_x', linewidth = 3.0)
axe1.plot(Velx_t,grid_z,color='k',label = 'V_x th', linewidth = \
        3.0)

temp = np.squeeze(np.mean(uy, axis = 1))
#temp = temp + (Vely_t[0] - temp[0] )

#axe2.plot(np.squeeze(np.mean(uy, axis = 1)),grid_z,color='b', \
#        label = 'V_y', linewidth = 3.0)
axe2.plot(temp,grid_z,color='b', \
        label = 'V_y', linewidth = 3.0)
axe2.plot(Vely_t,grid_z,color='k',label = 'V_y th', linewidth = \
        3.0)

axe3.plot(np.squeeze(np.mean(uz, axis = 1)),grid_z,color='b', \
        label = 'V_y', linewidth = 3.0)
axe3.plot(Velz_t,grid_z,color='k',label = 'V_y th', linewidth = \
        3.0)

axe1.set_xlabel(r'$u_x$',fontsize=60)
axe1.set_ylabel(r'$z$',fontsize=60)
axe2.set_xlabel(r'$u_y$',fontsize=60)
axe2.set_ylabel(r'$z$',fontsize=60)
axe3.set_xlabel(r'$u_z$',fontsize=60)
axe3.set_ylabel(r'$z$',fontsize=60)

fig1.show()
fig2.show()
fig3.show()

plt.show()
