import numpy as np
import matplotlib.pyplot as plt

resolutions = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)
ncores = resolutions[0]
NX = resolutions[1]*3/2
NY = resolutions[2]*3/2
NZ = resolutions[3]*3/2
ekman = input("what is ekman? (float64) > ")
tstr = 't' + raw_input("What timestep do you want displayed? > ").zfill(8)


#~~~~~~~~~~~~~  the grid  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5

#~~~~~~~~~~~~~  Read the data  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t_profile = np.zeros(NZ,dtype=np.float_)
for icore in range(ncores):
   t_profile+= np.fromfile('./Z_Profiles/var50_Xmean_Ymean_'+tstr+'_core'+str(icore).zfill(4)+'.dat', dtype=np.float_)

#~~~~~~~~~~~~~  Rescale the data  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t_profile /= NX
t_profile /= NY
t_profile *= ekman**(1.0/3.0)
print ('~~')
print ('The temperature fluctuations are bounded between '+str(t_profile.min())+' and '+str(t_profile.max())+'.')
t_profile -= grid_z


#~~~~~~~~~~~~~  dT/dz  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dT = np.diff(t_profile)
dz = np.diff(grid_z)
Tz = dT/dz

#~~~~~~~~~~~~~  Nusselts  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print ("nusselt top    is:"+ str(-Tz[-1]))
print ("nusselt bottom is:"+ str(-Tz[ 1]))

#~~~~~~~~~~~~~  Mid-depth gradients  ~~~~~~~~~~~~~~~~~~~~~~~~~
Tz_mean = np.mean(Tz[3*NZ/8:5*NZ/8])

d2T= np.diff(dT)
dz2= np.diff(dz)
T2z = d2T/dz2

mid_grid_1 = 0.5*(grid_z    [:-1]+grid_z    [1:])
mid_grid_2 = 0.5*(mid_grid_1[:-1]+mid_grid_1[1:])


#~~~~~~~~~~~~~  Figures  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt.figure()
plt.plot(grid_z, t_profile)
plt.plot(grid_z[3*NZ/8:5*NZ/8], -0.5*(Tz_mean+1)+grid_z[3*NZ/8:5*NZ/8]*Tz_mean)
plt.title('Temperature profile')
plt.xlabel('depth z')
plt.ylabel('Temperature')
plt.grid(True)

plt.figure()
plt.plot(mid_grid_1, Tz)
plt.plot(grid_z, grid_z*0+Tz_mean,'--', color='k')
plt.title('Temperature profile')
plt.xlabel('depth z')
plt.ylabel('Temperature')
plt.grid(True)
#~~~~~~~~~~~~~  Show  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.show()
