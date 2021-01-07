import numpy as np
import matplotlib.pyplot as plt

resolutions = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype = np.int32)
ncores = resolutions[0]
NX = resolutions[1]*3/2
NY = resolutions[2]*3/2
NZ = resolutions[3]*3/2

grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5

tstr = 't'+raw_input("What time should I display? >").zfill(8)
ekman =        input("Ekman number?               >")

t_profile = np.zeros(NZ,dtype=np.float_)
for icore in range(ncores):
   t_profile+= np.fromfile('./Z_Profiles/var50_Xmean_Ymean_'+tstr+'_core'+str(icore).zfill(4)+'.dat', dtype=np.float_)

t_profile /= NX
t_profile /= NY
t_profile *= ekman**(1.0/3.0)
print ('~~')
print ('The temperature fluctuations are bounded between '+str(t_profile.min())+' and '+str(t_profile.max())+'.')
plt.figure()
plt.plot(grid_z, t_profile)
plt.grid(True)
t_profile -= grid_z

plt.figure()
plt.plot(grid_z, t_profile)
plt.grid(True)
plt.show()

