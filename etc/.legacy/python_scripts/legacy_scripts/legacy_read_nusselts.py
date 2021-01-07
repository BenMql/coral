import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'
resolution = np.fromfile(path_2_files+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]
ekman = input("what is ekman? (float64) > ")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = np.fromfile(path_2_files+'Timeseries/heat_flux_top_core0000.dat', dtype = np.float_)
nusselt = np.copy(a)


nelems = nusselt.shape[0]
print(str(nelems))

for icore in range(ncores-1):
  nusselt += np.fromfile(path_2_files+'Timeseries/heat_flux_top_core'+str(icore+1).zfill(4)+'.dat', dtype = np.float_, count = nelems)

nusselt /= ncores
nusselt-= 1.0
nusselt*=ekman**(1./3.)
nusselt+= 1.0

time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)

plt.figure()
plt.plot(time, nusselt, color='C0',label='top')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('Heat Flux')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = np.fromfile(path_2_files+'Timeseries/heat_flux_bot_core0000.dat', dtype = np.float_)
nusselt = np.copy(a)


nelems = nusselt.shape[0]
print(str(nelems))

for icore in range(ncores-1):
  nusselt += np.fromfile(path_2_files+'Timeseries/heat_flux_bot_core'+str(icore+1).zfill(4)+'.dat', dtype = np.float_, count = nelems)

nusselt /= ncores
nusselt-= 1.0
nusselt*=ekman**(1./3.)
nusselt+= 1.0

time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)

plt.plot(time, nusselt, '--', color='C0',label='top')




plt.show()


