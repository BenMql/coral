import numpy as np
import matplotlib.pyplot as plt
import read_TFP_input_file as rif

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'
resolution = np.fromfile(path_2_files+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]
iekman = rif.read_iekman()
ekman = 1./iekman


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nusselt = np.fromfile(path_2_files+'Timeseries/heat_flux_top_full.dat', dtype = np.float_)
nelems = nusselt.shape[0]
print(str(nelems))

nusselt-= 1.0
nusselt*=ekman**(1./3.)
nusselt+= 1.0

time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)

plt.figure()
plt.plot(time, nusselt, color='C0',label='top')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('Heat Flux')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nusselt = np.fromfile(path_2_files+'Timeseries/heat_flux_bot_full.dat', dtype = np.float_)
nelems = nusselt.shape[0]
print(str(nelems))

nusselt-= 1.0
nusselt*=ekman**(1./3.)
nusselt+= 1.0

time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)

plt.plot(time, nusselt, '--', color='C0',label='top')




plt.show()


