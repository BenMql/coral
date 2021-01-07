import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'
resolution = np.fromfile(path_2_files+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


energy = np.fromfile(path_2_files+'Timeseries/KE_total_full.dat', dtype = np.float_)
nelems = energy.shape[0]
print(str(nelems))

time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)

plt.figure()
plt.semilogy(np.diff(time))
plt.title('Time-step size')

plt.figure()
plt.semilogy(time[1:], energy[1:], label='Total KE')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
energy = np.fromfile(path_2_files+'Timeseries/KE_barotropic_full.dat', dtype = np.float_)
nelems = energy.shape[0]
print(str(nelems))
time = np.fromfile(path_2_files+'Timeseries/time_core0000.dat', dtype = np.float_, count = nelems)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plt.semilogy(time[1:], energy[1:], label='barotropic KE')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('KE')
plt.legend(loc=2)
plt.title('Total and barotropic kinetic energy density')



plt.show()


