import numpy as np
import matplotlib.pyplot as plt

resolution = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype=np.int32)
ncores = resolution[0]
path_2_files = './'
a = np.fromfile(path_2_files+'kinetic_energy_core0000.dat', dtype = np.float_)
energy = np.copy(a)

nelems = energy.shape[0]
print(str(nelems))

for icore in range(ncores-1):
  energy += np.fromfile(path_2_files+'kinetic_energy_core'+str(icore+1).zfill(4)+'.dat', dtype = np.float_, count = nelems)

energy /= ncores 

time = np.fromfile(path_2_files+'time_core0000.dat', dtype = np.float_, count = nelems)

plt.figure()
plt.semilogy(time, energy)
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('Energy')



plt.figure()
plt.semilogy(np.diff(time))
plt.title('Time-step size')

plt.show()


