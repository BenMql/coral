import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


energy = np.fromfile(path_2_files+'Timeseries/KE_total_full.dat', dtype = np.float_)
nelems = energy.shape[0]
print(str(nelems))

time = np.fromfile(path_2_files+'Timeseries/time.dat', dtype = np.float_, count = nelems)

plt.figure()
plt.semilogy(np.diff(time), 'o')
plt.title('Time-step size')

plt.figure()
plt.semilogy(time[1:], energy[1:], label='Total KE')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('KE')
plt.title('Total kinetic energy density')



plt.show()


