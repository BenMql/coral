import numpy as np
import matplotlib.pyplot as plt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


energy = np.fromfile(path_2_files+'Timeseries/KE_total_full.dat', dtype = np.float_)
nelems = energy.shape[0]
print(str(nelems))

time = np.fromfile(path_2_files+'Timeseries/time.dat', dtype = np.float_, count = nelems)

T_top = np.fromfile(path_2_files+'Timeseries/temperature_top', dtype = np.float_, count = nelems)
T_mid = np.fromfile(path_2_files+'Timeseries/temperature_mid_depth', dtype = np.float_, count = nelems)
T_bot = np.fromfile(path_2_files+'Timeseries/temperature_bottom', dtype = np.float_, count = nelems)

plt.figure()
plt.semilogy(np.diff(time), 'o')
plt.title('Time-step size')

plt.figure()
plt.semilogy(time[1:], energy[1:], label='Total KE')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('KE')
plt.title('Total kinetic energy density')


plt.figure()
plt.plot(time[1:], T_top[1:], label='Temperature top')
plt.plot(time[1:], T_mid[1:], label='Temperature mid')
plt.plot(time[1:], T_bot[1:], label='Temperature bot')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('Temperature')
plt.title('Horizontally averaged temperatures')
plt.legend(loc=2)


plt.show()


