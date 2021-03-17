import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pa

ista = 0

time = np.fromfile('./test/Timeseries/time.dat',  dtype=np.float_)
uu   = np.fromfile('./test/Timeseries/uu_volAvg.dat',  dtype=np.float_)
vv   = np.fromfile('./test/Timeseries/vv_volAvg.dat',  dtype=np.float_)
ww   = np.fromfile('./test/Timeseries/ww_volAvg.dat',  dtype=np.float_)
plt.figure()
plt.plot(time,np.sqrt(uu+vv+ww), color='C2')
plt.ylim(0,250)
plt.xlabel('time')
plt.ylabel(r'$Rm$')
plt.title(r'Magnetic Reynolds number $Rm$')


time = np.fromfile('./test/Timeseries/time.dat', dtype=np.float_)
bxbx = np.fromfile('./test/Timeseries/bxbx_volAvg.dat', dtype=np.float_)
byby = np.fromfile('./test/Timeseries/byby_volAvg.dat', dtype=np.float_)
bzbz = np.fromfile('./test/Timeseries/bzbz_volAvg.dat', dtype=np.float_)
plt.figure()
plt.semilogy(time, 5.e-6*(bxbx+byby+bzbz), color='C2')
plt.xlabel('time')
plt.ylabel(r'$\Lambda$')
plt.title (r'Elsasser number $\Lambda$')

plt.show()


