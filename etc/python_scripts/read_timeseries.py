import numpy as np
import matplotlib.pyplot as plt

# =========================================================================
# Read timeseries
#..........................................................................
time = np.fromfile('Timeseries/time.dat',  dtype=np.float_)
uu = np.fromfile('Timeseries/uu_full.dat', dtype=np.float_)
vv = np.fromfile('Timeseries/vv_full.dat', dtype=np.float_)
ww = np.fromfile('Timeseries/ww_full.dat', dtype=np.float_)

# =========================================================================
# check that all data has same size 
#..........................................................................
dlength = min (uu.shape[0], 
               vv.shape[0], 
               ww.shape[0], 
               time.shape[0])
uu = uu[:dlength]
vv = vv[:dlength]
ww = ww[:dlength]
time = time[:dlength]

# =========================================================================
# plot                                             
#..........................................................................
plt.semilogy(time, uu, 'o', label = r'$<u^2>$')
plt.semilogy(time, vv, 'o', label = r'$<v^2>$')
plt.semilogy(time, ww, 'o', label = r'$<w^2>$')
plt.semilogy(time, uu+vv+ww,'o',  label = r'$E_K$')
plt.grid(True)
plt.grid(True, which='minor', linestyle='dotted')
plt.legend()
plt.xlabel(r'time')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.show()


