import numpy as np
import matplotlib.pyplot as plt

time = np.fromfile('./Timeseries/time.dat', dtype=np.float_)
u2   = np.fromfile('./Timeseries/var001_sum_full.dat', dtype=np.float_)
v2   = np.fromfile('./Timeseries/var002_sum_full.dat', dtype=np.float_)
w2   = np.fromfile('./Timeseries/var003_sum_full.dat', dtype=np.float_)
p2   = np.fromfile('./Timeseries/var004_sum_full.dat', dtype=np.float_)
t2   = np.fromfile('./Timeseries/var005_sum_full.dat', dtype=np.float_)
s2   = np.fromfile('./Timeseries/var006_sum_full.dat', dtype=np.float_)
plt.semilogy(time, u2)
plt.semilogy(time, v2)
plt.semilogy(time, w2)
plt.semilogy(time, u2+v2+w2)
plt.figure()
plt.semilogy(time, p2)
plt.figure()
plt.semilogy(time, t2)
plt.figure()
plt.semilogy(time, s2)
plt.show()


