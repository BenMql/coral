import numpy as np
import matplotlib.pyplot as plt
from sys import argv

plt.figure()
time = np.fromfile('./Timeseries/time.dat', dtype=np.float_)
for iarg in range(len(argv)-1):
    avgQty = np.fromfile('./Timeseries/'+argv[iarg+1] +'.bin', dtype=np.float_)
    plt.plot(time, avgQty [: time.shape[0]])
    del avgQty

plt.show()



