import numpy as np
import matplotlib.pyplot as plt
from sys import argv

plt.figure()
time = np.fromfile('./Timeseries/time.dat', dtype=np.float_)
for iarg in range(len(argv)-1):
    avgQty = np.fromfile('./Timeseries/'+argv[iarg+1] +'_volAvg.dat', dtype=np.float_)
    #plt.plot(time [:avgQty.shape[0]], avgQty, label=argv[iarg+1])
    plt.semilogy(time [:avgQty.shape[0]], avgQty, label=argv[iarg+1])

plt.legend()
plt.grid(True)
plt.grid(True, which='minor', linestyle='dotted')
plt.xlabel('Time')
plt.show()



