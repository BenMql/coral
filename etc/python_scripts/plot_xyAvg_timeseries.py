import numpy as np
import matplotlib.pyplot as plt
from sys import argv

plt.figure()
time = np.fromfile('./Timeseries/time.dat', dtype=np.float_)
for iarg in range(0,len(argv)-1,2):
    avgQty = np.fromfile('./Timeseries/'+argv[iarg+1]+'_XYavg_z'+str(int(argv[iarg+2])).zfill(5)+'.dat', dtype=np.float_)
    plt.plot(time, avgQty [: time.shape[0]])
    del avgQty

plt.show()



