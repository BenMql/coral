import numpy as np
import matplotlib.pyplot as plt
from sys import argv

icore = 0
domDecmp = np.fromfile('./Matrices/domDecmp.core'+str(icore).zfill(4), dtype=np.int32)

var_str = argv[1].zfill(3)
time_to_plot = argv[2].zfill(8)


phys = np.fromfile('Volumes/Vol_var'+var_str+'_t'+time_to_plot+'.phys', dtype = np.float_).reshape(domDecmp[5], domDecmp[4], domDecmp[3])
plt.pcolormesh(phys[1,:,:]); 
plt.figure()
plt.pcolormesh(phys[:,0,:]   )
plt.colorbar()   
plt.show()



