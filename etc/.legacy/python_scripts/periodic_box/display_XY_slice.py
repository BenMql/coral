import numpy as np
import matplotlib.pyplot as plt
from sys import argv

icore = 0
domDecmp = np.fromfile('./Matrices/domDecmp.core'+str(icore).zfill(4), dtype=np.int32)

var_str = 'var'+argv[1].zfill(4)
posn_str =  argv[2].zfill(5)
time_to_plot = argv[3].zfill(8)
kind_str = 'z'


phys = np.fromfile('XY_Slices/'+var_str+'_'+kind_str+posn_str+'_t'+time_to_plot+'_full.dat', dtype = np.float_).reshape(domDecmp[4], domDecmp[3])
plt.pcolormesh(phys[:,:]); 
plt.colorbar()   
plt.show()



