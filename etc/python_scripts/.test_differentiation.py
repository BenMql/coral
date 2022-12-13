import numpy as np
import matplotlib.pyplot as plt
from process_data import plane_layer_volume

a = plane_layer_volume(list_var_str=[])

a._xdiff_FD(pos_in_list=0)
plt.figure()
plt.plot(a.xgrid, a.dat[-1][:,0,0])
arg =  (a.xgrid+2*a.ygrid[0]) * 2 * np.pi / 10. +3.
plt.plot(a.xgrid, 2*np.pi/5.*np.cos(arg)*np.sin(arg)*np.exp(-np.cos(arg)**2), 'o')
a.xdiff(pos_in_list=0)
plt.plot(a.xgrid, a.dat[-1][:,0,0], '+')
print ('Relative error on the x derivative: '+str(np.sqrt(
       np.mean((a.dat[-1][:,0,0] - 
              2*np.pi/5.*np.cos(arg)*np.sin(arg)*np.exp(-np.cos(arg)**2))**2) /
       np.mean((a.dat[-1][:,0,0] )**2)
    )))


a.xdiff(pos_in_list=-1)
a.ydiff(pos_in_list= 0)
a.ydiff(pos_in_list=-1)
a.ydiff(pos_in_list= 0)
a.dat[-1] = a.dat[-2] + a.dat[-4]
a.lut[-1] = 'Lap a'
dxphi = np.pi/5
dyphi = np.pi/5*2.

a.inverseHorizontalLaplacian(pos_in_list=-1)

plt.figure()
plt.plot(a.xgrid, a.dat[-1][:,0,0])
plt.plot(a.xgrid, a.dat[ 0][:,0,0])
plt.show()

print ('Relative error on the inverse laplacian: ' +
        str(np.sqrt(np.mean(
            (a.dat[-1][:,:,0] - a.dat[0][:,:,0])**2
            )/ np.mean( a.dat[0][:,:,0]**2
                )
            ))
        )
