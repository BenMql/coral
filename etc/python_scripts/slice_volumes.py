import numpy as np
import matplotlib.pyplot as plt

time = 10000
var = 4 # this is temperature
#------------------------------------------------
NX = 256 # should match coral.parameters.in entry
NY = 512 # should match coral.parameters.in entry
NZ = 384 # should match coral.parameters.in entry
#------------------------------------------------

NXAA = int(NX /2 * 3)
NYAA = int(NY /2 * 3)
NZAA = int(NZ /2 * 3)


z = 0.5+ 0.5*np.cos ((
       np.linspace(0,NZAA, num=NZAA, endpoint=False)*2+1.)*np.pi/2./NZAA)
x = np.linspace(0,30, NXAA)
y = np.linspace(0,60, NYAA)

a = np.fromfile('Volumes/linear_var'+str(var).zfill(2)+'_time'+str(time).zfill(6)+'.bin', dtype=np.float_).reshape(NXAA, NYAA, NZAA)


plt.figure()
plt.pcolormesh(x,z,a[:,0,:].T)
plt.colorbar()
plt.title('temperature fluctuations, x-z slice at first Y-gridpoint')
plt.show()

# add the background
for ix in range(NXAA):
    for iy in range(NYAA):
        a[ix,iy,:] -= z

plt.figure()
plt.pcolormesh(x,z,a[:,0,:].T)
plt.colorbar()
plt.title('full temperature, x-z slice at first Y-gridpoint')
plt.show()


