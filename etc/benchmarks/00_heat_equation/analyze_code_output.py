import numpy as np
import matplotlib.pyplot as plt

dt=1.e-3
Prandtl=3.3e-1


center = 0.5
gap = 1.
NZ=192

grid_z = center + gap/2.*np.cos( (2.*np.arange(NZ) + 1.)  *np.pi / 2./NZ) #this is Gauss-Cheby

nu = 4
nv = 3
nw = 5
u0 = 1.34e-4
v0 = -2.3e-2
w0 = np.pi*1.e-3
decay_rate = {'u': Prandtl*((2*nu + 1)*np.pi/2.)**2,
              'v': Prandtl*((2*nv + 1)*np.pi/2.)**2,
              'w': Prandtl*(nw*np.pi)**2}



u_mean = u0*np.sin( (2*nu + 1) * (grid_z -1)*np.pi/2.)
v_mean = v0*np.sin( (2*nv + 1) * grid_z *np.pi/2.)
w_mean = w0*np.sin( nw * grid_z *np.pi)

it = 10

# ... 
u1 = np.fromfile('./Profiles/linear_var01_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
v1 = np.fromfile('./Profiles/linear_var02_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
w1 = np.fromfile('./Profiles/linear_var03_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
# ... 
it= 20
time2 = (it-1)*dt + 2.e-5
u2 = np.fromfile('./Profiles/linear_var01_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
v2 = np.fromfile('./Profiles/linear_var02_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
w2 = np.fromfile('./Profiles/linear_var03_time'+str(it).zfill(8)+'_full.dat', dtype=np.float_)
# ... 
plt.figure()
plt.plot(grid_z, u1, 'o-', color='C0')
plt.plot(grid_z, u2, 'v-', color='C1')
plt.plot(grid_z, u1*np.exp(-decay_rate['u']*10*dt), 'k')
plt.xlabel(r'z', fontsize=16)
plt.title(r'u(z)', fontsize=18)
plt.grid(True)
plt.figure()
plt.plot(grid_z, v1, 'o-', color='C0')
plt.plot(grid_z, v2, 'v-', color='C1')
plt.plot(grid_z, v1*np.exp(-decay_rate['v']*10*dt), 'k')
plt.xlabel(r'z', fontsize=16)
plt.title(r'v(z)', fontsize=18)
plt.grid(True)
plt.figure()
plt.plot(grid_z, w1, 'o-', color='C0')
plt.plot(grid_z, w2, 'v-', color='C1')
plt.plot(grid_z, w1*np.exp(-decay_rate['w']*10*dt), 'k')
plt.xlabel(r'z', fontsize=16)
plt.title(r'w(z)', fontsize=18)
plt.grid(True)
plt.show()

print ('relative error over 10 steps: '
       +str( np.sqrt( np.sum((u2 -    
                 u1*np.exp(-decay_rate['u']*10*dt))**2) / 
                      np.sum(u1**2)
                     )))

