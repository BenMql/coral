import numpy as np
import matplotlib.pyplot as plt
import read_TFP_input_file as rif

resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)
eta = 10.

NY = resolution[2]*3/2
NZ = resolution[3]*3/2

latitude = rif.read_latitude()
Ly_box = rif.read_Lbox_Y()
iekman = rif.read_iekman()                               
if (iekman>0.):
   ekman = 1./rif.read_iekman()                             
else:
    ekman = 1.
ymax = Ly_box*ekman**(1./3.)
eta = rif.read_H_over_ell()
radHeat_bool = rif.read_radiativeHeat_bool()

if (eta<1.e-8):
    H_over_ell_is_0 = True
else:
    H_over_ell_is_0 = False

tstr = 't'+raw_input("Display what time? > ").zfill(8)
vstr = 'var'+raw_input('Display what var?  > ').zfill(3)
xstr = 'X'+raw_input('Display what X?    > ').zfill(4)

grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5
grid_y = np.arange(NY)* ymax /NY

GY, GZ = np.meshgrid(grid_y, grid_z)

T = np.fromfile('./YZ_Slices/'+vstr+'_'+xstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NY,NZ).T
if (radHeat_bool):
    if (H_over_ell_is_0):
        T_back = GZ**2/2.0 - GZ
    else:
        T_back = GZ**2/2.0 + (GZ + 1./eta*np.exp(-eta*GZ))/(np.exp(-eta)-1.)
else:
    T_back = -GZ
#T += T_back

filter_mean = raw_input("Should I filter the mean? > ")
if (filter_mean=='yes'):
   Tmean1 = np.mean(T,axis = 1)
   todel, Tmean2 = np.meshgrid(grid_y, Tmean1)
   T-=Tmean2

f1 = plt.figure(figsize=(6,6))
if (ymax<1.0):
  # ax_width = 0.8*ymax
  # ax_shift = 0.5-ax_width/2.0
   ax_width = 0.5 
   ax_shift = 0.5-ax_width/2.0
   ax = f1.add_axes([ax_shift, 0.1, ax_width, 0.8])
else:            
   ax_width = 0.8/ymax
   ax_shift = 0.5-ax_width/2.0
   ax = f1.add_axes([0.1, ax_shift, 0.8, ax_width])
my_cmap = 'Paired'
my_cmap = 'RdYlBu_r'
ax.pcolormesh(GY,GZ,T, vmin = T.min(), vmax =T.max(), cmap = my_cmap)
print('T rms is: '+str(np.std(T)))
#~~ plot rotation axis ~~
if not (rif.read_upright_bool()):
   zvec = np.linspace(0,1)
   plt.plot(grid_y, grid_y*np.tan(latitude/180.*np.pi), '--', color='k')
ax.set_xlim([0, ymax])
ax.set_ylim([0, 1   ])
ax.set_xlabel(r'latitude $y$')
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]) 
ax.set_xticks([x for x in [0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] if x<ymax])
ax.set_ylabel(r'profondeur $z$', fontsize = 16)
ax.grid(True)
plt.show()
