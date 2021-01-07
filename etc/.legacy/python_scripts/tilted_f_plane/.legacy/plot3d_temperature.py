import numpy as np
from mayavi import mlab

resolution = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype = np.int32)

NX = resolution[1]*3/2
NY = resolution[2]*3/2
NZ = resolution[3]*3/2
Lx_box = 2.0
Ly_box = 2.0

ekman = input       ('Ekman number?      > ')              
tstr = 't'+raw_input("Display what time? > ").zfill(8)

grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5
grid_x = np.arange(NX)*Lx_box/NX
grid_y = np.arange(NY)*Lx_box/NY

GX, GY, GZ = np.meshgrid(grid_x, grid_y, grid_z)

T = np.fromfile('./Volumes/Vol_var50_'+tstr+'_full.dat', dtype=np.float_).reshape(NY,NX,NZ)*ekman**(1.0/3.0)

total_temperature_bool = raw_input("Should I display the full field (yes) or only the fluctuation (no)? > ")
if (total_temperature_bool=='yes'):
   T-=GZ

print (T.max())
print (T.min())

mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()
#source = mlab.pipeline.scalar_field(GX, GY, GZ, T)
source = mlab.pipeline.scalar_field(T)
min = T.min()
max = T.max()
#vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
                                   vmax=min + 0.95* (max - min))
#iso = mlab.pipeline.iso_surface(source, vmin=min,vmax=max, colormap = 'RdYlBu')
#mlab.contour3d(GX, GY, GZ, T)
#mlab.view(132, 54, 45, [21, 20, 21.5])
mlab.show()


