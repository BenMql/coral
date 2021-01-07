import numpy as np
iekman = 1.e4
NZ = T.shape[2]
grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5

GX, GY, GZ = np.meshgrid(np.arange(T.shape[0]), np.arange(T.shape[1]), grid_z)

Tfull = T/iekman**(1.0/3.0) - grid_z.astype(np.float32)
