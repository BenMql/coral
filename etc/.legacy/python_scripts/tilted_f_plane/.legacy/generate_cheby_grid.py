import numpy as np

resolutions = np.fromfile("./CheckPoints/TFP_MPI_config.sav", dtype=np.int32)

NZ = resolutions[3]*3/2

grid_z = np.arange(NZ)
grid_z = np.cos((grid_z*2+1.0)*np.pi/(2*NZ))
grid_z += 1.0
grid_z *= 0.5

NX = resolutions[1]*3/2
grid_x = np.arange(NX,dtype=np.float_)/NX

NY = resolutions[2]*3/2
grid_y = np.arange(NY,dtype=np.float_)/NY

np.savetxt('xcoords.txt', grid_x)
np.savetxt('ycoords.txt', grid_y)
np.savetxt('zcoords.txt',-grid_z)
