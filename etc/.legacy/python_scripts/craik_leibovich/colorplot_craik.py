import numpy as np
import matplotlib.pyplot as plt
import read_CLR_input_file as rif

def vertical_vorticity(u,v,Lx,Ly,NX,NY):
   kx = fft.fftfreq(NX,d=Lx/NX)
   ky = fft.fftfreq(NY,d=Ly/NY)
   KX,KY = np.meshgrid(kx,ky)
   uf = fft.fft2(u)
   vf = fft.fft2(v)
   vortf = uf*1.0j*kx + vf*1.0j*ky
   vort  = fft.ifft2(vortf).real
   return vort

   

NX =rif.read_NX()/2*3
NZ =rif.read_NZ()/2*3
depth = rif.read_depth()
lbox_x= rif.read_Lbox_X()

t_str = raw_input('Display what time? >').zfill(8)
v_str = raw_input('Display what var.? >').zfill(3)
#t_str = str(3000).zfill(8)

T = np.fromfile('./XZ_Slices/var'+v_str+'_t'+t_str+'_full.dat', dtype=np.float_).reshape(NZ, NX)

grid_z = np.arange(NZ, dtype=np.float_)*2+1
grid_z = np.cos(grid_z*np.pi/2./NZ)
grid_z*= 0.5
grid_z-= 0.5
grid_z*= depth

grid_x = np.linspace(0.,lbox_x, num=NX)

plt.pcolormesh(grid_x, grid_z, T, cmap='RdYlBu_r')
plt.colorbar()
plt.show()




