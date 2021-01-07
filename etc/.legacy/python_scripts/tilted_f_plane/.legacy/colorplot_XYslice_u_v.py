import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

   
def vertical_vorticity(u,v,Lx,Ly,NX,NY):
   kx = fft.fftfreq(NX,d=Lx/NX)
   ky = fft.fftfreq(NY,d=Ly/NY)
   KX,KY = np.meshgrid(kx,ky)
   uf = fft.fft2(u)
   vf = fft.fft2(v)
   vortf = uf*1.0j*kx + vf*1.0j*ky
   vort  = fft.ifft2(vortf).real
   return vort
   
resolution = np.fromfile('./CheckPoints/TFP_MPI_config.sav', dtype=np.int32)

path2file = './XY_Slices/'

NXAA = resolution[1]*3/2
NYAA = resolution[2]*3/2
Zstr = 'Z' + raw_input("What height do you want displayed?").zfill(4)
tstr = 't' + raw_input("What timestep do you want displayed?").zfill(8)

u = np.fromfile(path2file+'var01_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NYAA, NXAA)
v = np.fromfile(path2file+'var02_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NYAA, NXAA)


vort = vertical_vorticity(u,v,2.0,2.0,NXAA,NYAA)

plt.figure()
plt.pcolormesh(u)
plt.title('Mean velocity u (East/West)')
plt.figure()
plt.pcolormesh(v)
plt.title('Mean velocity v (South/North)')
plt.figure()
plt.pcolormesh(vort)
plt.title('Mean vertical vorticity')
plt.show()


