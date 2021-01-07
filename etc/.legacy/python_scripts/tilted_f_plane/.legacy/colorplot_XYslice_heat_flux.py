import numpy as np
import read_TFP_input_file as rif
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
Zstr = 'Z' + raw_input("What height do you want displayed?   > ").zfill(4)
tstr = 't' + raw_input("What timestep do you want displayed? > ").zfill(8)
ekman = 1./rif.read_iekman()

nu = np.fromfile(path2file+'var61_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NYAA, NXAA)
nu-=1
nu*=ekman**(1.0/3.0)
nu+=1


print ("Nu max is "+str(nu.max()))
print ("Nu min is "+str(nu.min()))
print ("Nu average is "+str(np.mean(nu)))

plt.figure()
plt.pcolormesh(nu, cmap = 'inferno')
plt.title('local Nusselt')
plt.colorbar()
plt.figure()
plt.plot(nu[0,:])
plt.plot(nu[0,NXAA/2:],'--')
plt.show()


