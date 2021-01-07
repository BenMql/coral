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
w = np.fromfile(path2file+'var03_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NYAA, NXAA)
T = np.fromfile(path2file+'var50_'+Zstr+'_'+tstr+'_full.dat', dtype=np.float_).reshape(NYAA, NXAA)

print('u mean is: '+str(np.mean(u)))
print('v mean is: '+str(np.mean(v)))
print('w mean is: '+str(np.mean(w)))
print('u std  is: '+str(np.std(u)))
print('v std  is: '+str(np.std(v)))
print('w std  is: '+str(np.std(w)))

vort = vertical_vorticity(u,v,2.0,2.0,NXAA,NYAA)

# ~~~~ DISPLAY VORTICITY ~~~~~~
plt.figure()
plt.pcolormesh(vort)
plt.title('Vertical vorticity')
# ~~~~ DISPLAY VELOCITY W ~~~~~
plt.figure()
plt.pcolormesh(w)
plt.title('Vertical velocity')
# ~~~~ DISPLAY TEMPERATURE ~~~~
plt.figure()
T-=np.mean(T)
tmax = np.abs(T).max()
#plt.pcolormesh(T, vmin = -tmax, vmax = tmax, cmap='PuOr')
plt.pcolormesh(T, cmap='RdYlGn')
plt.title('Temperature Fluctuation')
plt.show()




