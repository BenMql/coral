import numpy as np
import matplotlib.pyplot as plt

'''
We open the timeseries produced by running the stress-free
Rayleigh-Benard benchmark run (benchmark_01). We compute the 
growthrate obtained numerically during the linear growth 
stage to the analytical prediction
'''

# Open the files, construct the kinetic energy, and compute
# its time-derivative.

argv=['uu','vv','ww']

plt.figure()
time = np.fromfile('./Timeseries/time.dat', dtype=np.float_)
uu = np.fromfile('./Timeseries/uu_volAvg.dat', dtype=np.float_)
vv = np.fromfile('./Timeseries/vv_volAvg.dat', dtype=np.float_)
ww = np.fromfile('./Timeseries/ww_volAvg.dat', dtype=np.float_)
KE = uu + vv + ww

def centered_time_derivative(signal, t):
   left_Derivative = (signal[1:-1] - signal[0:-2])/(t[1:-1] - t[0:-2])
   rightDerivative = (signal[2:]   - signal[1:-1])/(t[2:]   - t[1:-1])
   dfdt = left_Derivative + 0.5*(t[1:-1] - t[0:-2]) * (
                   (rightDerivative-left_Derivative) / ((t[2:]- t[0:-2])*0.5))
   return dfdt
                                                                

dEdt = centered_time_derivative(KE, time)

# plot the (numerical) growthrate
plt.semilogy(time[1:-1], dEdt/KE[1:-1], 'C0', linewidth=1, label=r'$\frac{1}{E}\frac{dE}{dt}$')
plt.semilogy(time[1:-1],-dEdt/KE[1:-1], 'C0', linewidth=1,linestyle='dashed',label=r'$-\frac{1}{E}\frac{dE}{dt}$')
plt.title('growthrate')


# Now evaluate (analytically) the maximum growthrate, for:
rayleigh = 2000.
prandtl = 7.
Lx = 10.           
Ly = 10.           

NX = 96*3//2
NY = 96*3//2

kx, ky = np.meshgrid(np.linspace(0,2*np.pi/Lx*NX,num=NX, endpoint=False),
                     np.linspace(0,2*np.pi/Lx*NY,num=NY, endpoint=False))


D2 = kx**2 + ky**2 + np.pi**2
k  = np.sqrt(kx**2 + ky**2)
s = -D2*(prandtl+1.) + np.sqrt((prandtl+1)**2*D2**2 - 4.*(prandtl*D2**2 - rayleigh*prandtl* k**2 /D2))
s*= 0.5

print('maximum growthrate s = '+ str(s.max()))
print('On this grid: (exact values may differ on different/finer grids)')
print('optimal wavenumber k = '+str(np.sqrt(D2[s==s.max()][0]-np.pi**2)))
print('                   kx= '+str(        kx[s==s.max()][0]          ))
print('                   ky= '+str(        ky[s==s.max()][0]          ))
print('optimal wavenumber l = '+str(2*np.pi/np.sqrt(D2[s==s.max()][0]-np.pi**2)))

# add this information on the energy growthrate plot
# (bearing in mind that energy grows twice as fast as amplitude)

plt.semilogy([0.,1.],[s.max()*2., s.max()*2], 'k', linewidth=1,linestyle='-.', label=r'Analytic $\gamma$')
plt.xlim([0.,3.])
plt.ylim([1.e-3,1.e2])
plt.legend()
plt.grid(True)
plt.grid(True, which='minor', linestyle='dotted')

plt.figure()
plt.plot(np.sqrt(kx**2 + ky**2), s, '+', color='C0')
# mark the optimum
plt.plot(np.sqrt(D2[s==s.max()][0]-np.pi**2), s.max(), 'o', color='C1')
plt.xlabel(r'horizontal wavenumber $k=\sqrt{k_x+k_y}$')
plt.ylabel('growthrate s(k)')
plt.xlim([0, 6.5])
plt.ylim([-20, 30])
plt.grid(True)
plt.show()

plt.show()

