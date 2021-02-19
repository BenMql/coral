import numpy as np
import matplotlib.pyplot as plt

# =========================================================================
# Read timeseries
#..........................................................................
time = np.fromfile('Timeseries/time.dat',  dtype=np.float_)
Ttop = np.fromfile('Timeseries/tFull_XYavg_z00001.dat', dtype=np.float_)[:-2]
Tbot = np.fromfile('Timeseries/tFull_XYavg_z00072.dat', dtype=np.float_)[:-2]

quarterLength = time.shape[0]//4

# =========================================================================
# compute dz
#..........................................................................
dz = 1.-0.5 - 0.5*np.cos(np.pi/2/72.)

nusseltTop_1stQuarter = np.sum(np.diff(time)[0*quarterLength:1*quarterLength]*
                                        Ttop[0*quarterLength:1*quarterLength]
                              ) / np.sum(np.diff(time)[0*quarterLength:1*quarterLength])/dz
nusseltBot_1stQuarter =-np.sum(np.diff(time)[0*quarterLength:1*quarterLength]*
                                        Tbot[0*quarterLength:1*quarterLength]
                              ) / np.sum(np.diff(time)[0*quarterLength:1*quarterLength])/dz
nusseltTop_2ndQuarter = np.sum(np.diff(time)[1*quarterLength:2*quarterLength]*
                                        Ttop[1*quarterLength:2*quarterLength]
                              ) / np.sum(np.diff(time)[1*quarterLength:2*quarterLength])/dz
nusseltBot_2ndQuarter =-np.sum(np.diff(time)[1*quarterLength:2*quarterLength]*
                                        Tbot[1*quarterLength:2*quarterLength]
                              ) / np.sum(np.diff(time)[1*quarterLength:2*quarterLength])/dz
nusseltTop_3rdQuarter = np.sum(np.diff(time)[2*quarterLength:3*quarterLength]*
                                        Ttop[2*quarterLength:3*quarterLength]
                              ) / np.sum(np.diff(time)[2*quarterLength:3*quarterLength])/dz
nusseltBot_3rdQuarter =-np.sum(np.diff(time)[2*quarterLength:3*quarterLength]*
                                        Tbot[2*quarterLength:3*quarterLength]
                              ) / np.sum(np.diff(time)[2*quarterLength:3*quarterLength])/dz
nusseltTop_4thQuarter = np.sum(np.diff(time)[3*quarterLength:4*quarterLength]*
                                        Ttop[3*quarterLength:4*quarterLength]
                              ) / np.sum(np.diff(time)[3*quarterLength:4*quarterLength])/dz
nusseltBot_4thQuarter =-np.sum(np.diff(time)[3*quarterLength:4*quarterLength]*
                                        Tbot[3*quarterLength:4*quarterLength]
                              ) / np.sum(np.diff(time)[3*quarterLength:4*quarterLength])/dz


mean_nusselt = (nusseltTop_2ndQuarter+
                nusseltBot_2ndQuarter+
                nusseltTop_3rdQuarter+
                nusseltBot_3rdQuarter+
                nusseltTop_4thQuarter+
                nusseltBot_4thQuarter)
mean_nusselt /= 6.
std_nusselt = np.sqrt(1./6.*(
                (nusseltTop_2ndQuarter-mean_nusselt)**2+
                (nusseltBot_2ndQuarter-mean_nusselt)**2+
                (nusseltTop_3rdQuarter-mean_nusselt)**2+
                (nusseltBot_3rdQuarter-mean_nusselt)**2+
                (nusseltTop_4thQuarter-mean_nusselt)**2+
                (nusseltBot_4thQuarter-mean_nusselt)**2))

print( '===================================================')
print( '===================================================')
print( 'nusselt top (1/4):    '+str(nusseltTop_1stQuarter))
print( 'nusselt bottom (1/4): '+str(nusseltBot_1stQuarter))
print( 'nusselt top (2/4):    '+str(nusseltTop_2ndQuarter))
print( 'nusselt bottom (2/4): '+str(nusseltBot_2ndQuarter))
print( 'nusselt top (3/4):    '+str(nusseltTop_3rdQuarter))
print( 'nusselt bottom (3/4): '+str(nusseltBot_3rdQuarter))
print( 'nusselt top (4/4):    '+str(nusseltTop_4thQuarter))
print( 'nusselt bottom (4/4): '+str(nusseltBot_4thQuarter))
print( '===================================================')
print( '===================================================')
print( 'Neglecting the first quarter of the timeseries,')
print( 'the mean Nusselt number is: '+str(mean_nusselt) )
print( 'and the standard deviation: '+str(std_nusselt)  )

plt.plot(time[:-1], Ttop/dz, label=r'$\left.dT/dz\right\|_{z=1}$')
plt.plot(time[:-1],-Tbot/dz, label=r'$\left.dT/dz\right\|_{z=0}$')
plt.xlabel('diffusive time')
plt.title ('Instantaneous heat flux')
plt.legend()
plt.grid(True)
plt.ylim((0,8.))
plt.show()
