import numpy as np
import matplotlib.pyplot as plt

def compute_a_clean_mean(le_signal, le_temps):
   s_dt = le_signal[:-1] * np.diff(le_temps)
   return np.sum(s_dt)/(le_temps[-1] - le_temps[0])
prandtl  = 1.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_2_files = './'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p2f= path_2_files + 'Timeseries/'

time = np.fromfile(p2f+'time.dat', dtype = np.float_)[:-1]
nelems= time.shape[0]
ista = nelems/4


full_KE = np.fromfile(p2f+'KE_total_full.dat', dtype = np.float_, count=nelems) / prandtl **2

T_top = np.fromfile(p2f+'temperature_top_full.dat',        dtype = np.float_, count = nelems)
T_bot = np.fromfile(p2f+'temperature_bottom_full.dat',     dtype = np.float_, count = nelems)
T_mid = np.fromfile(p2f+'temperature_mid_depth_full.dat',  dtype = np.float_, count = nelems)
T_avg = np.fromfile(p2f+'temperature_average_full.dat',    dtype = np.float_, count = nelems)
T_msq = np.fromfile(p2f+'temperature_mean_square_full.dat',dtype = np.float_, count = nelems)

u2avg = np.fromfile(p2f+'u2mean_full.dat', dtype=np.float_, count=nelems) / prandtl **2
v2avg = np.fromfile(p2f+'v2mean_full.dat', dtype=np.float_, count=nelems) / prandtl **2
w2avg = np.fromfile(p2f+'w2mean_full.dat', dtype=np.float_, count=nelems) / prandtl **2
u2top = np.fromfile(p2f+'u2top_full.dat',  dtype=np.float_, count=nelems) / prandtl **2
v2top = np.fromfile(p2f+'v2top_full.dat',  dtype=np.float_, count=nelems) / prandtl **2

#plt.figure()
#plt.semilogy(np.diff(time), 'o')
#plt.title('Time-step size')

plt.figure()
plt.semilogy(time, full_KE, label='Total KE')
plt.semilogy(time, u2avg, label='<u2> vol')
plt.semilogy(time, v2avg, label='<v2> vol')
plt.semilogy(time, w2avg, label='<w2> vol')
plt.semilogy(time, u2top, label='<u2> top')
plt.semilogy(time, v2top, label='<v2> top')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('KE')
plt.legend()
plt.title('Kinetic energy densities')
plt.grid(True)
plt.grid(True, which='minor', linestyle=':')


plt.figure()
plt.plot(time, T_top, label='Temperature top')
plt.plot(time, T_mid, label='Temperature mid')
plt.plot(time, T_bot, label='Temperature bot')
plt.plot(time, T_avg, label='Temperature avg')
plt.plot(time, np.sqrt(T_msq), label='Temperature rms')
plt.xlabel('(Rescaled) Viscous Time')
plt.ylabel('Temperature')
plt.title('Averaged temperatures')
plt.legend()
plt.grid(True)
#plt.legend(loc=2)

print(' ===================================================' )
print('              *** FIRST QUARTER ***')
print(' ===================================================' )
print('tmid - tbot:   '+ str(compute_a_clean_mean(T_mid[:ista]-T_bot[:ista], time[:ista])))
print('u2 + v2 top:   '+ str(compute_a_clean_mean(u2top[:ista]+v2top[:ista], time[:ista])))
print('bulk KE:       '+ str(compute_a_clean_mean(full_KE[:ista], time[:ista])))
print('Reynolds top:  '+ str(np.sqrt(compute_a_clean_mean(u2top[:ista]+v2top[:ista], time[:ista]))))
print('Reynolds vol:  '+ str(np.sqrt(compute_a_clean_mean(full_KE[:ista], time[:ista]))))
print(' ===================================================' )
print('              *** SECOND QUARTER ***')
print(' ===================================================' )
print('tmid - tbot:   '+ str(compute_a_clean_mean(T_mid[ista:2*ista]-T_bot[ista:2*ista], time[ista:2*ista])))
print('u2 + v2 top:   '+ str(compute_a_clean_mean(u2top[ista:2*ista]+v2top[ista:2*ista], time[ista:2*ista])))
print('bulk KE:       '+ str(compute_a_clean_mean(full_KE[ista:2*ista], time[ista:2*ista])))
print('Reynolds top:  '+ str(np.sqrt(compute_a_clean_mean(u2top[ista:2*ista]+v2top[ista:2*ista], time[ista:2*ista]))))
print('Reynolds vol:  '+ str(np.sqrt(compute_a_clean_mean(full_KE[ista:2*ista], time[ista:2*ista]))))



print(' ===================================================' )
print('              *** THIRD QUARTER ***')
print(' ===================================================' )
print('tmid - tbot:   '+ str(compute_a_clean_mean(T_mid[2*ista:3*ista]-T_bot[2*ista:3*ista], time[2*ista:3*ista])))
print('u2 + v2 top:   '+ str(compute_a_clean_mean(u2top[2*ista:3*ista]+v2top[2*ista:3*ista], time[2*ista:3*ista])))
print('bulk KE:       '+ str(compute_a_clean_mean(full_KE[2*ista:3*ista], time[2*ista:3*ista])))
print('Reynolds top:  '+ str(np.sqrt(compute_a_clean_mean(u2top[2*ista:3*ista]+v2top[2*ista:3*ista], time[2*ista:3*ista]))))
print('Reynolds vol:  '+ str(np.sqrt(compute_a_clean_mean(full_KE[2*ista:3*ista], time[2*ista:3*ista]))))
print(' ===================================================' )
print('              *** FOURTH QUARTER ***')
print(' ===================================================' )
print('tmid - tbot:   '+ str(compute_a_clean_mean(T_mid[3*ista:4*ista]-T_bot[3*ista:4*ista], time[3*ista:4*ista])))
print('u2 + v2 top:   '+ str(compute_a_clean_mean(u2top[3*ista:4*ista]+v2top[3*ista:4*ista], time[3*ista:4*ista])))
print('bulk KE:       '+ str(compute_a_clean_mean(full_KE[3*ista:4*ista], time[3*ista:4*ista])))
print('Reynolds top:  '+ str(np.sqrt(compute_a_clean_mean(u2top[3*ista:4*ista]+v2top[3*ista:4*ista], time[3*ista:4*ista]))))
print('Reynolds vol:  '+ str(np.sqrt(compute_a_clean_mean(full_KE[3*ista:4*ista], time[3*ista:4*ista]))))

plt.show()


