import numpy as np
import read_matrices as rm
import read_TFP_input_file as rif
import matplotlib.pyplot as plt

prandtl = rif.read_prandtl()
iekman  = rif.read_iekman()
rayleigh = rif.read_rayleigh()*iekman**(4./3.)

path2mat = '/data/plane_layer/E5_R40_P1_A02/Matrices'
i_core =  1
i_kx   =  4 
i_ky   =  1

L = rm.coral_matrix(path2mat+'/stif_matrix_L_mx'
                    + str(i_kx).zfill(4) +
              '_ky' + str(i_ky).zfill(4) + 
           '.core' + str(i_core).zfill(4) , 'z')

M = rm.coral_matrix(path2mat+'/mass_matrix_M_mx'
                    + str(i_kx).zfill(4) +
              '_ky' + str(i_ky).zfill(4) + 
           '.core' + str(i_core).zfill(4) , 'z')

L.compute_full_spectrum(other=M)
L.plot_full_spectrum()

# Now compute the analytic spectrum

table_of_kx = np.fromfile(path2mat+'/dx_arr_core'+str(i_core).zfill(4)+'.dat', dtype=np.float_)
table_of_ky = np.fromfile(path2mat+'/dy_arr_core'+str(i_core).zfill(4)+'.dat', dtype=np.float_)

kx = table_of_kx[i_kx-1]*iekman**(1./3.)
ky = table_of_ky[i_ky-1]*iekman**(1./3.)

m2 = kx**2 + ky**2
kz =np.arange(20)+1
D2s = m2 + np.pi**2*kz**2

s = []
for i in range(D2s.shape[0]):
   D2 = D2s[i]
   poly = np.empty(4)
   poly[0] = prandtl**(-2.)*D2
   poly[1] = D2**2/prandtl**2+ 2.* D2**2/prandtl
   poly[2] =(2. /  prandtl +1.)* D2**3 - rayleigh * m2 / prandtl + iekman**2*np.pi**2*kz[i]**2
   poly[3] = D2**4 - D2 * rayleigh * m2 + iekman**2*np.pi**2*D2*kz[i]**2

   s.append(np.roots(poly))
   plt.plot(s[-1].real/iekman**(2./3.), s[-1].imag/iekman**(2./3.),'+', color='k')

plt.show()



