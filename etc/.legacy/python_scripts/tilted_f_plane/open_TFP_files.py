import numpy as np
import scipy.fftpack as fft

def compute_a_clean_mean(le_signal, le_temps):
   s_dt = le_signal[:-1] * np.diff(le_temps)
   return np.sum(s_dt)/(le_temps[-1] - le_temps[0])

class plane_layer_volume:
   def __init__(self, path_to_run, list_var_int, time_int, is_padded = False):
      if is_padded:
         padding = 2
      else:
         padding = 0
      NX = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[1]
      NY = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[2]
      NZ = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[3]
      self.NX = 3*NX/2
      self.NY = 3*NY/2
      self.NZ = 3*NZ/2
      self.dat=[]
      self.profiles=[]
      self.lut=[]
      self.lut_profiles=[]
      for ivar in list_var_int:
         self.dat.append(np.fromfile(path_to_run+'Volumes/Vol_var'+str(ivar).zfill(3)
                  +'_t'+str(time_int).zfill(8)+'_full.dat',
                  dtype=np.float_).reshape(self.NY, self.NX+padding, self.NZ)[:, :-padding, :])
         self.lut.append('var'+str(ivar).zfill(3))
      self.xgrid = np.linspace(0,1,self.NX, endpoint=False)
      self.ygrid = np.linspace(0,1,self.NY, endpoint=False)
      self.zgrid = 0.5+0.5*np.cos((2.*np.arange(self.NZ)+1.)*np.pi/(2.*self.NZ))

   def spectrum_FFC(self, pos_in_list):
      ''' this *function* returns the power spectum of self.dat[pos_in_list], i.e.
      a rank 3 array that contains the square modulus of the DFT-DFT-DChebyT of the input.
      todo: for a proper fourier-fourier-cosine transform, an interpolation on an 
      equispaced grid is needed.'''
      from scipy.fftpack import fft2, ifft2, dct, idct
      field_phys = self.dat[pos_in_list] # *points* to data
      field_aux  = dct(field_phys, axis=2, type=2)
      field_aux [:,:,0] *= 0.5 # normalize the first Chebyshev properly
      return np.abs(fft2(field_aux, axes=(0,1)))**2
     
      
      

   def xdiff_FD(self, pos_in_list):
      my_var = np.copy(self.dat[pos_in_list])
      self.dat.append  (np.diff(my_var,axis=1))
      self.lut.append('[d/dx]'+self.lut[pos_in_list])

   def xdiff(self, pos_in_list):
      field = np.copy(self.dat[pos_in_list])
      for iy in range(self.NY):
         for iz in range(self.NZ):
            field[iy,:,iz] = fft.diff(self.dat[pos_in_list][iy,:,iz], period=1.)
      self.dat.append  (field)
      self.lut.append('[d/dx]'+self.lut[pos_in_list])

   def ydiff(self, pos_in_list):
      field = np.copy(self.dat[pos_in_list])
      for ix in range(self.NX):
         for iz in range(self.NZ):
            field[:,ix,iz] = fft.diff(self.dat[pos_in_list][:,ix,iz], period=1.)
      self.dat.append  (field)
      self.lut.append('[d/dy]'+self.lut[pos_in_list])

   def ydiff_FD(self, pos_in_list):
      my_var = np.copy(self.dat[pos_in_list])
      self.dat.append  (np.diff(my_var,axis=0))
      self.lut.append('[d/dy]'+self.lut[pos_in_list])


   def zdiff_FD(self, pos_in_list):
      my_var = np.copy(self.dat[pos_in_list])
      my_df_full   = np.diff(my_var,axis=2)
      my_df_quart  = np.diff(my_var,axis=2)[:self.NY,:self.NX,:]
      for iy in range(self.NY):
         for ix in range(self.NX):
            my_df_quart[iy,ix,:]/=np.diff(self.zgrid)
      self.quart.append(my_df_quart)
      for iy in range(2*self.NY):
         for ix in range(2*self.NX):
            my_df_full [iy,ix,:]/=np.diff(self.zgrid)
      self.dat.append(my_df_full)
      self.lut.append('[d/dz]'+self.lut[pos_in_list])

   def zdiff(self, pos_in_list):
      from cheby_tools import chebyshev_elementary_integration
      from scipy.sparse.linalg import factorized
      import scipy.sparse as sp
      cheb_I = chebyshev_elementary_integration(N=self.NZ, center = 0.5, gap=1.)
      field_spec = 2*fft.dct(self.dat[pos_in_list], axis=2, type=2)
      #field_spec[:,:,0]/=2.
      deriv_spec = np.zeros(field_spec.shape, dtype=np.float_)
      solve = factorized(cheb_I.todense()[1:,:-1])
      for iy in range(self.NY):
         for ix in range(self.NX):
            deriv_spec[iy,ix,:-1] = 0.5*solve(field_spec[iy,ix,1:])
      deriv_spec[:,:,0]*=2.
      deriv_phys = fft.idct(deriv_spec, axis=2)/self.NZ
      self.dat.append(deriv_phys)
      self.lut.append('[d/dz]'+self.lut[pos_in_list])

   def zint(self, pos_in_list, gap = 1.0, center = 0.5):
      from cheby_tools import chebyshev_elementary_integration
      CEI = chebyshev_elementary_integration(gap = gap, center = center, N = self.NZ)
      field_spec = 2*fft.dct(self.profiles[pos_in_list], type=2)
      field_spec/= 0.5
      field_spec[0]/=2.
      integral_spec = np.zeros(field_spec.shape, dtype=np.float_)
      CEIdense = CEI.todense()
      CEIdense[0,:] = 0.
      integral_spec = CEIdense.dot(field_spec)[0,:]
      integral_spec*=0.5
      integral_spec[0]*=2.
      integral_phys = fft.idct(integral_spec)/self.NZ
      self.profiles.append(integral_phys)


   def compute_diUj2(self, pos_u=0, pos_v=1, pos_w=2):
      init_length = len(self.dat)
      self.xdiff(pos_in_list=pos_u)
      self.xdiff(pos_in_list=pos_v)
      self.xdiff(pos_in_list=pos_w)
      self.ydiff(pos_in_list=pos_u)
      self.ydiff(pos_in_list=pos_v)
      self.ydiff(pos_in_list=pos_w)
      self.zdiff(pos_in_list=pos_u)
      self.zdiff(pos_in_list=pos_v)
      self.zdiff(pos_in_list=pos_w)
      aux = np.zeros(self.dat[1].shape, dtype=np.float_)
      for i in range(9):
         aux += self.dat[init_length+i]**2
      self.dat.append(np.copy(aux))
      self.lut.append('Sum{i,j=1..3} [dU_i / dx_j]**2')
      self.profiles.append(np.sum(aux, axis=(0,1))/self.NY/self.NX)
      self.lut_profiles.append('Int dx dy Sum{i,j=1..3} [dU_i / dx_j]**2 = f(z)')
      prof_length = len(self.profiles)
      self.zint(prof_length-1)
      self.lut_profiles.append('Int dx dy dz (Sum{i,j=1..3} [dU_i / dx_j]**2)')

   def horizontally_average(self, pos_in_list):
      self.profiles.append(np.sum(self.dat[pos_in_list], axis=(0,1))/self.NY/self.NX)
      self.lut_profiles.append('Int dx dy '+self.lut[pos_in_list])




class cube_volume:
   def __init__(self, path_to_run, list_var_int, time_int):
      NX = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[1]
      NY = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[2]
      NZ = np.fromfile(path_to_run+'CheckPoints/TFP_MPI_config.sav', dtype=np.int32)[3]
      self.NX = 3*NX/4
      self.NY = 3*NY/4
      self.NZ = 3*NZ/2
      self.quart=[]
      self.lut=[]
      for ivar in list_var_int:
         self.quart.append(np.fromfile(path_to_run+'Volumes/Vol_var'+str(ivar).zfill(3)
                  +'_t'+str(time_int).zfill(8)+'_quart.dat',
                  dtype=np.float_).reshape(self.NY, self.NX, self.NZ))
         self.lut.append('var'+str(ivar).zfill(3))
      self.xgrid = np.linspace(0,1,self.NX, endpoint=False)
      self.ygrid = np.linspace(0,1,self.NY, endpoint=False)
      self.zgrid = 0.5+0.5*np.cos((2.*np.arange(self.NZ)+1.)*np.pi/(2.*self.NZ))



class TFP_timeseries:
   def __init__(self,list_of_repositories):
      time_list = []
      u2_list = []
      v2_list = []
      w2_list = []
      EK_list = []
      Ttop_list = []
      Tbot_list = []
      Tmid_list = []
      for irepo in range(len(list_of_repositories)):
         time_list.append(np.fromfile(list_of_repositories[irepo]+'/Timeseries/time.dat',
                      dtype=np.float_))
         count = time_list[irepo].shape[0]
         u2_list.append(np.fromfile(list_of_repositories[irepo]+'/Timeseries/u2mean_full.dat',
                      dtype=np.float_, count=count))
         v2_list.append(np.fromfile(list_of_repositories[irepo]+'/Timeseries/v2mean_full.dat',
                      dtype=np.float_, count=count))
         w2_list.append(np.fromfile(list_of_repositories[irepo]+'/Timeseries/w2mean_full.dat',
                      dtype=np.float_, count=count))
         EK_list.append(np.fromfile(list_of_repositories[irepo]+'/Timeseries/KE_total_full.dat',
                      dtype=np.float_, count=count))
         Ttop_list.append(np.fromfile(list_of_repositories[irepo]+
                     '/Timeseries/temperature_top_full.dat',  dtype=np.float_, count=count))
         Tbot_list.append(np.fromfile(list_of_repositories[irepo]+
                     '/Timeseries/temperature_bottom_full.dat',  dtype=np.float_, count=count))
         Tmid_list.append(np.fromfile(list_of_repositories[irepo]+
                     '/Timeseries/temperature_mid_depth_full.dat',  dtype=np.float_, count=count))
      self.time = np.concatenate(tuple(time_list))
      self.u2   = np.concatenate(tuple(u2_list))
      self.v2   = np.concatenate(tuple(v2_list))
      self.w2   = np.concatenate(tuple(w2_list))
      self.EK   = np.concatenate(tuple(EK_list))
      self.Ttop = np.concatenate(tuple(Ttop_list))
      self.Tbot = np.concatenate(tuple(Tbot_list))
      self.Tmid = np.concatenate(tuple(Tmid_list))


   def compute_means(self):
      nsignal = self.time.shape[0]
      nslice  = nsignal/4
      for islice in range(4):
         print('============================================')
         print('       Quarter '+str(islice)+'/4')
         print('--------------------------------------------')
         print(' Mean TOTAL Kinetic Energy: '
               +str(compute_a_clean_mean(self.EK  [islice*nslice: (islice+1)*nslice],
                                         self.time[islice*nslice: (islice+1)*nslice])))
         print(' Mean HORIZONTAL Kinetic Energy: '
               +str(compute_a_clean_mean(self.u2  [islice*nslice: (islice+1)*nslice]+
                                         self.v2  [islice*nslice: (islice+1)*nslice],
                                         self.time[islice*nslice: (islice+1)*nslice])))
         print(' Mean VERTICAL Kinetic Energy: '
               +str(compute_a_clean_mean(self.w2  [islice*nslice: (islice+1)*nslice],
                                         self.time[islice*nslice: (islice+1)*nslice])))
         print(' Average Reynolds number:'
               +str(np.sqrt(compute_a_clean_mean(self.EK  [islice*nslice: (islice+1)*nslice],
                                                 self.time[islice*nslice: (islice+1)*nslice]))))
         print(' < T_bottom - T_mid > = '
               +str(compute_a_clean_mean(self.Tbot[islice*nslice: (islice+1)*nslice] -
                                         self.Tmid[islice*nslice: (islice+1)*nslice],
                                         self.time[islice*nslice: (islice+1)*nslice])))
         print(' < T_top    - T_mid > = '
               +str(compute_a_clean_mean(self.Ttop[islice*nslice: (islice+1)*nslice] -
                                         self.Tmid[islice*nslice: (islice+1)*nslice],
                                         self.time[islice*nslice: (islice+1)*nslice])))

   def plot(self):
      import matplotlib.pyplot as plt
      plt.figure()
      plt.plot(self.time, self.Ttop-self.Tmid, label=r'T_{top} - T_{mid}')
      plt.plot(self.time, self.Tbot-self.Tmid, label=r'T_{bot} - T_{mid}')
      plt.xlabel('Time')
      plt.ylabel('Temperature')
      plt.legend()
      plt.figure()
      plt.plot(self.time, np.sqrt(self.EK),           label=r'KE, total')
      plt.plot(self.time, np.sqrt(self.u2 + self.v2), label=r'KE, horizontal')
      plt.plot(self.time, np.sqrt(self.w2),           label=r'KE, vertical')
      plt.xlabel('Time')
      plt.ylabel('Reynolds')
      plt.legend()
      plt.show()

