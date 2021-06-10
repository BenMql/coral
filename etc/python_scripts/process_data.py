import numpy as np
import scipy.fftpack as fft
import get_coral_params as gcp
import matplotlib.pyplot as plt

def compute_a_clean_mean(le_signal, le_temps):
   s_dt = le_signal[:-1] * np.diff(le_temps)
   return np.sum(s_dt)/(le_temps[-1] - le_temps[0])

class plane_layer_volume:
   def __init__(self, list_var_str, time_int,path_to_run = './') :
      ddecomp = gcp.domDecomp(0) 
      NX = ddecomp.NX                                                                     
      NY = ddecomp.NY                                                                     
      NZ = ddecomp.NZ                                                                     
      self.dat=[]
      self.profiles=[]
      self.verticalSums=[]
      self.lut=[]
      self.lut_profiles=[]
      self.lut_verticalSums=[]
      self.lut=[]
      svar = list_var_str[0]
      geom = gcp.coralGeom(0) 
      self.lx = geom.Lx                           
      self.ly = geom.Ly                             
      self.lz = geom.gap                           
      try:
         self.NX = NX
         self.NY = NY
         self.NZ = NZ
         padding = 0
         suffix = '_twoThirds.dat'
         myArr = np.fromfile(path_to_run+'Volumes/'+list_var_str[0]+'_time'
                             +str(time_int).zfill(8)+suffix,
                  dtype=np.float_).reshape(
                        self.NY, self.NX, self.NZ)
      except IOError:
         self.NX = 3*NX//2
         self.NY = 3*NY//2
         self.NZ = 3*NZ//2
         suffix = '_full.dat'
         myArr = np.fromfile(path_to_run+'Volumes/'+list_var_str[0]+'_time'
                                +str(time_int).zfill(8)+suffix,
                     dtype=np.float_).reshape(
                           self.NX, self.NY, self.NZ)
      del myArr
      for ivar in range(len(list_var_str)):
         self.dat.append(
              np.fromfile(path_to_run+'Volumes/'+list_var_str[ivar]+'_time'
                             +str(time_int).zfill(8)+suffix,
                  dtype=np.float_).reshape(
                        self.NX, self.NY, self.NZ))
         self.lut.append(list_var_str[ivar])
      self.xgrid = np.linspace(0,self.lx,self.NX, endpoint=False)
      self.ygrid = np.linspace(0,self.ly,self.NY, endpoint=False)
      self.zgrid = 0.5+0.5*np.cos((2.*np.arange(self.NZ)+1.)*np.pi/(2.*self.NZ))*self.lz

   def plot_slice(self, varInt= 0, x=-100000, y=-100000, z=-100000):
       '''
       Display a slice of self.dat[varInt] taken at a given x, or y, or z.
       Call to this routine must specify only one x, y, or z location (as suitable for a slice).
       '''
       # check that the call has the proper form:
       if ((x==-100000) and (y==-100000)):
           if not(-self.NZ < z < self.NZ):
              print ('incorrect input')
              return 0
           else: 
              myFig=plt.figure()
              plt.pcolormesh(self.xgrid, self.ygrid, self.dat[varInt][:,:,z].T, shading='auto')
              plt.title(self.lut[varInt]+', index z='+str(z))
              return myFig
       # check that the call has the proper form:
       elif ((x==-100000) and (z==-100000)):
           if not(-self.NY < y < self.NY):
              print ('incorrect input')
              return 0
           else: 
              myFig=plt.figure()
              plt.pcolormesh(self.xgrid, self.zgrid, self.dat[varInt][:,y,:].T)
              plt.title(self.lut[varInt]+', index y='+str(y))
              return myFig
       elif ((y==-100000) and (z==-100000)):
           if not(-self.NX < x < self.NX):
              print ('incorrect input')
              return 0
           else: 
              myFig=plt.figure()
              plt.pcolormesh(self.ygrid, self.zgrid, self.dat[varInt][x,:,:].T)
              plt.title(self.lut[varInt]+', index x='+str(x))
              return myFig
       else:
              print ('incorrect input')
              return 0


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
      self.dat.append  (np.diff(my_var,axis=1)/self.lx*self.NX)
      self.lut.append('[d/dx]'+self.lut[pos_in_list])

   def xdiff(self, pos_in_list):
      field = np.copy(self.dat[pos_in_list])
      for iy in range(self.NY):
         for iz in range(self.NZ):
            field[iy,:,iz] = fft.diff(self.dat[pos_in_list][iy,:,iz], period=self.lx)
      self.dat.append  (field)
      self.lut.append('[d/dx]'+self.lut[pos_in_list])

   def ydiff(self, pos_in_list):
      field = np.copy(self.dat[pos_in_list])
      for ix in range(self.NX):
         for iz in range(self.NZ):
            field[:,ix,iz] = fft.diff(self.dat[pos_in_list][:,ix,iz], period=self.ly)
      self.dat.append  (field)
      self.lut.append('[d/dy]'+self.lut[pos_in_list])

   def ydiff_FD(self, pos_in_list):
      my_var = np.copy(self.dat[pos_in_list])
      self.dat.append  (np.diff(my_var,axis=0)/self.ly*self.NY)
      self.lut.append('[d/dy]'+self.lut[pos_in_list])


   def zdiff(self, pos_in_list):
      from cheby_tools import chebyshev_elementary_integration
      from scipy.sparse.linalg import factorized
      import scipy.sparse as sp
      cheb_I = chebyshev_elementary_integration(N=self.NZ, center = 0.5, gap=self.lz)
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

   def zint_Vol(self, pos_in_list):
      from cheby_tools import chebyshev_elementary_integration
      from scipy.sparse.linalg import factorized
      import scipy.sparse as sp
      cheb_I = chebyshev_elementary_integration(N=self.NZ, center = 0.5, gap=self.lz)
      field_spec = 2*fft.dct(self.dat[pos_in_list], axis=2, type=2)
      field_spec[:,:,0]/=2.
      integral_spec = np.zeros(field_spec.shape, dtype=np.float_)
      CEIdense = cheb_I.todense()
      CEIdense[0,:] = 0.
      for iy in range(self.NY):
         for ix in range(self.NX):
            integral_spec[iy,ix,:] = CEIdense.dot(field_spec[iy,ix,:])
      # we compute a definite integral: g(z) = int_{z'=z_bottom}^{z'= z} dz' f(z')
      # hence g(z_bottom=0)
      integral_spec[:,:,0] = 0.
      for i in range(integral_spec.shape[2]-1):
         integral_spec[:,:,0] += integral_spec[:,:,i+1] * (-1.)**(i)
      full_int = np.sum(integral_spec, axis=2)
      #normalization
      integral_spec[:,:,0]*=2.
      #integral_spec[:,:,0]/=np.sqrt(self.NZ)
      integral_phys = fft.idct(integral_spec, axis=2)/self.NZ
      self.dat.append(integral_phys)
      self.lut.append('[Int dz]'+self.lut[pos_in_list])
      return full_int/self.NZ/4.

   def zint(self, pos_in_list):
      gap = self.lz
      center = self.lz/2.
      from cheby_tools import chebyshev_elementary_integration
      CEI = chebyshev_elementary_integration(gap = gap, center = center, N = self.NZ)
      field_spec = 2*fft.dct(self.profiles[pos_in_list], type=2)
      field_spec[0]/=2.
      integral_spec = np.zeros(field_spec.shape, dtype=np.float_)
      CEIdense = CEI.todense()
      CEIdense[0,:] = 0.
      integral_spec = CEIdense.dot(field_spec)[0,:]
      copy_iSpec = np.copy(integral_spec).reshape(np.prod(integral_spec.shape))
      del integral_spec
      integral_spec = copy_iSpec
      # we compute a definite integral: g(z) = int_{z'=z_bottom}^{z'= z} dz' f(z')
      # hence g(z_bottom=0)
      integral_spec[0] = 0.
      for i in range(integral_spec.shape[0]-1):
         integral_spec[0] += integral_spec[i+1] * (-1.)**(i)
      self.verticalSums.append(np.sum(integral_spec)/4./self.NZ)
      self.lut_verticalSums.append('Int dz '+self.lut_profiles[pos_in_list])
      integral_spec[0]*=2.
      integral_phys = fft.idct(integral_spec)/self.NZ/8.
      self.profiles.append(integral_phys)
      self.lut_profiles.append('Int dz '+self.lut_profiles[pos_in_list])

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
      self.lut_profiles.append('<'+self.lut[pos_in_list]+ '>_{x,y}')

