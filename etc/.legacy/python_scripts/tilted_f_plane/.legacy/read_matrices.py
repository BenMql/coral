import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib.pyplot as plt

class coral_matrix:
   def __init__(self, path_n_prefix, type_char):
      type_dict = {'d': np.float64,
                   'z': np.complex128 }
      self.dtype= type_dict[type_char]
      my_dat  = np.fromfile(path_n_prefix+'.dat',     dtype=self.dtype)
      indices = np.fromfile(path_n_prefix+'.indices', dtype=np.int32)
      nelems = indices[0]
      nrow = indices[1]
      ncol = indices[2]
      rows = indices[3:3+nrow+1]
      cols = indices[3+nrow+1:3+nrow+1+nelems]
      self.csr = sp.csr_matrix((my_dat,cols-1, rows-1), shape=(nrow,ncol))
       

   def compute_full_spectrum(self, other=None):
      my_L =  self.csr.toarray()
      if other is not None:
         my_M = other.csr.toarray()
      else:
         M=None
      self.eVal, self.eVec = la.eig(my_L,my_M, check_finite = False)

   def plot_full_spectrum(self, other=None):
      if not hasattr(self, 'eVal'):
          self.compute_full_spectrum(other=other)
      plt.figure()
      plt.plot(self.eVal.real, self.eVal.imag, 'o')
      plt.xlabel('Real part')
      plt.ylabel('Imaginary part')
      plt.title ('Spectrum in the complex plane')
     

  
