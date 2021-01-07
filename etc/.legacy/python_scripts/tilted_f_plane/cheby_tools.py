import numpy as np
import scipy.sparse as sp

def chebyshev_elementary_multiplication(N, gap, center):
   sub__diag = np.ones((1,N), dtype=np.float_)*gap/4.
   main_diag = np.ones((1,N), dtype=np.float_)*center
   superdiag = np.ones((1,N), dtype=np.float_)*gap/4.
   sub__diag[0,0]=gap/2.
   matrix_entries = np.concatenate((sub__diag, main_diag, superdiag), axis = 0)
   offsets = np.array([-1,0,+1])
   cheby_M_mat = sp.dia_matrix((matrix_entries, offsets),shape=(N,N),dtype=np.float_)
   return cheby_M_mat

def chebyshev_elementary_integration(N, gap, center):
   sub__diag = 1./(np.arange(N).reshape(1,N) + 1.)
   sub__diag[0,0]=2.
   superdiag = (np.arange(N).reshape(1,N) - 1.)
   superdiag[0,0:1] = 1.
   superdiag = -1./superdiag
   superdiag[0,0:1] = 0.
   matrix_entries = np.concatenate((sub__diag, superdiag), axis = 0)*gap/2
   offsets = np.array([-1,1])
   cheby_I_mat = sp.dia_matrix((matrix_entries, offsets),shape=(N,N),dtype=np.float_)
   return cheby_I_mat


