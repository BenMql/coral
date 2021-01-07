import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.linalg as la

indices = np.fromfile('../../../build/mass.core0001.indices', dtype=np.int32)
nelems = indices[0]
nrow = indices[1]
ncol = indices[2]
row  = indices[3:3+nrow+1]-1
col  = indices[3+nrow+1:]-1
data = np.fromfile('../../../build/mass.core0001.dat', dtype=np.complex_)
mass = sp.csr_matrix((data, col, row), shape=(nrow, ncol))

indices = np.fromfile('../../../build/stif.core0001.indices', dtype=np.int32)
nelems = indices[0]
nrow = indices[1]
ncol = indices[2]
row  = indices[3:3+nrow+1]-1
col  = indices[3+nrow+1:]-1
data = np.fromfile('../../../build/stif.core0001.dat', dtype=np.complex_)
stif = sp.csr_matrix((data, col, row), shape=(nrow, ncol))

stifF = stif.todense()
massF = mass.todense()

val,vec = la.eig(stifF, b=massF)

plt.plot(val.real, val.imag, 'o')
plt.show()


