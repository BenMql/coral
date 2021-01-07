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
kx = 81.681408993334628  
ky =  4.1887902047863905     

plt.plot(-kx**2 -ky**2 - (np.arange(1,64) *np.pi)**2, np.zeros((63)), 'x')
vval = np.sort(val.real)
errRel = np.zeros((60))
for iz in range(60):
   errRel[iz] = vval[-1-iz]/(kx**2 + ky**2 + (iz+1)**2 * np.pi**2 ) + 1 

plt.figure()
plt.semilogy(np.abs(errRel))
plt.show()


