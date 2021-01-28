
import numpy as np
import matplotlib.pyplot as plt
from get_coral_params import domDecomp, coralGeom
from sys import argv, version_info



dd = domDecomp  (0)
geom = coralGeom(0)


NXAA = dd.NXAA
NYAA = dd.NYAA
NZAA = dd.NZAA

z = geom.center+ 0.5*geom.gap*np.cos ((
       np.linspace(0,NZAA, num=NZAA, endpoint=False)*2+1.)*np.pi/2./NZAA)
x = np.linspace(0,geom.Lx, NXAA)
y = np.linspace(0,geom.Ly, NYAA)

reso1 = {'x': NYAA, 'y': NXAA, 'z' : NXAA}
reso2 = {'x': NZAA, 'y': NZAA, 'z' : NYAA}
grid1 = {'x': y, 'y': x, 'z' : x}
grid2 = {'x': z, 'y': z, 'z' : y}
varDict = {"l" : "linear", "q" : "quadratic",
           "L" : "linear", "Q" : "quadratic"}

posStr = 'y00001'
if version_info[0]<3:
  timeStr   = raw_input('Time? [integer]  >> ')
else:
  timeStr   = input('Time? [integer]  >> ')

kindOfVar='l'


varNum   ='1'
u = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

varNum   ='3'
w = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

varNum   ='2'
v = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

for ix in range(v.shape[0]):
   v[ix,:] += -239.*z + 3900./z


plt.figure()
try:
  plt.pcolormesh(grid2[posStr[0]], grid1[posStr[0]],v, cmap=argv[1])
except IndexError:
  plt.pcolormesh(grid2[posStr[0]], grid1[posStr[0]],v)
plt.colorbar()
z_equi = np.linspace(z[-1], z[0], num=200)
u_equi = np.zeros((u.shape[0], z_equi.shape[0]), dtype=np.float_)
w_equi = np.zeros((u.shape[0], z_equi.shape[0]), dtype=np.float_)
for ix in range(u_equi.shape[0]):
   u_equi[ix,:]=np.interp(z_equi, z[::-1], u[ix,::-1])
   w_equi[ix,:]=np.interp(z_equi, z[::-1], w[ix,::-1])
plt.streamplot(z_equi, x, u_equi, w_equi, color='k')
#plt.streamplot(w[:,::-1].T, u[:,::-1].T, color='k')
#plt.title('temperature fluctuations, x-z slice at first Y-gridpoint')
#
#
#############################
#
#############################
#
#############################
#

posStr = 'z00032'

varNum   ='3'
w = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

varNum   ='2'
v = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

varNum   ='1'
u = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])

plt.figure()
try:
  plt.pcolormesh(grid2[posStr[0]], grid1[posStr[0]],v, cmap=argv[1])
except IndexError:
  plt.pcolormesh(grid2[posStr[0]], grid1[posStr[0]],v)
plt.colorbar()
plt.streamplot(grid2[posStr[0]], grid1[posStr[0]],v, w, color='k')

plt.show()
