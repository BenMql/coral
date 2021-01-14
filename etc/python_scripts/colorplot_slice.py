
import numpy as np
import matplotlib.pyplot as plt
from get_coral_params import domDecomp, coralGeom
from sys import argv


dd = domDecomp  (0)
geom = coralGeom(0)


NXAA = dd.NXAA
NYAA = dd.NYAA
NZAA = dd.NZAA

z = geom.center+ 0.5*geom.gap*np.cos ((
       np.linspace(0,NZAA, num=NZAA, endpoint=False)*2+1.)*np.pi/2./NZAA)
x = np.linspace(0,geom.Lx, NXAA)
y = np.linspace(0,geom.Ly, NYAA)

time = 10000
var = 4 # this is temperature

kindOfVar = input('Linear (l) or Quadratic (q) variable? [l/q] >> ')
varNum    = input('Var. number? [integer]                      >> ')
posStr    = input('Slice position? [eg x00001, y00512, z00064] >> ')
timeStr   = input('Time? [integer]                             >> ')

varDict = {"l" : "linear", "q" : "quadratic",
           "L" : "linear", "Q" : "quadratic"}

reso1 = {'x': NYAA, 'y': NXAA, 'z' : NXAA}
reso2 = {'x': NZAA, 'y': NZAA, 'z' : NYAA}
grid1 = {'x': y, 'y': x, 'z' : x}
grid2 = {'x': z, 'y': z, 'z' : y}


a = np.fromfile('Slices/'+varDict[kindOfVar]+'_var'+varNum.zfill(2)+'_'+posStr+'_time'+timeStr.zfill(8)+'_full.dat', dtype=np.float_).reshape(reso1[posStr[0]], reso2[posStr[0]])


plt.figure()
try:
  plt.pcolormesh(grid1[posStr[0]], grid2[posStr[0]],a.T, cmap=argv[1])
except IndexError:
  plt.pcolormesh(grid1[posStr[0]], grid2[posStr[0]],a.T)
plt.colorbar()
#plt.title('temperature fluctuations, x-z slice at first Y-gridpoint')
plt.show()

