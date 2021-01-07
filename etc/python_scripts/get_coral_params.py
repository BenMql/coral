import numpy as np

class domDecomp():
   def __init__(self, iCore):
       # read a List of Integers (LoI)
       LoI= np.fromfile('Geometry/domDecmp.core'+str(iCore).zfill(4), dtype=np.int32)
       # sort that list into a class
       self.NX = LoI[0]
       self.NY = LoI[1]
       self.NZ = LoI[2]
       self.NXAA = LoI[3]
       self.NYAA = LoI[4]
       self.NZAA = LoI[5]
       self.spec_local_NX = LoI[6]
       self.spec_local_NY = LoI[7]
       self.nVar          = LoI[8]
       self.phys_iStart   = LoI[9]
       self.phys_iSize    = LoI[10]
       self.phys_iEnd     = LoI[11]
       self.spec_iStart   = LoI[12]
       self.spec_iSize    = LoI[13]
       self.spec_iEnd     = LoI[14]

class coralGeom():
   def __init__(self, iCore):
       # read a List of Floats64
       LoF= np.fromfile('Geometry/geometry.core'+str(iCore).zfill(4), dtype=np.float64)
       # sort that list into a class
       self.Lx  = LoF[0]
       self.Ly  = LoF[1]
       self.gap = LoF[2]
       self.center = LoF[3]
       self.kx_box = LoF[4]
       self.ky_box = LoF[5]
       self.spec_local_NX = np.fromfile('Geometry/domDecmp.core'+str(iCore).zfill(4), dtype=np.int32)[6]
       self.spec_local_NY = np.fromfile('Geometry/domDecmp.core'+str(iCore).zfill(4), dtype=np.int32)[7]
       self.kx = LoF[6:6+self.spec_local_NX]
       self.ky = LoF[6+self.spec_local_NX:6+self.spec_local_NX+self.spec_local_NY]


