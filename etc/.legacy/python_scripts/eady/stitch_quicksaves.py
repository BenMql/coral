import numpy as np
from os import walk

path2data = './Restart/'

full_loq = []
for (dirpath, dirnames, filenames) in walk(path2data):
   full_loq.extend(filenames)

loq = [s for s in full_loq if ('core' in s)]

loq.sort(key=str.lower)

#first remove dt and atime

for fname in loq:
  a = np.fromfile(path2data+fname, dtype=np.float_)
  a[:-2].tofile(path2data+fname)
a[-2:].tofile(path2data+'timing_details.dat')

with open(path2data+'monolithic_quicksave.dat', 'w') as outfile:
   for fname in loq:
      with open(path2data+fname) as infile:
         outfile.write(infile.read())

