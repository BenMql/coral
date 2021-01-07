import numpy as np

Tmean = np.mean(T,axis=(0,1))
d1, d2, Tm3d = np.meshgrid(np.arange(T.shape[0]), np.arange(T.shape[1]), Tmean)

toto = T - Tm3d
Tfluc = toto.astype(np.float32)
