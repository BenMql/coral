import numpy as np
import matplotlib.pyplot as plt

err111 = [2.159e-6,8.578e-6,3.383e-5,0.0002027, 0.0007561, 0.002631, 0.01088]
err222 = [1.153e-10, 9.161e-10, 7.234e-9,
          1.087e-7,  8.144e-7,  5.717e-6,
          6.037e-5
          ]
err443 = [3.85e-14, 6.246e-13, 9.817e-12,
          3.681e-10, 5.499e-9, 7.671e-8,
          1.986e-6
          ]
dt = [1.e-5, 2.e-5, 4.e-5,
      1.e-4, 2.e-4, 4.e-4,
      1.e-3]

plt.figure()
plt.loglog(dt, err111, 'vk', markerfacecolor='xkcd:blue', label='ARS111')
plt.loglog([1.e-5, 1.e-3], [3.e-6, 3.e-2], 'k--')
plt.loglog(dt, err222, 'vk', markerfacecolor='xkcd:sage', label='ARS222')
plt.loglog([1.e-5, 1.e-3], [3.e-10, 3.e-4], 'k--')
plt.loglog(dt, err443, 'vk', markerfacecolor='xkcd:salmon', label='ARS443')
plt.loglog([1.e-5, 1.e-3], [5.e-14, 5.e-6], 'k--')
plt.grid(True)
plt.grid(True, which='minor', linestyle='dashed')
plt.xlabel(r'timestep $\delta t$', fontsize=16)
plt.title(r'$L_2$ of error after 10 steps', fontsize=18)
plt.legend(fontsize=16)
plt.show()
