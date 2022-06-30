import matplotlib.pyplot as plt
import numpy as np
x,y,h = np.loadtxt("out.txt", unpack=True)

CS=plt.contour(x.reshape(11,21), y.reshape(11,21), h.reshape(11,21))
#plt.tricontour(x,y,z)
plt.clabel(CS, CS.levels, inline=True, fmt='%1.3f', fontsize=10)
plt.savefig('tend.png')
plt.show()