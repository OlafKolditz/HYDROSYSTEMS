import matplotlib.pyplot as plt
import numpy as np
x,y,h = np.loadtxt("out.txt", unpack=True)
plt.contour(x.reshape(33,43), y.reshape(33,43), h.reshape(33,43))
#plt.tricontour(x,y,z)
plt.show()