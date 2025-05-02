import numpy as np
#import matplotlib as plt

import matplotlib.pyplot as plt
from numpy import arange
from numpy import meshgrid

#
# Plot for paper
#

a, b = 1, 0

bnd1 = lambda x, y : (a + x)**2  + (b + y)**2 - (a - x)**2 - (b - y)**2
bnd2 = lambda x, y : (((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**2).real - abs((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**4


delta = 0.005
xvals = np.arange(-2, 2, delta)
yvals = np.arange(-2.0, 2.0, delta)
x, y = np.meshgrid(xvals, yvals)

fig, ax = plt.subplots()

ax.imshow( ((bnd1(x,y)<=0) & (bnd2(x,y)<=0)).astype(int), extent=(x.min(),x.max(),y.min(),y.max()), origin="lower", cmap="Greys", alpha=0.4)
ax.contour(x, y, bnd1(x,y), [0])
ax.contour(x, y, bnd2(x,y), [0], linestyles = 'dashed')
#ax.plot(-a, -b, 'o', color = 'black')


#bnd1 = (a + X)**2  + (b + Y)**2 - (a - X)**2 - (b - Y)**2
#bnd2 = (((a + X + (b + Y)*1j)/(X - a + (Y - b)*1j))**2).real - abs((a + X + (b + Y)*1j)/(X - a + (Y - b)*1j))**4

#plt.contour(X, Y, bnd1, [0])
#plt.contour(X, Y, bnd2, [0])

#fig, ax = plt.subplots()
#ax.contour(X, Y, bnd1, [0])
#ax.contour(X, Y, bnd2, [0])
#ax.fill_between(X, bnd1, bnd2, alpha=0.2)

plt.show()