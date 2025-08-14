import numpy as np
import matplotlib.pyplot as plt
from numpy import arange
from numpy import meshgrid

#
# Code for Figure 1
#

# Configure Matplotlib to use LaTex for text rendering
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Computer Modern"

# Define the real (a) and imaginary (b) parts of a known eiegenvalue
a, b = 1, 0

# Define two constraint functions that determine the admissible region for the second eigenvalue
# The real and imaginary parts of the second eigenvalue are x and y, respectively
bnd1 = lambda x, y : (a + x)**2  + (b + y)**2 - (a - x)**2 - (b - y)**2
bnd2 = lambda x, y : (((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**2).real - abs((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**4

# Set up the plotting grid
delta = 0.005
xvals = np.arange(-2, 2, delta)
yvals = np.arange(-2.0, 2.0, delta)
x, y = np.meshgrid(xvals, yvals)

# Use imshow to plot the region defined by the constraints
plt.imshow( ((bnd1(x,y)<=0) & (bnd2(x,y)<=0)).astype(int), extent=(x.min(),x.max(),y.min(),y.max()), origin="lower", cmap="Greys", alpha=0.4)
# Plot boundaries of admissible region
plt.contour(x, y, bnd1(x,y), [0], colors='black')
plt.contour(x, y, bnd2(x,y), [0], linestyles = 'dashed', colors='black')

# Plot the point corresponding to -a - bi since this is admissible
plt.plot(-a, -b, "o", markersize=5, color = 'black')

# Label the axes
plt.xlabel("Real")
plt.ylabel("Imaginary")
plt.grid(True, linestyle=':', alpha=0.6)

# Display the final plot
plt.show()