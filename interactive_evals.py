import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Configure Matplotlib to use LaTex for text rendering
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Computer Modern"

# Define two constraint functions that determine the admissible region for the second eigenvalue
# The real and imaginary parts of the second eigenvalue are x and y, respectively
bnd1 = lambda a, b, x, y : (a + x)**2  + (b + y)**2 - (a - x)**2 - (b - y)**2
bnd2 = lambda a, b, x, y : (((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**2).real - abs((a + x + (b + y)*1j)/(x - a + (y - b)*1j))**4

# Set up the plotting grid
delta = 0.005
x_min, x_max = -2, 2
y_min, y_max = -2.0, 2.0
xvals = np.arange(x_min, x_max, delta)
yvals = np.arange(y_min, y_max, delta)
x, y = np.meshgrid(xvals, yvals)

# Create the main figure and axes for the plot
fig, ax = plt.subplots(figsize=(10,10))
# Adjust the main plot to make room for the sliders
fig.subplots_adjust(bottom=0.25)

# Define the position of the real part slider axis: [left, bottom, width, height]
ax_real = plt.axes([0.25, 0.1, 0.6, 0.03])
# Make a horizontal slider to control real part.
real_slider = Slider(
    ax = ax_real,
    label = 'Real part',
    valmin = -5,
    valmax = 5,
    valinit = 1,
)

# Define the position of the imaginary part slider axis: [left, bottom, width, height]
ax_imag = plt.axes([0.25, 0.05, 0.6, 0.03])
# Make a horizontal slider to control imaginary part.
imag_slider = Slider(
    ax = ax_imag,
    label = 'Imaginary part',
    valmin = -5,
    valmax = 5,
    valinit = 0,
)

# The function to be called anytime a slider's value changes
def update(val):
    a = real_slider.val
    b = imag_slider.val

    ax.clear()

    # Re-set all axis properties (since ax.clear() removes them)
    ax.set_xlabel("Real")
    ax.set_ylabel("Imaginary")
    ax.set_title('Admissible Eigenvalues')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Avoid division by 0
    if a !=0 or b!=0:

        ax.imshow(((bnd1(a,b,x,y)<=0) & (bnd2(a,b,x,y)<=0)).astype(int), 
                extent=(x.min(),x.max(),y.min(),y.max()), 
                origin="lower", 
                cmap="Greys", 
                alpha=0.4
        )

        ax.contour(x, y, bnd1(a, b, x, y), [0], colors='black')
        ax.contour(x, y, bnd2(a, b, x, y), [0], linestyles = 'dashed', colors='black')

        if x_min <= -a <= x_max and y_min <= -b <= y_max:
            ax.plot([-a], [-b], "o", markersize=5, color='black')
        else:
            pass
    else:
        pass

    fig.canvas.draw_idle()

# register the update function with each slider
real_slider.on_changed(update)
imag_slider.on_changed(update)

# Initial call to update to draw the plot for the first time with initial values
update(None)

# Display the plot
plt.show()