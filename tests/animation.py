#!/usr/bin/python3

import math
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import animation
import calpost

from levels import gen_levels

# File path to calpuff output file
#file='./data/conc_mumbai.dat'
#file='./data/conc_losangeles.dat'
file='./data/conc_ducson.dat'
#file='./data/conc_waterloo.dat'

# Open and parse a CALPUFF output file
dat = calpost.read_file(file)

# Print access parsed information
dat.info()

# Get specie name
pollut=dat.species[0]

# Get discrete receptor concentrations.
X,Y   = dat.get_coordinates()

## Get gridded concentrations.
C     = dat.get_data(pollut)

#If not gridded receptors, then interpolate to regular grid.
if ( dat.ndrec > dat.ngrec ): # if discrete receptor is
    # Create a regular grid covering the domain of X and Y
    xi = np.arange(dat.x0, np.max(X), dat.dx)
    yi = np.arange(dat.y0, np.max(Y), dat.dy)
    X_grid, Y_grid = np.meshgrid(xi, yi)
    C_grid=np.zeros([dat.run_length, len(xi), len(yi)])
    points = np.column_stack((X, Y))

    # Perform bilinear (linear) interpolation
    for t in range(dat.run_length):
        C_grid[t,:,:] = griddata(points, C[t,:], (X_grid, Y_grid), method='linear')

    X=X_grid; Y=Y_grid
    C=C_grid*1e6

nlev=10
minval=max(1e-12, np.nanmin(C))
maxval=np.nanmax(C)

levels=gen_levels(minval,maxval,nlev)

fig, ax = plt.subplots()
contourf = [ax.contourf(X, Y, C[0,:,:], levels=levels, cmap='viridis',norm = LogNorm())]
cbar = plt.colorbar(contourf[0], ax=ax)

def update(t):
    for coll in contourf[0].collections:
        coll.remove()  # Remove old contours
    contourf[0] = ax.contourf(X, Y, C[t,:,:], levels=levels, cmap='viridis',norm = LogNorm())
    ax.set_title(f"Frame {t}")
    return contourf[0].collections

ani = animation.FuncAnimation(
    fig, update, frames=C.shape[0], blit=False
)

# Save as GIF using Pillow
ani.save("animation_"+pollut+".gif", writer='pillow', fps=5)

plt.close()
print("GIF saved as animation.gif")

