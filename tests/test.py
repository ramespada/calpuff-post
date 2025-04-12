#!/usr/bin/python3

import numpy as np
from scipy.interpolate import griddata
from matplotlib import pyplot as plt
import calpost

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
units=dat.units[0]

# Get discrete receptor concentrations.
X,Y   = dat.get_coordinates()

## Get gridded concentrations.
C     = dat.get_data(pollut)

## Compute time-averaged concentrations for specific pollutant.
C_1hr = dat.get_time_avg_max(pollut, interval=24, rank=1)

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

    C_grid_1hr=griddata(points, C_1hr, (X_grid, Y_grid), method='linear')

    X=X_grid; Y=Y_grid 
    C=C_grid
    C_1hr=C_grid_1hr

#Check units so we get ug/m3
if ( dat.units[0] == 'g/m3' ):
    C_1hr=C_1hr*1e6 #g/m3 to ug/m3

levels=[0.02,0.03,0.05,0.08,0.10,0.3,0.5,0.8,1.0,1.86]

pf=plt.contourf(X,Y,C_1hr, alpha=0.6,cmap="Spectral_r",levels=levels) #,norm = LogNorm())
cbar=plt.colorbar(pf, shrink=0.7)
pc=plt.contour(X,Y,C_1hr, alpha=0.8, colors='black', linewidths=0.5, levels=levels)
plt.clabel(pc, fmt='%.2f', fontsize=6)
plt.show()

