#!/usr/bin/python3

from matplotlib import pyplot as plt
import calpost

# Open and parse a CALPUFF output file
#calpuff_dat = calpost.read_file('./data/conc_mumbai.dat')
#calpuff_dat = calpost.read_file('./data/conc_losangeles.dat')
calpuff_dat = calpost.read_file('./data/conc_ducson.dat')
#calpuff_dat = calpost.read_file('./data/conc_waterloo.dat')

# Print access parsed information
calpuff_dat.info()           

print(calpuff_dat.model_version    )

pollut=calpuff_dat.species[0]

# Get discrete receptor concentrations.
X = calpuff_dat.x_r
Y = calpuff_dat.y_r
C = calpuff_dat.get_discrete_data(pollut)

plt.scatter(X, Y, c=C[1,:],         # Color values
    cmap='viridis',                 # Colormap (options: 'viridis', 'plasma', 'cool', 'jet', etc.)
    alpha=0.8,                      # Transparency (0=transparent, 1=opaque)
    edgecolors='k',                 # Edge color of markers ('k'=black)
    s=100                           # Marker size
)
plt.show()

## Get gridded concentrations.
#X,Y = calpuff_dat.get_coordinates()
#C   = calpuff_dat.get_data(pollut)
#
#plt.contourf(X,Y,C[4,:,:])
#plt.show()
#
#### Ask to compute time-averaged concentrations for specific pollutant.
##
#C_1hr = calpuff_dat.time_avg_max(calpuff_dat.species[0], 1)
# 
#plt.imshow(C_1hr)
#plt.show()
#
### Get ranked values table for a given pollutant:
# 
 
