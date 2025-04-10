#!/usr/bin/python3

import calpost
from matplotlib import pyplot as plt

# Open and parse a CALPUFF output file
calpuff_dat = calpost.read_file("conc.dat")

# Print access parsed information
calpuff_dat.info()           

# Get gridded concentrations.
X,Y = calpuff_dat.get_coordinates()
C=calpuff_dat.get_data("SO2           1")

plt.contourf(X,Y,C[4,:,:])
plt.show()

## Ask to compute time-averaged concentrations for specific pollutant.
#
C_1hr = calpuff_dat.time_avg_max(calpuff_dat.species[0], 1)
 
plt.imshow(C_1hr)
plt.show()

## Get ranked values table for a given pollutant:
#
 
 
