#!/usr/bin/env python3

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.io.img_tiles as cimgt
from osgeo import gdal
import geopandas as gpd
from matplotlib.colors import LogNorm

import calpost

import math

def gen_levels(minval, maxval, nlevels):
    if minval <= 0 or maxval <= 0:
        raise ValueError("Both minval and maxval must be positive numbers")
    if minval >= maxval:
        raise ValueError("minval must be less than maxval")
    if nlevels < 1:
        raise ValueError("nlevels must be at least 1")

    # Calculate logarithmic range
    log_min = math.log10(minval)
    log_max = math.log10(maxval)

    # Generate initial candidates in log space
    log_candidates = []

    # Always include the endpoints
    log_candidates.append(log_min)
    log_candidates.append(log_max)

    # Generate key numbers (1 and 5 in each decade)
    current_decade = math.floor(log_min)
    while current_decade <= math.ceil(log_max):
        # Add 1 and 5 in this decade
        log1 = current_decade
        log5 = current_decade + math.log10(5)

        if log1 > log_min and log1 < log_max:
            log_candidates.append(log1)
        if log5 > log_min and log5 < log_max:
            log_candidates.append(log5)

        current_decade += 1

    # Add evenly spaced levels in log space (if needed to reach nlevels)
    if len(log_candidates) < nlevels + 1:
        additional_points = nlevels + 1 - len(log_candidates)
        step = (log_max - log_min) / (additional_points + 1)
        for i in range(1, additional_points + 1):
            log_point = log_min + i * step
            if log_point not in log_candidates:
                log_candidates.append(log_point)

    # Convert back to linear space and round to proper form
    candidates = [10**x for x in log_candidates]

    # Round to have only one non-zero digit
    def round_to_single_digit(x):
        if x == 0:
            return 0
        magnitude = 10 ** math.floor(math.log10(x))
        first_digit = round(x / magnitude)
        return first_digit * magnitude

    rounded = [round_to_single_digit(x) for x in candidates]

    # Remove duplicates and sort
    unique_levels = sorted(list(set(rounded)))

    # Ensure minval and maxval are exactly represented
    unique_levels[0] = minval
    unique_levels[-1] = maxval

    # If we still don't have enough levels, add more while maintaining constraints
    while len(unique_levels) < nlevels + 1:
        # Find the largest gap
        max_gap = 0
        gap_index = 0
        for i in range(len(unique_levels) - 1):
            gap = unique_levels[i+1] - unique_levels[i]
            if gap > max_gap:
                max_gap = gap
                gap_index = i

        # Add midpoint in log space
        log_mid = (math.log10(unique_levels[gap_index]) +
                   math.log10(unique_levels[gap_index+1])) / 2
        new_val = round_to_single_digit(10**log_mid)
        unique_levels.append(new_val)
        unique_levels = sorted(list(set(unique_levels)))

    return sorted(unique_levels)

#puff=calpost.read_file('conc_losangeles.dat')
#puff=calpost.read_file('conc_ducson.dat')
puff=calpost.read_file('conc_waterloo.dat')
#puff=calpost.read_file('conc_mumbai.dat')

X,Y = puff.get_coordinates()

xmin=np.min(X) 
xmax=np.max(X)
ymin=np.min(Y)
ymax=np.max(Y) 

projection = ccrs.UTM(zone=puff.proj["zone"], southern_hemisphere=(puff.proj['hemis']=='S'))
request=cimgt.GoogleTiles(style='satellite')

nlev=10
borde=-20
extent=[xmin-borde,xmax+borde,ymin-borde,ymax+borde]

polluts=puff.species
periods=[1]
for i in range(len(polluts)):
    pollutid=polluts[i]
    period=periods[i]
    unit=puff.units[i] #"ug/m3"
    
    C=puff.time_avg_max(pollutid,period)
    C=C#*1e6
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    
    ax.set_extent(extent, crs=projection)
    ax.set_xticks(np.linspace(xmin,xmax,3).round(0), crs=projection)
    ax.set_yticks(np.linspace(ymin,ymax,4).round(0), crs=projection)
    ax.tick_params(labelsize=8)
    ax.add_image(request,12)
    
    minval=max(1e-3, np.nanmin(C))
    maxval=1.0*np.nanmax(C)
    levels=gen_levels(minval,maxval,nlev)

    pf=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap="Spectral_r",levels=levels,norm = LogNorm())
    cbar=fig.colorbar(pf, shrink=0.7)

    pc=ax.contour(X,Y,C[:,:],alpha=0.8,colors='black', linewidths=0.5, levels=levels)
    ax.clabel(pc, fmt='%.2f', fontsize=6)

    cbar.ax.set_title(pollutid+"_"+str(period)+"\n("+unit+")")
    ax.set_title("", fontsize=9)
    
    #dominio.plot(ax=ax, color='none',edgecolor='red',alpha=0.95)
    #predio.plot(ax=ax, color='#f0f0f0',edgecolor='black',alpha=0.75)
    #edif.plot(ax=ax, color='#fafcc2',edgecolor='black',alpha=0.5)
    #calles.plot(ax=ax,color='#ebe7ca',linewidth=0.5)
    
    fig.set_size_inches(10,6)
    plt.savefig('plot_'+pollutid[0:12].strip()+'_'+str(period)+'.svg',dpi=1400,
                bbox_inches = 'tight',pad_inches = 0.1)

