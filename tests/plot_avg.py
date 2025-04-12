#!/usr/bin/env python3

import math
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.io.img_tiles as cimgt
from osgeo import gdal
import geopandas as gpd

import calpost
from levels import gen_levels

# Path to calpuff output file:
file='./data/conc_losangeles.dat'
#file='./data/conc_ducson.dat'
#file='./data/conc_waterloo.dat'
#file='./data/conc_mumbai.dat'

#Read file:
puff=calpost.read_file(file)

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
    
    C=puff.get_time_avg_max(pollutid,period)
    C=C*1e6
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    
    ax.set_extent(extent, crs=projection)
    ax.set_xticks(np.linspace(xmin,xmax,3).round(0), crs=projection)
    ax.set_yticks(np.linspace(ymin,ymax,4).round(0), crs=projection)
    ax.tick_params(labelsize=8)
    ax.add_image(request, 12)
    
    minval=max(1e-9, np.nanmin(C))
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

