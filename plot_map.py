#!/usr/bin/env python3

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.io.img_tiles as cimgt
from osgeo import gdal
import geopandas as gpd
from matplotlib.colors import LogNorm

import calpost

#polluts=("CO", "CO", "NOx", "NOx", "SO2", "SO2", "PM10", "PM10", "PM25", "PM25") #, "H2S")
#periods=("1", "8", "1", "ANNUAL", "1", "24", "24", "ANNUAL", "24", "ANNUAL") #, "8")


puff=calpost.read_file('conc.dat')

X,Y=puff.get_coordinates()

xmin=np.min(X) 
xmax=np.max(X)
ymin=np.min(Y)
ymax=np.max(Y) 

print(puff.proj['hemis']=='S')
projection = ccrs.UTM(zone=puff.proj["zone"], southern_hemisphere=(puff.proj['hemis']=='S'))
request=cimgt.GoogleTiles(style='satellite')

borde=-20
extent=[xmin-borde,xmax+borde,ymin-borde,ymax+borde]

polluts=puff.species
periods=[1]
for i in range(len(polluts)):
    pollutid=polluts[i]
    period=periods[i]
    unit=puff.units[i] #"ug/m3"
    
    C=puff.time_avg_max(pollutid,period)
    C=C*1e6
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    
    ax.set_extent(extent, crs=projection)
    ax.set_xticks(np.linspace(xmin,xmax,3).round(0), crs=projection)
    ax.set_yticks(np.linspace(ymin,ymax,4).round(0), crs=projection)
    ax.tick_params(labelsize=8)
    ax.add_image(request,12)
    
    minval=0.0
    maxval=1.0*np.nanmax(C)
    #levels=np.array([0,0.001,0.005,0.01,0.02,0.1,0.15,0.2,0.5,0.7,0.9,1,3,5,7,10,15,20,50,80,90,100,150,200])#, 25, 50, 100]) ç
    #levels=np.linspace(minval,maxval,9)
    levels=np.array([0.7,1.0,2.0,5.0,7.0,10.0,20.0,50.0,70.0,72.1])


    #p=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap='viridis',levels=levels)    
    pf=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap="Spectral_r",levels=levels)#,cmap='RdBu_r',levels=levels,norm = LogNorm())    
    cbar=fig.colorbar(pf, shrink=0.7)

    pc=ax.contour(X,Y,C[:,:],alpha=0.6,colors='black', linewidths=0.5, levels=levels)
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

