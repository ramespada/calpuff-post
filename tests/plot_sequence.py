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

proyecto="test"
# Path to calpuff output file:
#file='./data/conc_losangeles.dat'
#file='./data/conc_ducson.dat'
#file='./data/conc_waterloo.dat'
file='./data/conc_mumbai.dat'

#Read file:
puff=calpost.read_file(file)

X,Y = puff.get_coordinates()

pollut=puff.species[0]
unit=puff.units[0]
periods=1

xmin=np.min(X)
xmax=np.max(X)
ymin=np.min(Y)
ymax=np.max(Y)

projection = ccrs.UTM(zone=puff.proj["zone"], southern_hemisphere=(puff.proj['hemis']=='S'))
request=cimgt.GoogleTiles(style='satellite')

C=puff.get_data(pollut)

if ( puff.ndrec > puff.ngrec ): # if discrete receptor is
    # Create a regular grid covering the domain of X and Y
    xi = np.arange(puff.x0, np.max(X), puff.dx)
    yi = np.arange(puff.y0, np.max(Y), puff.dy)
    X_grid, Y_grid = np.meshgrid(xi, yi)
    C_grid=np.zeros([puff.run_length, len(xi), len(yi)])
    points = np.column_stack((X, Y))

    # Perform bilinear (linear) interpolation
    for t in range(puff.run_length):
        C_grid[t,:,:] = griddata(points, C[t,:], (X_grid, Y_grid), method='linear')

    X=X_grid; Y=Y_grid
    C=C_grid

#Check units so we get ug/m3
if ( puff.units[0] == 'g/m3' ):
    C=C*1e6 #g/m3 to ug/m3
    unit ='ug/m3'

nlev=10
borde=-20
extent=[xmin-borde,xmax+borde,ymin-borde,ymax+borde]

times = list(range(0,puff.run_length,3))
rows=2
cols=len(times)//rows

print("cols:",cols,"rows:",rows)

nlevs=10
minval=np.quantile(C[C>0],0.2)
maxval=np.nanmax(C)*1.1
levels = gen_levels(minval,maxval,nlevs)

fig,axes = plt.subplots(rows, cols, figsize=(12, 6), sharex=True, sharey=True, gridspec_kw = {'wspace':0.2, 'hspace':0.007}, subplot_kw={'projection': projection, "aspect": 2})

for t in range(len(times)):

    col=t%cols; row=t//cols
    print("col=",col,"row=",row)
    ax = axes[row,col]
    ax.set_extent(extent, crs=projection)
    #ax.add_image(request,13)
    ax.add_image(request,9)
    ax.tick_params(labelsize=8)

    if ( col == 0      ):
        ax.set_yticks(np.linspace(ymin,ymax,4).round(0), crs=projection)

    if ( row == rows-1 ):
        ax.set_xticks(np.linspace(xmin,xmax,3).round(0), crs=projection)

    pf=ax.contourf(X,Y,C[t,:,:],alpha=0.6,cmap="RdBu_r", levels=levels, transform=projection , norm = LogNorm() )  
    cbar=fig.colorbar(pf, ax=ax, shrink=0.5)

    #p=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap="RdBu_r",levels=levels,norm = LogNorm(),transform=projection)#,cmap='RdBu_r',levels=levels    

    cbar.ax.set_title("Conc. \n ["+unit+"]",loc="left",  fontsize=10)
    ax.set_title(pollut+" "+str( times[t]) ) #,fontsize=6)

#plt.show()
#plt.tight_layout()
plt.savefig('conc_'+pollut+'.svg',dpi=1400, bbox_inches = 'tight',pad_inches = 0.1)

