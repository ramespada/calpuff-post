#!/usr/bin/env python3

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import cartopy.io.img_tiles as cimgt
import geopandas as gpd
from osgeo import gdal
from matplotlib.colors import LogNorm

def getData(filepath):
    ds=gdal.Open(filepath, gdal.GA_ReadOnly)
    for x in range(1, ds.RasterCount + 1):
        #band = ds.GetRasterBand(x)
        C = ds.ReadAsArray()
    nx=ds.RasterXSize
    ny=ds.RasterYSize
    ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
    lrx = ulx + (nx * xres)
    lry = uly + (ny * yres)
    C[C<=0]=np.NaN
    C=np.nan_to_num(C)
    X,Y=np.meshgrid(np.linspace(ulx,lrx,nx),np.linspace(uly,lry,ny))
    return X,Y,C

directory="/home/ram/github/gauss-ma/aermod/runs/45-PAE-25"
proyecto="CPAE"
epsg="32720"
utmzone=20

dominioPath=directory+"/GIS/dominio.csv" #_cropped.csv"
dominio = gpd.read_file(dominioPath)
dominio.crs = 'epsg:'+epsg

predioPath=directory+"/GIS/predio.csv"
predio= gpd.read_file(predioPath) 
predio.crs = 'epsg:'+epsg

edifPath=directory+"/GIS/edificios.csv"
edif= gpd.read_file(edifPath) 
edif.crs = 'epsg:'+epsg  

polluts=("CO", "CO","NO2",   "NO2","SO2","SO2","PM10", "PM10" )#, "PM25",  "PM25" #,"NHEX", "NHEPT", "O3") "IPA", 
periods=( "1", "8" ,"1"  ,"ANNUAL",  "1", "24",  "24","ANNUAL")#,  "24" ,"ANNUAL" #,  "24",     "8",  "8")  "24", 

unit="ug/m3"

#request = cimgt.OSM()
projection = ccrs.UTM(zone=utmzone,southern_hemisphere=True)
request=cimgt.GoogleTiles(style='satellite')

borde=10
xmin,ymin,xmax,ymax=dominio.geometry.total_bounds
extent=[xmin-borde,xmax+borde,ymin-borde,ymax+borde]

cols=len(set(polluts))
rows=2
print("cols:",cols,"rows:",rows)
fig,axes = plt.subplots(rows,cols, figsize=(12, 6), sharex=True, sharey=True, gridspec_kw = {'wspace':0.2, 'hspace':0.007}, subplot_kw={'projection': projection, "aspect": 2})

for i in range(len(polluts)):
    print(i)
    pollutid=polluts[i]
    period=periods[i]
    filepath=directory+'/out/mapas/'+proyecto+'_'+pollutid+'_'+period+'.tif'
    X,Y,C=getData(filepath)

    nlevs=6
    minval=max(1e-20,np.nanmin(C))
    maxval=np.nanmax(C)*1.1
    
    levels = np.linspace(minval,maxval,nlevs)
    #levels = np.logspace(np.log10(minval),np.log10(maxval),nlevs)
    #levels = np.vectorize(lambda x: float(f"{x:.2g}"))(levels)

    col=i//rows; row=i%rows
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

    p=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap="RdBu_r", levels=levels, transform=projection)# , norm = LogNorm(vmin=levels[0], vmax=levels[nlevs-1],clip=True))   
    #p=ax.contourf(X,Y,C[:,:],alpha=0.6,cmap="RdBu_r",levels=levels,norm = LogNorm(),transform=projection)#,cmap='RdBu_r',levels=levels    
    dominio.plot(ax=ax, color='none'   ,edgecolor='red'  ,alpha=0.95)
    predio.plot(ax=ax,  color='#f0f0f0',edgecolor='black',alpha=0.75)

    cbar=fig.colorbar(p, ax=ax, shrink=0.5)
    cbar.ax.set_title("Conc. \n ["+unit+"]",loc="left",  fontsize=10)
    ax.set_title(pollutid+" ("+period+")") #,fontsize=6)

#plt.show()
#plt.tight_layout()
plt.savefig('conc_'+proyecto+'.svg',dpi=1400, bbox_inches = 'tight',pad_inches = 0.1)

