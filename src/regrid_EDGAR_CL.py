#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:50:35 2020

This file imports and regrids EDGAR N2O emissions

@author: elizaharris

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import xarray as xr
import cartopy.crs as ccrs
from datetime import datetime as dt
import rasterio as rio
from scipy.interpolate import griddata
import os

parentdir = os.getcwd()
os.chdir(parentdir+"/isotone_arcticBranch") # Change to the isotone repo as working directory
#from src.utils import plot_map, reducesize, nc_to_list
#import configs as configs

# define function for plotting    
def plot_map(longi,lati,gridval,title="title",vminmax=(np.nan,np.nan),cmap="viridis",filename="figs/testfig") :
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(title)
    if np.isnan(vminmax[0]) :
        pr = plt.contourf(longi,lati,gridval,cmap=cmap)
    else : 
        #gridval = gridval.clip(vminmax[0],vminmax[1])
        pr = plt.contourf(longi,lati,gridval,vmin=vminmax[0],vmax=vminmax[1],cmap=cmap)
    cbar = plt.colorbar(pr,fraction=0.016, pad=0.04)
    fig.tight_layout()
    plt.savefig(filename+".png") 
    plt.savefig(filename+".pdf") 
    fig.show() 


# function to subsample means of various numbers of gridcells (otherwise regridding too slow)
def reducesize(lon,lat,arrdata,Xred) :
    lon = lon.values
    lat = lat.values
    arrdata = arrdata.values.squeeze()
    output = np.zeros((round(arrdata.shape[0]/Xred)*round(arrdata.shape[1]/Xred),3))
    z = 0
    for n in range(0,round(arrdata.shape[0]/Xred)) :
        for i in range(0,round(arrdata.shape[1]/Xred)) :
            nmid = int(n*Xred+Xred/2) # mid point for lat
            imid = int(i*Xred+Xred/2) # mid point for lon
            tmp = (lon[imid],lat[nmid])
            output[z,0] = tmp[0]
            output[z,1] = tmp[1]
            # take mean of all values within xred window
            nlen = len(range(int(nmid-Xred/2),int(nmid+Xred/2)))
            ilen = len(range(int(imid-Xred/2),int(imid+Xred/2)))
            nvals = np.repeat(range(int(nmid-Xred/2),int(nmid+Xred/2)),ilen) # select for lat: repeats each element X times
            ivals = np.tile(range(int(imid-Xred/2),int(imid+Xred/2)),nlen) # select for lon: repeats whole range X times
            output[z,2] = np.nanmean(arrdata[nvals,ivals])
            z = z+1
    return output

def nc_to_list(loni,lati,data) :
    output = np.zeros((loni.shape[0]*lati.shape[0],3))
    output[:,0] = np.tile(loni,lati.shape[0])
    output[:,1] = np.repeat(lati,loni.shape[0])
    output[:,2] = data.flatten()
    return output

configs = {"res": {
    "ystart": 1800, # Years to model
    "yend": 2023,
    "lonstart": -180, # Longitude range and resolution
    "lonend": 180,
    "lonres": 0.05,
    "latstart": 30, # Latitude range and resolution
    "latend": 90,
    "latres": 0.05,
    }
    }

# set resolution and grid for model
resolution = configs["res"]["ystart"]
lat_out = np.arange(configs["res"]["latstart"], configs["res"]["latend"], configs["res"]["latres"]) # min max spacing
lon_out = np.arange(configs["res"]["lonstart"], configs["res"]["lonend"], configs["res"]["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

#%% import and regrid the edgar data
    
years = np.arange(1970,2016)
cats = ("AGS","TRO","WWT","IDE","ENE","CHE")
# AGS = Agricultural soils: 4C+4D1+4D2+4D4 / 3C2+3C3+3C4+3C7
# TRO = Road transportation: 1A3b / 1A3b
# WWT = Waste water handling: 6B / 4D
# IDE = Indirect emissions from NOx and NH3: 7B+7C / 5A
# ENE = Energy, 1A1
# CHE = Chemical processes 2B
#cats = ['AGS']

for i in cats:
    print(i,"of",cats)
    EDGAR_griddata = np.zeros((len(years),LAT.shape[0],LON.shape[1]))
    for n in range(0,years.shape[0]) :
        print(years[n])
        if (i == "ENE") | (i=="CHE"):
            filename = "../isotone-rawdata/EDGAR/"+i+"_nc/v6.0_N2O_"+str(years[n])+"_"+i+".0.1x0.1.nc" 
        else: 
            filename = "../isotone-rawdata/EDGAR/"+i+"_nc/v50_N2O_"+str(years[n])+"_"+i+".0.1x0.1.nc"    
        EDGAR = xr.open_dataset(filename)
        EDGAR_list = reducesize(lon = EDGAR.lon, lat = EDGAR.lat, arrdata = EDGAR.emi_n2o, Xred = 5) # reduce size and convert to list
        EDGAR_list[EDGAR_list[:,0]>179.5,0] = EDGAR_list[EDGAR_list[:,0]>179.5,0]-360
        # regrid the data 
        EDGAR_grid = griddata((EDGAR_list[:,0], EDGAR_list[:,1]), EDGAR_list[:,2], (LON,LAT), method='linear')
        # change the units from kg N2O /m2 s1 into g N2O-N m-2 y-1
        EDGAR_grid2 = EDGAR_grid*1000*60*60*24*365
        EDGAR_griddata[n,:,:] = EDGAR_grid2
    exec(i + "_griddata = EDGAR_griddata")

# plot if needed (last year)
#plot_map(LON,LAT,AGS_griddata[n,:,:],"AGS emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_AGS")
#plot_map(LON,LAT,TRO_griddata[n,:,:],"TRO emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_TRO")
#plot_map(LON,LAT,IDE_griddata[n,:,:],"IDE emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_IDE")
#plot_map(LON,LAT,WWT_griddata[n,:,:],"WWT emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_WWT")
#plot_map(LON,LAT,ENE_griddata[n,:,:],"ENE emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_ENE")
#plot_map(LON,LAT,CHE_griddata[n,:,:],"CHE emissions (g N2O-N/m2/year): "+str(years[n]),filename="figs/input_data/EDGAR_CHE")

#%% Extrapolate each grid cell to 1850-2050

from sklearn.linear_model import LinearRegression

# Use linear regression on each grid cell for each category to extrapolate to periods outside the EDGAR range
startyear = 1850 # EDGAR emissions will be 0 in this year (this can be after the model period starts)
nyears = np.arange(startyear,2050)
for c in cats:
    print(c,"of",cats)
    exec("EDGAR_griddata = "+ c + "_griddata")
    EDGAR_extrap = np.zeros((len(nyears),EDGAR_griddata.shape[1],EDGAR_griddata.shape[2]))+np.nan
    for n in np.arange(0,EDGAR_griddata.shape[1]) :
        for i in np.arange(0,EDGAR_griddata.shape[2]) :
            if sum(~np.isnan(EDGAR_griddata[:,n,i]))>2 :
                # Extrapolate early data using the first 15 years
                tmp = years <= min(years)+15
                m = LinearRegression(fit_intercept=False).fit(years[tmp].reshape(-1, 1)-startyear,EDGAR_griddata[tmp,n,i])
                fit = (nyears-startyear)*m.coef_ + m.intercept_
                fit[fit<0] = 0
                EDGAR_extrap[:,n,i] = fit
                # Extrapolate current data using the last 15 years
                tmp = years >= max(years)-15
                m = LinearRegression().fit(years[tmp].reshape(-1, 1)-startyear,EDGAR_griddata[tmp,n,i])
                fit = (nyears-startyear)*m.coef_ + m.intercept_
                fit[fit<0] = 0
                tmp = nyears >= max(years)-15
                EDGAR_extrap[tmp,n,i] = fit[tmp]
                # Use the real data where available
                a = np.where(nyears == years[0])[0][0]
                EDGAR_extrap[np.arange(a,(a+len(years))),n,i] = EDGAR_griddata[:,n,i] 
    exec(c + "_extrap = EDGAR_extrap")

    # Also plot
    print(c,"of",cats)
    exec("EDGAR_griddata = "+ c + "_griddata")
    exec("EDGAR_extrap = "+ c + "_extrap")
    fig, ax = plt.subplots(1,1)
    ax.plot(nyears,np.nansum(np.nansum(EDGAR_extrap[:,:,:],axis=1),axis=1))
    ax.plot(years,np.nansum(np.nansum(EDGAR_griddata[:,:,:],axis=1),axis=1))
    ax.set_ylabel(c)
    fig.tight_layout()
    filename = "figs/input_data/EDGAR_"+c+"_extrap"
    plt.savefig(filename+".png") 
    plt.savefig(filename+".pdf") 
    #fig.show() 

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/EDGAR_data_extrap.nc','w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
ncout.createDimension('time',len(nyears))
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
timevar = ncout.createVariable('time','f4',('time'))
timevar[:]=nyears
timevar.setncattr('details','annual total for the given year')
AGS = ncout.createVariable('AGS','f4',('time','lat','lon'))
AGS.setncattr('units','g N2O-N/m2/year')
AGS.setncattr('source','EDGAR category AGS = Agricultural soils: 4C+4D1+4D2+4D4 / 3C2+3C3+3C4+3C7')
AGS[:,:,:] = AGS_extrap[:,:,:]
TRO = ncout.createVariable('TRO','f4',('time','lat','lon'))
TRO.setncattr('units','g N2O-N/m2/year')
TRO.setncattr('source','EDGAR category TRO = Road transportation: 1A3b / 1A3b')
TRO[:,:,:] = TRO_extrap[:,:,:]
WWT = ncout.createVariable('WWT','f4',('time','lat','lon'))
WWT.setncattr('units','g N2O-N/m2/year')
WWT.setncattr('source','EDGAR category WWT = Waste water handling: 6B / 4D')
WWT[:,:,:] = WWT_extrap[:,:,:]
IDE = ncout.createVariable('IDE','f4',('time','lat','lon'))
IDE.setncattr('units','g N2O-N/m2/year')
IDE.setncattr('source','EDGAR category IDE = Indirect emissions from NOx and NH3: 7B+7C / 5A')
IDE[:,:,:] = IDE_extrap[:,:,:]
CHE = ncout.createVariable('CHE','f4',('time','lat','lon'))
CHE.setncattr('units','g N2O-N/m2/year')
CHE.setncattr('source','EDGAR category CHE = Chemical industry: 2B')
CHE[:,:,:] = CHE_extrap[:,:,:]
ENE = ncout.createVariable('ENE','f4',('time','lat','lon'))
ENE.setncattr('units','g N2O-N/m2/year')
ENE.setncattr('source','EDGAR category ENE = Power industry: 1A1')
ENE[:,:,:] = ENE_extrap[:,:,:]
ncout.close()

# test the file!
f = nc4.Dataset('data/EDGAR_data_extrap.nc','r')
n = 100
plot_map(LON,LAT,f.variables["AGS"][n,:,:],"AGS (g N2O-N/m2/year): "+str(f.variables["time"][n]))
#plot_map(LON,LAT,f.variables["TRO"][n,:,:],"TRO (g N2O-N/m2/year): "+str(f.variables["time"][n]))
#plot_map(LON,LAT,f.variables["IDE"][n,:,:],"IDE (g N2O-N/m2/year): "+str(f.variables["time"][n]))
#plot_map(LON,LAT,f.variables["WWT"][n,:,:],"WWT (g N2O-N/m2/year): "+str(f.variables["time"][n]))
#plot_map(LON,LAT,f.variables["CHE"][n,:,:],"CHE (g N2O-N/m2/year): "+str(f.variables["time"][n]))
#plot_map(LON,LAT,f.variables["ENE"][n,:,:],"ENE (g N2O-N/m2/year): "+str(f.variables["time"][n]))
