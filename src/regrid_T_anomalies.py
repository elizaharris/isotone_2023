#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:50:35 2020

This file imports and regrids temperature anomaly data and saves as netcdf

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
from src.utils import plot_map, reducesize, nc_to_list
import configs as configs

# set resolution and grid for model
resolution = configs.res["ystart"]
lat_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lon_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

# Years that Hadley data is available (doesn't need to be same as model years)
years = np.arange(1850,2019)
Tanom_griddata = np.zeros((len(years),LAT.shape[0],LON.shape[1]))
                          
# get Hadley centre data
landT_anom = xr.open_dataset("../isotone-rawdata/HadClim/CRUTEM.4.6.0.0.anomalies.nc")
lon = landT_anom.longitude.values
for n in range(0,years.shape[0]) :
    print(years[n])    
    data = landT_anom.temperature_anomaly.sel(time=str(years[n]))
    data_mean = np.nanmean(data,axis=0)
    for i in range(0,data_mean.shape[0]) :
        inan = np.where(np.isnan(data_mean[i,:]))
        inotnan = np.where(~np.isnan(data_mean[i,:]))
        if (len(data_mean[i,:])-len(inan[0])>0) :
            x = np.interp(lon[inan],lon[inotnan],data_mean[i,inotnan].flatten())
            data_mean[i,inan] = x
    data_list = nc_to_list(landT_anom.longitude, landT_anom.latitude, data_mean)  
    data_grid = griddata((data_list[:,0], data_list[:,1]), data_list[:,2], (LON,LAT), method='linear') 
    Tanom_griddata[n,:,:]  = data_grid

# plot if needed (last year)
plt.contourf(data_mean)
plot_map(LON,LAT,Tanom_griddata[n,:,:],"T anomaly: "+str(years[n]),filename="figs/input_data/T_anomalies")

# Normalise to the beginning year, so that we have T anomalies relative to the first model year
Tnorm_griddata = np.zeros(Tanom_griddata.shape)
Tstart = Tanom_griddata[0,:,:]
Tstart[np.isnan(Tstart)] = 0
for n in range(0,len(years)) :
    print(years[n])    
    Tnorm_griddata[n,:,:] = Tanom_griddata[n,:,:] - Tstart
    Tnorm_griddata[np.isnan(Tnorm_griddata)] = 0 
plot_map(LON,LAT,Tnorm_griddata[150,:,:],"Temp anomaly: "+str(years[150]),filename="figs/input_data/T_anomalies_norm")
mean_T_rise = np.nanmean(np.nanmean(Tnorm_griddata,axis=1),axis=1)

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/largeData/Tanom_data.nc','w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
ncout.createDimension('time',len(years))
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
timevar = ncout.createVariable('time','f4',('time'))
timevar[:]=years
timevar.setncattr('details','annual total for the given year')
Tanom = ncout.createVariable('Tanom','f4',('time','lat','lon'))
Tanom.setncattr('units','degC')
Tanom.setncattr('source','Hadley centre temperature anomalies')
Tanom[:,:,:] = Tanom_griddata[:,:,:]
Tnorm = ncout.createVariable('Tnorm','f4',('time','lat','lon'))
Tnorm.setncattr('units','degC')
Tnorm.setncattr('source','Hadley centre temperature anomalies, normalised to T in start year')
Tnorm[:,:,:] = Tnorm_griddata[:,:,:]
T_rise = ncout.createVariable('mean_T_rise','f4',('time'))
T_rise.setncattr('units','degC')
T_rise.setncattr('source','Hadley centre temperature anomalies, mean per year across grid cells (not weighted by size)')
T_rise[:] = mean_T_rise[:]
ncout.close()

# test the file!
f = nc4.Dataset('data/largeData/Tanom_data.nc','r')
plot_map(LON,LAT,f.variables["Tanom"][n,:,:],"Temp anomaly: "+str(f.variables["time"][5]))
plot_map(LON,LAT,f.variables["Tnorm"][n,:,:],"Temp anomaly: "+str(f.variables["time"][5]))
fig, ax = plt.subplots(1,1)
plt.plot(years,f.variables["mean_T_rise"][:],"o")
plt.show()
plt.close("all")