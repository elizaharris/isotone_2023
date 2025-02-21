#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:50:35 2020

This file imports and regrids WFPS from Dan (not with other inputs; added later and annual)

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
import os
from scipy.interpolate import griddata

parentdir = os.getcwd()
os.chdir(parentdir+"/isotone_arcticBranch") # Change to the isotone repo as working directory
from src.utils import plot_map
from src.utils import rast_to_list
import configs as configs

# set resolution and grid for model
resolution = configs.res["ystart"]
lat_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lon_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

#%% Global AI for ocean mask... (not needed, no ocean mask required, copy again from NH3 regrid if needed)

#%% Get the WFPS and check
# For some reason there are lots of zeros/nans in parts, that shouldn't be there... need to ask Dan
wfps_rast = rio.open("../isotone-rawdata/wfps_0-7cm/wfps0007_globalmap_2010.tif")
plt.imshow(wfps_rast.read(1), cmap='pink')
plt.show()

wfps_rast = rio.open("../isotone-rawdata/wfps_7-28cm/wfps0728_globalmap_2010.tif")
plt.imshow(wfps_rast.read(1), cmap='pink')
plt.show()

wfps_list = rast_to_list(wfps_rast,2,dataclip=(0.05,np.inf),takemean="Y") # data as list
# Because of the lines, using less than 5-10 here gives problems
# regrid the data from 1D back to 2D and to a rougher grid (too many values for now otherwise)
wfps_grid = griddata((wfps_list[:,0],wfps_list[:,1]), wfps_list[:,2], (LON,LAT), method='linear')
# wfps_grid = (wfps_grid).clip(0,100)*ocean_mask # not needed
plot_map(LON,LAT,wfps_grid,title="wfps",filename="figs/input_data/wfps_0-7cm_2020",show=1)

#%% Get file names

from os import listdir
from os.path import isfile, join
wfps_filenames_007 = [f for f in listdir("../isotone-rawdata/wfps_0-7cm/") if isfile(join("../isotone-rawdata/wfps_0-7cm/", f))]
wfps_filenames_728 = [f for f in listdir("../isotone-rawdata/wfps_7-28cm/") if isfile(join("../isotone-rawdata/wfps_7-28cm/", f))]
nyears = len(wfps_filenames_007)

# Make sure ordering correct
years_007 = np.array([ int(s.split("_")[2].split(".")[0]) for s in wfps_filenames_007 ])
order = np.argsort(years_007)
years_007 = years_007[order]
wfps_filenames_007 = np.array(wfps_filenames_007)[order]
years_728 = np.array([ int(s.split("_")[2].split(".")[0]) for s in wfps_filenames_728 ])
order = np.argsort(years_728)
years_728 = years_728[order]
wfps_filenames_728 = np.array(wfps_filenames_728)[order]

# Check years
if len(wfps_filenames_728)!=len(wfps_filenames_007):
    print("Different numbers of files for soil depths!")
if np.nanmin(years_007)!=np.nanmin(years_728):
    print("Different start years for different depths")
years = years_007.copy()

#%% Import and regrid each file
wfps_007 = np.zeros((wfps_grid.shape[0],wfps_grid.shape[1],nyears))*np.nan
for n,s in enumerate(wfps_filenames_007):
    print(s)
    wfps_rast = rio.open("../isotone-rawdata/wfps_0-7cm/"+s)
    wfps_list = rast_to_list(wfps_rast,2,dataclip=(0,np.inf),takemean="Y") # data as list 
    wfps_grid = griddata((wfps_list[:,0],wfps_list[:,1]), wfps_list[:,2], (LON,LAT), method='linear')
    wfps_007[:,:,n] = wfps_grid
    if years[n]==2020:
        plot_map(LON,LAT,wfps_007[:,:,n],title="wfps",filename="figs/input_data/wfps_0-7cm_2020")
wfps_007_mean = np.nanmean(wfps_007,axis=2)
plot_map(LON,LAT,wfps_007_mean,title="wfps",filename="figs/input_data/wfps_0-7cm_mean",show=1)

wfps_728 = np.zeros((wfps_grid.shape[0],wfps_grid.shape[1],nyears))*np.nan
for n,s in enumerate(wfps_filenames_728):
    print(s)
    wfps_rast = rio.open("../isotone-rawdata/wfps_7-28cm/"+s)
    wfps_list = rast_to_list(wfps_rast,2,dataclip=(0,np.inf),takemean="Y") # data as list 
    wfps_grid = griddata((wfps_list[:,0],wfps_list[:,1]), wfps_list[:,2], (LON,LAT), method='linear')
    wfps_728[:,:,n] = wfps_grid
    if years[n]==2020:
        plot_map(LON,LAT,wfps_728[:,:,n],title="wfps",filename="figs/input_data/wfps_7-28cm_2020")
wfps_728_mean = np.nanmean(wfps_728,axis=2)
plot_map(LON,LAT,wfps_728_mean,title="wfps",filename="figs/input_data/wfps_7-28cm_mean")

plot_map(LON,LAT,wfps_007_mean-wfps_728_mean,title="wfps, surface - depth (difference)",filename="figs/input_data/wfps_surface-depth")
wfps_mean = (wfps_007_mean+3*wfps_728_mean)/4 # mean, weighted by depth span
plot_map(LON,LAT,wfps_mean,title="wfps",filename="figs/input_data/wfps_new_mean",show=1)

#%% Save each year of the data as a netcdf file

import netCDF4 as nc4

os.remove('data/largeData/wfps_new_data.nc','w','NETCDF4')
ncout = nc4.Dataset('data/largeData/wfps_new_data.nc','w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
ncout.createDimension('year',nyears)
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
yearvar = ncout.createVariable('year','f4',('year'))
yearvar[:] = years
yearvar.setncattr('units','year')
wfps_007var = ncout.createVariable('wfps_007','f4',('lat','lon','year'))
wfps_007var.setncattr('units','percentage (0-7 cm)')
wfps_007var[:,:,:] = wfps_007[:,:,:]
wfps_728var = ncout.createVariable('wfps_728','f4',('lat','lon','year'))
wfps_728var.setncattr('units','percentage (7-28 cm)')
wfps_728var[:,:,:] = wfps_728[:,:,:]
wfps_mean007var = ncout.createVariable('wfps_mean007','f4',('lat','lon'))
wfps_mean007var.setncattr('units','percentage (0-7 cm), mean')
wfps_mean007var[:,:] = wfps_007_mean[:,:]
wfps_mean728var = ncout.createVariable('wfps_mean728','f4',('lat','lon'))
wfps_mean728var.setncattr('units','percentage (7-28 cm), mean')
wfps_mean728var[:,:] = wfps_728_mean[:,:]
wfps_meanvar = ncout.createVariable('wfps_mean','f4',('lat','lon'))
wfps_meanvar.setncattr('units','percentage (0-28 cm), mean')
wfps_meanvar[:,:] = wfps_mean[:,:]
ncout.close()

# test the file!
f = nc4.Dataset('data/largeData/wfps_new_data.nc','r')
plot_map(LON,LAT,f.variables["wfps_mean"][:,:],"wfps mean",show=1)
plot_map(LON,LAT,f.variables["wfps_007"][:,:,20],"wfps",show=1)
