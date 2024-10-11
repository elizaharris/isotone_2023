#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:50:35 2020

This file imports and regrids NH3 from Bai (not with other inputs as added later)

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

#%% Global AI for ocean mask...

AI_rast = rio.open("../isotone-rawdata/PET_AI/ai_et0.tif")
AI_list = rast_to_list(AI_rast,50,dataclip=(0,np.inf),takemean="N") # data as list
# regrid the data from 1D back to 2D and to a rougher grid (too many values for now otherwise)
AI_grid = griddata((AI_list[:,0],AI_list[:,1]), AI_list[:,2], (LON,LAT), method='linear')
ocean_mask = AI_grid/AI_grid

#%% Get the fNH3 from Bai et al.
    
fNH3 = np.genfromtxt(fname = "data/Bai2012_Files/Other/fNH3.txt",delimiter=",",skip_header=1) 
# Columns: X,X,lon,lat,...,fNH3
fNH3[fNH3<-900] = np.nan
fNH3_grid = griddata((fNH3[:,2],fNH3[:,3]),fNH3[:,-1], (LON,LAT), method='linear') # lin gives best results!
fNH3_grid = (fNH3_grid).clip(0,1)*ocean_mask
plot_map(LON,LAT,fNH3_grid,title="f of N emitted as NH3",filename="figs/input_data/bai_input_NH3")

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/fNH3_data.nc','w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
fNH3 = ncout.createVariable('fNH3','f4',('lat','lon'))
fNH3.setncattr('units','fraction')
fNH3[:,:] = fNH3_grid[:,:]
ncout.close()

# test the file!
f = nc4.Dataset('data/fNH3_data.nc','r')
plot_map(LON,LAT,f.variables["fNH3"][:,:],"fNH3")
