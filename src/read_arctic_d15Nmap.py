### Read in the Arctic data 

# Author: Eliza Harris
# Created on Wed Oct 2

#%% Import packages, set up grid

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import scipy as scipy
import math as math
import os
import rasterio as rio
import rasterio.plot as rioplot
from scipy.interpolate import griddata

parentdir = os.getcwd()
os.chdir(parentdir+"/isotone_arcticBranch") # Change to the isotone repo as working directory
import src.utils
from src.utils import plot_map
import configs as configs

#%% Import the gridded dataset from Dan Kou

arctic_d15N_rast = rio.open("data/largeData/202409_pred_soil15N_rast_permGlobal.tif")
fig = plt.figure(figsize=(8,4))
rioplot.show(arctic_d15N_rast)

# set resolution and grid for this dataset
lat_res = 180/arctic_d15N_rast.shape[0]*10
lon_res = 360/arctic_d15N_rast.shape[1]*10
lat_out = np.arange(60, 90, lat_res) # min max spacing
lon_out = np.arange(-180, 180, lon_res)
(LON,LAT) = np.meshgrid(lon_out,lat_out)

arctic_d15N_list = src.utils.rast_to_list(arctic_d15N_rast,10,dataclip=(0,np.inf),takemean="N") # data as list (can reduce resolution here if desired)
# Regrid the data 
# Currently reduces resolution just to speed things up during testing
arctic_d15N_grid = griddata((arctic_d15N_list[:,0],arctic_d15N_list[:,1]), arctic_d15N_list[:,2], (LON,LAT), method='linear')

# Create the arctic mask
arctic_mask = arctic_d15N_grid/arctic_d15N_grid

# Arbitrary uncertainty for now (set at 5 permil)
arctic_d15N_unc = arctic_d15N_grid/arctic_d15N_grid*5

# plot the data
plot_map(LON,LAT,arctic_d15N_grid,"Arctic soil d15N",filename="figs/input_data/arctic_soil_d15N")

##% Create an arctic mask that can be used with the main isoscape, with exactly the same resolution

# set resolution and grid for main d15N data
latm_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lonm_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LONm,LATm) = np.meshgrid(lonm_out,latm_out)
arctic_d15N_coarsegrid = griddata((arctic_d15N_list[:,0],arctic_d15N_list[:,1]), arctic_d15N_list[:,2], (LONm,LATm), method='linear')

plot_map(LON,LAT,arctic_d15N_coarsegrid,"Arctic soil d15N",filename="figs/input_data/arctic_soil_d15N_coarse")

# Create the arctic mask
arctic_mask_coarse = arctic_d15N_coarsegrid/arctic_d15N_coarsegrid

##% Save as a netcdf for easy use

import netCDF4 as nc4

ncout = nc4.Dataset('data/arctic_d15N_data.nc','w','NETCDF4'); # using netCDF3 for output format 
# Save high res
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
soild15N_arctic = ncout.createVariable('soild15N_arctic','f4',('lat','lon'))
soild15N_arctic.setncattr('units','permil')
soild15N_arctic[:,:] = arctic_d15N_grid
arcticmask = ncout.createVariable('arcticmask','f4',('lat','lon'))
arcticmask.setncattr('units','NA')
arcticmask[:,:] = arctic_mask
soild15N_arctic_Unc = ncout.createVariable('soild15N_arctic_Unc','f4',('lat','lon'))
soild15N_arctic_Unc.setncattr('units','permil')
soild15N_arctic_Unc[:,:] = arctic_d15N_unc 
# Save coarse res
ncout.createDimension('lon_coarse',LONm.shape[1])
ncout.createDimension('lat_coarse',LATm.shape[0])
lonvar_coarse = ncout.createVariable('lon_coarse','f4',('lon_coarse'))
lonvar_coarse[:] = LONm[0,:]
lonvar_coarse.setncattr('units','degrees east')
latvar_coarse = ncout.createVariable('lat_coarse','f4',('lat'))
latvar_coarse[:] = LATm[:,0]
latvar_coarse.setncattr('units','degrees north')
soild15N_arctic_coarse = ncout.createVariable('soild15N_arctic_coarse','f4',('lat_coarse','lon_coarse'))
soild15N_arctic_coarse.setncattr('units','permil')
soild15N_arctic_coarse[:,:] = arctic_d15N_coarsegrid
arcticmask_coarse = ncout.createVariable('arcticmask_coarse','f4',('lat_coarse','lon_coarse'))
arcticmask_coarse.setncattr('units','NA')
arcticmask_coarse[:,:] = arctic_mask_coarse
ncout.close()

# test the file!
f = nc4.Dataset('data/climate_input_data.nc','r')
plot_map(LON,LAT,f.variables["soild15N_arctic"][:,:],"soild15N_arctic")
plot_map(LON,LAT,f.variables["arcticmask"][:,:],"arcticmask")
plot_map(LON,LAT,f.variables["arctic_d15N_unc"][:,:],"arctic_d15N_unc")
plot_map(LONm,LATm,f.variables["soild15N_arctic_coarse"][:,:],"soild15N_arctic_coarse")
plot_map(LONm,LATm,f.variables["arctic_mask_coarse"][:,:],"arctic_mask_coarse")