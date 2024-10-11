#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Data preprocessing for the IsoTONE model: Generate d15N isoscape

# Author: Eliza Harris
# Created on Thu Jan 23 13:26:42 2020

# Tasks:
# Import MAP, MAT, C and pH data and regrid
# Import d15N data and calculate dependence on MAP, MAP, C and pH using Bayesian MCMC
# Calculate gridded d15N from MAP, MAT, C and pH
# Compare to original d15N values

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

# set resolution and grid for model
resolution = configs.res["ystart"]
lat_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lon_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

#%% Global PET and AI 
# Large dataset, takes a while!
# Use this to generate ocean mask

AI_rast = rio.open("../isotone-rawdata/PET_AI/ai_et0.tif")
fig = plt.figure(figsize=(8,4))
rioplot.show(AI_rast)

AI_list = src.utils.rast_to_list(AI_rast,50,dataclip=(0,np.inf),takemean="N") # data as list
# regrid the data from 1D back to 2D and to a rougher grid (too many values for now otherwise)
AI_grid = griddata((AI_list[:,0],AI_list[:,1]), AI_list[:,2], (LON,LAT), method='linear')

# Create and use the ocean mask
ocean_mask = AI_grid/AI_grid 
AI_grid = AI_grid*ocean_mask

# plot the data
plot_map(LON,LAT,AI_grid,"Aridity index",filename="figs/input_data/aridity_index")

#%% Global soil organic C data
# organic carbon (t ha-1) for the topsoil (0 – 30cm)

globalC_rast = rio.open("../isotone-rawdata/soil_carbon/HWSDa_OC_Dens_Top_5min.rst")
fig = plt.figure(figsize=(8,4))
rioplot.show(globalC_rast)
fig.show() 

globalC_list = src.utils.rast_to_list(globalC_rast,10,takemean="N") # raster data as list
# regrid the data from 1D back to 2D and to a rougher grid (too many values for now otherwise)
globalC_grid = griddata((globalC_list[:,0],globalC_list[:,1]), globalC_list[:,2], (LON,LAT), method='linear')
# convert to mg / g
globalC_grid = globalC_grid*1000*1000/10000 # T / ha to kg / m-2 (same as Zinke dataset)
#globalC_grid = globalC_grid/(10000*100*100*30) # mg / cm-3 for 30 cm depth
#globalC_grid = globalC_grid/(1.33) # mg g-1 considering average BulkDity of 1.33 g cm-3 
globalC_grid = globalC_grid*ocean_mask

# plot the data
plot_map(LON,LAT,globalC_grid,title="organic carbon density (t ha-1) for the topsoil (0 – 30cm)",filename="figs/input_data/soil_carbon")

#%% HADCRU data (from https://crudata.uea.ac.uk/cru/data/temperature/)
# Mean temperature and temperature anomalies

def plot_land(data,filename="figs/testfig"): # how to plot
    fig = plt.figure(figsize=(8,4))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.contourf(landT_anom.longitude, landT_anom.latitude, data)
    fig.tight_layout()
    plt.savefig(filename+".png") 
    plt.savefig(filename+".pdf") 
    fig.show() 

# get T anomalies
landT_anom = xr.open_dataset("../isotone-rawdata/HadClim/CRUTEM.4.6.0.0.anomalies.nc")
data = landT_anom.temperature_anomaly.sel(time="2012-04").values.squeeze() # how to select a month
plot_land(data,filename="figs/input_data/landT_anomaly_HADCRU_2012-04")

# get absolute mean T
landT_abs = xr.open_dataset("../isotone-rawdata/HadClim/absolute.nc")
data = landT_abs.tem.sel(time=12).values.squeeze() # how to select a month
plot_land(data,filename="figs/input_data/landT_abs_HADCRU_12")

#%% Climatology from Hadley centre
# https://crudata.uea.ac.uk/cru/data/hrg/tmc/ 
# https://crudata.uea.ac.uk/cru/data/hrg/

### precipitation
precip = np.loadtxt(fname = "../isotone-rawdata/HadClim/grid_10min_pre.dat")
# lon, lat, 12 monthly values, 12 CV values in %
MAP = np.sum(precip[:,2:14],1) # calc the mean annual precip from the values
# regrid the data from 1D to 2D and to a rougher grid (too many values for now otherwise)
MAP_grid = griddata((precip[:,1],precip[:,0]), MAP, (LON,LAT), method='linear')
MAP_grid = MAP_grid*ocean_mask
plot_map(LON,LAT,MAP_grid,"MAP (mm)",filename="figs/input_data/MAP_HadClim",vminmax=(0,4000))

### temperature
temp = np.loadtxt(fname = "../isotone-rawdata/HadClim/grid_10min_tmp.dat")
# lon, lat, 12 monthly values
MAT = np.mean(temp[:,2:14],1) # calc the mean annual from the values
# regrid the data from 1D to 2D and to a rougher grid (too many values for now otherwise)
MAT_grid = griddata((temp[:,1],temp[:,0]), MAT, (LON,LAT), method='linear')
MAT_grid = MAT_grid*ocean_mask
plot_map(LON,LAT,MAT_grid,"MAT (degC)",filename="figs/input_data/MAT_HadClim")

#%% Global pH data
    
pH = np.genfromtxt(fname = "../isotone-rawdata/soil_pH/pH_10cm.csv",delimiter=",",skip_header=1) 
# Database,Pedon_ID,Latitude,Longitude,Year,pH
pH_grid = griddata((pH[:,3],pH[:,2]),pH[:,5], (LON,LAT), method='linear') # lin gives best results!
pH_grid = pH_grid*ocean_mask

# plot the regridded data
plot_map(LON,LAT,pH_grid,"pH",filename="figs/input_data/soil_pH")

#%% Global soil bulk density data
    
BulkD = np.genfromtxt(fname = "data/Bai2012_Files/other/Bulk_density.txt",delimiter=",",skip_header=1) 
# X,X,lon,lat,BulkD
BulkD_grid = griddata((BulkD[:,2],BulkD[:,3]),BulkD[:,4], (LON,LAT), method='linear') # lin gives best results!
BulkD_grid = (BulkD_grid*ocean_mask).clip(0.7,4)
plot_map(LON,LAT,BulkD_grid,"Bulk Density (g cm-3)",filename="figs/input_data/bulk_density")

#%% Global N data
    
CN = np.genfromtxt(fname = "../isotone-rawdata/soil_N_C_Zinke/zinke_soil_edit.csv",delimiter=",",skip_header=1) 

# C, N, lat, lon, elev
N_grid = griddata((CN[:,3],CN[:,2]),CN[:,1], (LON,LAT), method='linear') # lin gives best results!
tmp = np.where(np.isnan(N_grid))
N_grid[tmp] = np.mean(N_grid[np.where(~np.isnan(N_grid))]) # make NAs to mean to not interrupt ML
tmp = np.where(N_grid==0)
N_grid[tmp] = np.mean(N_grid[np.where(~np.isnan(N_grid))]) # make zeros to mean to not interrupt ML
N_grid = N_grid*ocean_mask
C_grid = griddata((CN[:,3],CN[:,2]),CN[:,0], (LON,LAT), method='linear') # lin gives best results!
C_grid = C_grid*ocean_mask

# plot the regridded data
plot_map(LON,LAT,C_grid,"soil C",filename="figs/input_data/soil_C_zinke")
plot_map(LON,LAT,N_grid,"soil N",filename="figs/input_data/soil_N_zinke")
 
#%% Soil 15N data from Craine et al and other sources
    
soil15N_1 = np.genfromtxt(fname = "../isotone-rawdata/Soil15Ndata/CraineSoil15N_edit.csv",delimiter=",",skip_header=1)  
soil15N_2 = np.genfromtxt(fname = "../isotone-rawdata/Soil15Ndata/15NDataSummarized_edit.csv",delimiter=",",skip_header=0)  
soil15N_3 = np.genfromtxt(fname = "../isotone-rawdata/Soil15Ndata/New_d15N_data.csv",delimiter=",",skip_header=1)     
# columns: 0-4 Latitude	Longitude	MAT (°C)	    MAP (mm)	    AverageDepth (cm)	
# 5-11 d15N (‰)	[N] (mg g-1)	   [C] (mg g-1)	   C:N	   %Sand	   %Silt	   %Clay
# 12-16 pH,Bulk BulkD,Elevation,Reference,Site 

# deal with craine data only first
soil15N_old = np.append(soil15N_1,soil15N_2,axis=0) 
tmp = np.where(~np.isnan(soil15N_old[:,[0,1,5]]).any(axis=1)) # find rows where there are no nans in lon, lat or d15N for interpolation
soil15N_old = soil15N_old[tmp[0],:]
np.savetxt("../isotone-rawdata/Soil15Ndata/CraineDataset_export.txt", soil15N_old)

# regrid the data from 1D to 2D 
soil15N_old_grid = griddata((soil15N_old[:,1],soil15N_old[:,0]), soil15N_old[:,5], (LON,LAT), method='linear') # lin gives best results!
soil15N_old_grid = soil15N_old_grid*ocean_mask

# plot without function to set both points and contour same colour range
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines() 
ax.set_title("soil d15N, old data (permil)")
cont = plt.contourf(LON,LAT,soil15N_old_grid,vmin=-5,vmax=15)
points = plt.scatter(soil15N_old[:,1],soil15N_old[:,0],c=soil15N_old[:,5],s=15,edgecolor="black",linewidth=0.5,vmin=-5,vmax=15)
cbar = plt.colorbar(points,fraction=0.016, pad=0.04) # add colour bar
fig.tight_layout()
plt.savefig("figs/input_data/d15N_craine.png") 
plt.savefig("figs/input_data/d15N_craine.pdf") 
fig.show() 

# Add the rest of the data
soil15N = np.append(soil15N_old,soil15N_3,axis=0) #  all data
tmp = np.where(~np.isnan(soil15N[:,[0,1,5]]).any(axis=1)) # find rows where there are no nans in lon, lat or d15N for interpolation
soil15N = soil15N[tmp[0],:]

# regrid the data from 1D to 2D 
soil15N_grid = griddata((soil15N[:,1],soil15N[:,0]), soil15N[:,5], (LON,LAT), method='linear') # lin gives best results!
soil15N_grid = soil15N_grid*ocean_mask

# plot without function to set both points and contour same colour range
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_title("soil d15N, all data (permil)")
cont = plt.contourf(LON,LAT,soil15N_grid,vmin=-5,vmax=15)
points = plt.scatter(soil15N[:,1],soil15N[:,0],c=soil15N[:,5],s=15,edgecolor="black",linewidth=0.5,vmin=-5,vmax=15)
cbar = plt.colorbar(points,fraction=0.016, pad=0.04) # add colour bar
fig.tight_layout()
plt.savefig("figs/input_data/d15N_all.png") 
plt.savefig("figs/input_data/d15N_all.pdf") 
fig.show() 

soil15N = np.column_stack([soil15N, np.log(soil15N[:,3]), np.log(soil15N[:,7])])
soil15N_names = ["lat","lon","MAT","MAP","depth","d15N","N","C","CtoN","sand","silt","clay","pH","BulkD","elevation","ref","site","logMAP","logC"]

#%% Get WFPS data also (need this for the N2O:N2:NO partitioning in the model)
    
WFPS = np.genfromtxt(fname = "data/Bai2012_Files/Other/WFPS.txt",delimiter=",",skip_header=1) 
# X,X,lon,lat,...,WFPSmean
WFPS_grid = griddata((WFPS[:,2],WFPS[:,3]),WFPS[:,-1], (LON,LAT), method='linear') # lin gives best results!
WFPS_grid = (WFPS_grid*ocean_mask).clip(0,1)
plot_map(LON,LAT,WFPS_grid,"WFPS",filename="figs/input_data/WFPS")

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/largeData/climate_input_data_arctic.nc','w','NETCDF4'); # using netCDF3 for output format 
ncout.createDimension('lon',LON.shape[1])
ncout.createDimension('lat',LAT.shape[0])
lonvar = ncout.createVariable('lon','f4',('lon'))
lonvar[:] = LON[0,:]
lonvar.setncattr('units','degrees east')
latvar = ncout.createVariable('lat','f4',('lat'))
latvar[:] = LAT[:,0]
latvar.setncattr('units','degrees north')
soilC = ncout.createVariable('soilC','f4',('lat','lon'))
soilC.setncattr('units','kg/m2')
soilC[:,:] = globalC_grid[:,:]
MAPrecip = ncout.createVariable('MAPrecip','f4',('lat','lon'))
MAPrecip.setncattr('units','mm')
MAPrecip[:,:] = MAP_grid
MATemp = ncout.createVariable('MATemp','f4',('lat','lon'))
MATemp.setncattr('units','mm')
MATemp[:,:] = MAT_grid
AridIndex = ncout.createVariable('AridIndex','f4',('lat','lon'))
AridIndex.setncattr('AridityIndex','unitless; see documentation')
AridIndex[:,:] = AI_grid
WFPS = ncout.createVariable('WFPS','f4',('lat','lon'))
WFPS.setncattr('WFPS','fraction')
WFPS[:,:] = WFPS_grid
BulkD = ncout.createVariable('BulkD','f4',('lat','lon'))
BulkD.setncattr('BulkD','g cm-3')
BulkD[:,:] = BulkD_grid
soilpH = ncout.createVariable('soilpH','f4',('lat','lon'))
soilpH[:,:] = pH_grid
soilN = ncout.createVariable('soilN','f4',('lat','lon'))
soilN.setncattr('units','N: g m-2')
soilN[:,:] = N_grid
soild15N_craine = ncout.createVariable('soild15N_craine','f4',('lat','lon'))
soild15N_craine.setncattr('units','permil')
soild15N_craine[:,:] = soil15N_grid 
ncout.close()

# test the file!
f = nc4.Dataset('data/largeData/climate_input_data_arctic.nc','r')
plot_map(LON,LAT,f.variables["soilC"][:,:],"globalC_grid")
plot_map(LON,LAT,f.variables["MAPrecip"][:,:],"MAP")
plot_map(LON,LAT,f.variables["MATemp"][:,:],"MAT")
plot_map(LON,LAT,f.variables["soilpH"][:,:],"pH")
plot_map(LON,LAT,f.variables["soilN"][:,:],"N")
plot_map(LON,LAT,f.variables["soild15N_craine"][:,:],"d15N")
plot_map(LON,LAT,f.variables["WFPS"][:,:],"WFPS")
