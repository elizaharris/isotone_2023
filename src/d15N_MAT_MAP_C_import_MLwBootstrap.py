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

# plot the raw data
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines() 
ax.set_title("pH, points")
points = plt.scatter(pH[:,3],pH[:,2],c=pH[:,5],s=15,vmin=3.2,vmax=9.6)
cbar = plt.colorbar(points,fraction=0.016, pad=0.04) # add colour bar
fig.tight_layout()
plt.savefig("figs/input_data/soil_pH_Points.png") 
plt.savefig("figs/input_data/soil_pH_Points.pdf") 
fig.show() 

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

# plot the raw data
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines() 
ax.set_title("soil N, points")
points = plt.scatter(CN[:,3],CN[:,2],c=CN[:,1],s=15,vmin=0,vmax=8000)
cbar = plt.colorbar(points,fraction=0.016, pad=0.04) # add colour bar
fig.tight_layout()
plt.savefig("figs/input_data/soil_N_Points.png") 
plt.savefig("figs/input_data/soil_N_Points.pdf") 
fig.show() 
 
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

#%% Compare dataset values of MAT, MAP, C_rast, C_Zinke, N, CN, pH to gridded dataset values

# fix C and N grid units (based on last run of this code - Zinke dataset is currently kg/m2 and g/m2, not mg/g like the Craine data!)
N_grid = N_grid/100
C_grid = C_grid*10
globalC_grid = globalC_grid*500/15000

soil15N_ancdata = soil15N+np.nan
for n in range(0,soil15N.shape[0]) :
    # find matching index
    tmp1 = abs(lat_out - soil15N[n,0])
    tmp1 = np.where(tmp1 == min(tmp1))
    if len(tmp1[0]) > 1: tmp1 = tmp1[0][0]
    tmp2 = abs(lon_out - soil15N[n,1])
    tmp2 = np.where(tmp2 == min(tmp2))
    if len(tmp2[0]) > 1: tmp2 = tmp2[0][0]
    # match MAT
    soil15N_ancdata[n,2] = MAT_grid[tmp1,tmp2]
    # match MAP
    soil15N_ancdata[n,3] = MAP_grid[tmp1,tmp2]
    # match C, N, CN to Zinke dataset
    soil15N_ancdata[n,6] = N_grid[tmp1,tmp2]
    soil15N_ancdata[n,7] = C_grid[tmp1,tmp2]
    soil15N_ancdata[n,8] = soil15N_ancdata[n,7]/soil15N_ancdata[n,6]
    # match C to raster dataset
    soil15N_ancdata[n,7] = globalC_grid[tmp1,tmp2] # overwrites C from Zinke!
    # match pH
    soil15N_ancdata[n,12] = pH_grid[tmp1,tmp2]
    # match AI (put in column that has depth in the d15N dataset)
    soil15N_ancdata[n,4] = AI_grid[tmp1,tmp2]
    # match BulkD
    soil15N_ancdata[n,13] = BulkD_grid[tmp1,tmp2]

#%% Look at matching between dataset and gridded ancdata

import scipy.stats as stats
cols = (2,3,6,7,8,12,13) # Choose cols to compare between gridded and point data
fig, ax = plt.subplots(3,3,figsize=(10,8))
x,y = (0,0)
corrs = np.zeros((len(cols),6)) # slope, int, R2, p, MAD, MRD
for n in range(0,len(cols)) :
    im = ax[x,y].scatter(soil15N[:,cols[n]],soil15N_ancdata[:,cols[n]],c=soil15N[:,0]) # Coloured by lat
    ax[x,y].set_ylabel(soil15N_names[cols[n]]+", grid")
    ax[x,y].set_xlabel(soil15N_names[cols[n]]+", point")
    # set lims and plot 1:1 line
    tmp = np.vstack((soil15N[:,cols[n]],soil15N_ancdata[:,cols[n]]))
    tmp[np.isinf(tmp)] = np.nan
    low = np.nanmin(tmp)
    high = np.nanmax(tmp)
    ax[x,y].plot((low-0.1*(high-low),high+0.1*(high-low)),(low-0.1*(high-low),high+0.1*(high-low)),"k:") # 1:1 line
    ax[x,y].set_xlim((low-0.1*(high-low),high+0.1*(high-low)))
    ax[x,y].set_ylim((low-0.1*(high-low),high+0.1*(high-low)))
    if soil15N_names[cols[n]] == "CtoN": 
        ax[x,y].set_xlim((-10,400))
        ax[x,y].set_ylim((-10,400))
    # look at correlation
    mask = ~np.isnan(soil15N[:,cols[n]]) & ~np.isnan(soil15N_ancdata[:,cols[n]])
    tmp = stats.linregress(soil15N[mask,cols[n]],soil15N_ancdata[mask,cols[n]])
    corrs[n,0:4] = tmp[0:4]
    if tmp[3] < 0.01: ax[x,y].plot(soil15N[:,cols[n]],soil15N[:,cols[n]]*tmp[0]+tmp[1],"r-")
    if tmp[3] < 0.05: ax[x,y].plot(soil15N[:,cols[n]],soil15N[:,cols[n]]*tmp[0]+tmp[1],"r--")
    corrs[n,4] = np.nanmean(abs(soil15N[mask,cols[n]]-soil15N_ancdata[mask,cols[n]]))
    corrs[n,5] = corrs[n,4]/np.nanmean(soil15N[mask,cols[n]])
    if n==0: 
        ax[x,y].legend((["data","1:1","fit, p<0.01","fit, p<0.05"]),fontsize=5)
        fig.colorbar(im)
    x = x+1
    if x==3: x,y=(0,y+1)
fig.tight_layout()
plt.savefig("figs/input_data/comparison_point-grid.png") 
plt.savefig("figs/input_data/comparison_point-grid.pdf") 
fig.show() 

# Summarise comparison 
corrs = pd.DataFrame(corrs)
corrs.columns = ["slope", "int", "R2", "p", "MAD", "MRD"]
corrs.index = [ soil15N_names[c] for c in cols ]
print(corrs)
    
# make a final dataset with best data
soil15N_final = soil15N.copy() # Use point measurement data as a starting point
soil15N_final[:,2] = soil15N_ancdata[:,2] # MAT: Good match: Use grid data for consistency
soil15N_final[:,3] = soil15N_ancdata[:,3] # MAP: Good match: Use grid data for consistency
tmp = (np.isnan(soil15N[:,6])) | (soil15N[:,6] > 60) # N: Poor match: Only use grid data for gaps
soil15N_final[tmp,6] = soil15N_ancdata[tmp,6] 
tmp = np.isnan(soil15N[:,7])  # C: Fairly poor match: Only use grid data for gaps
soil15N_final[tmp,7] = soil15N_ancdata[tmp,7]
soil15N_final[:,8] = soil15N_ancdata[:,7]/soil15N_ancdata[:,6]  # C:N: Poor match: Calc from filled C and N
soil15N_final[:,12] = soil15N_ancdata[:,12] # pH: Good match: Use grid data for consistency
soil15N_final[:,13] = soil15N_ancdata[:,13] # BulkD: Good match: Use grid data for consistency
soil15N_final[:,16] = soil15N_ancdata[:,4] # aridity index (no point data)
soil15N_final[:,17] = np.log(soil15N_final[:,6]) # Log transformed values
soil15N_final[:,18] = np.log(soil15N_final[:,7])
soil15N_final = np.hstack( (soil15N_final, np.zeros((soil15N_final.shape[0], 1)) ) )
soil15N_final[:,19] = abs(soil15N_final[:,1].copy()) # Absolute longitude
soil15N_names = np.array(("lat","lon","MAT","MAP","depth","d15N","N","C","CtoN","sand","silt","clay","pH","BulkD","elevation","ref","AI","logN","logC","lon_rescaled"))
    
#%% Set up to model d15N isoscape using NN with Keras (weights, data reorganisation, uncertainty...)

import keras as k
from keras.models import Sequential
from keras.layers import Activation,Dense,BatchNormalization
from keras.optimizers import Adam
from keras.metrics import mean_squared_error
from keras.metrics import mean_absolute_error
from keras.callbacks import EarlyStopping,ModelCheckpoint
from keras import regularizers
from keras.models import load_model
from sklearn.linear_model import LinearRegression
from numpy import mean, random

# Select features for the ML 
xcols = ("lon_rescaled","lat","MAP","MAT","N","C","CtoN","pH","BulkD","AI"); xrefs = np.zeros(len(xcols))

# Create also a weighting term by lat and lon clustering
from sklearn.cluster import KMeans
tmp = np.where((soil15N_names == "lat") | (soil15N_names == "lon"))[0]
data = soil15N_final[:,tmp]
kmeans = KMeans(n_clusters=30, random_state=0).fit(data)
# Calculate weights based on inverse number of entries in each cluster
w = np.zeros(len(kmeans.labels_))
for n in np.unique(kmeans.labels_):
    tmp = np.where(kmeans.labels_ == n)[0]
    w[tmp] = 1/(len(tmp)**0.5)
w = w/np.nanmax(w)
# Plot the clusters
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_title("k means cluster")
plt.scatter(data[:,1],data[:,0],c=kmeans.labels_,cmap='rainbow')
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/k-means_clusters.png") 
plt.savefig("figs/isoscape_submodel/k-means_clusters.pdf") 
fig.show() 
# Plot the weights
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_title("weight from k means cluster members")
plt.scatter(data[:,1],data[:,0],c=w,cmap='cool')
plt.colorbar()
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/k-means_weights.png") 
plt.savefig("figs/isoscape_submodel/k-means_weights.pdf") 
fig.show() 

# Prep data for the isoscape grid prediction; all ancdata on grid shape, but vectorised
gridtovec = np.where(~np.isnan(MAT_grid))
inputmatrix = np.zeros((len(gridtovec[0]),len(xcols)))-1
for n in range(0,len(xcols)):
    # select the appropriate grid
    if xcols[n] == "MAT" : grid = MAT_grid
    elif xcols[n] == "MAP" : grid = MAP_grid
    elif xcols[n] == "AI" : grid = AI_grid
    elif xcols[n] == "logC" : grid = np.log(globalC_grid)
    elif xcols[n] == "logN" : grid = np.log(N_grid) 
    elif xcols[n] == "BulkD" : grid = BulkD_grid
    elif xcols[n] == "N" : grid = N_grid  
    elif xcols[n] == "C" : grid = globalC_grid 
    elif xcols[n] == "CtoN" : 
        grid = globalC_grid/N_grid
        grid[np.where(abs(grid)>1e4)] = np.nan
        grid[np.where(abs(grid)>1e3)] = 1e3 # this is the maximum order of mag that the model was trained on
    elif xcols[n] == "pH" : grid = pH_grid
    elif xcols[n] == "lon" : grid = LON
    elif xcols[n] == "lon_rescaled" : grid = abs(LON)
    elif xcols[n] == "lat" : grid = LAT
    else : print("param not found!!",xcols[n])
    # vectorise
    vec = grid[gridtovec]
    # add to matrix
    inputmatrix[:,n] = vec

# Look at the standard deviation of sites in a grid cell to get a measure of how representative d15N is
# Could also consider at multiple temporal samplings at any site to estimate this
df = pd.DataFrame({"lon":[0],"lat":[0],"mean":[0],"std":[0],"n":[0]})
grid_representativity = df
for i in np.arange(LON.shape[0]):
    for j in np.arange(LON.shape[1]):
        lon = (LON[i,j]-0.25,LON[i,j]+0.25)
        lat = (LAT[i,j]-0.25,LAT[i,j]+0.25)
        tmp1 = (soil15N_final[:,1] >= lon[0]) & (soil15N_final[:,1] <= lon[1])
        tmp2 = (soil15N_final[:,0] >= lat[0]) & (soil15N_final[:,0] <= lat[1])
        tmp = np.where(tmp1 & tmp2)[0]
        if len(tmp)>2:
            df["lon"] = LON[i,j]
            df["lat"] = LAT[i,j]
            df["mean"] = np.nanmean(soil15N_final[tmp,5])
            df["std"] = np.nanstd(soil15N_final[tmp,5])
            df["n"] = len(tmp)
            grid_representativity = grid_representativity._append(df)
d15N_reperr = np.nanmean(grid_representativity["std"]) # 1.295 permil, based on 461 grid cells
print(str(grid_representativity.shape[0])+" cells with >3 measurements")
print("Average of",str(round(np.nanmean(grid_representativity["n"]),0))+" measurements per gridcell")
print("Mean std of",str(round(d15N_reperr,3))+" across gridcell")
print("Mean std of",str(round(np.nanmean(soil15N_final[:,5]),3))+" across all data")

#%% Train model of d15N isoscape using NN with Keras 

# Select cols to use
xrefs = [5] + [ np.where(soil15N_names == c)[0][0] for c in xcols ] # Cols corresponding to d15N then the chosen input vars
soil15N_input = soil15N_final[:,xrefs]

# Remove nans
soil15N_input[np.isinf(soil15N_input)] = np.nan
notnans = np.sum(np.isnan(soil15N_input),axis=1)==0
print(str(soil15N_input.shape[0]),"values;",str(sum(notnans)),"rows without nans")
soil15N_input = soil15N_input[notnans,:]
w_input = w[notnans]

# Reserve 15% of the data for testing
np.random.seed(20)
tsplit = np.random.uniform(0,1,len(soil15N_input))<0.15 # True = test data, False = training data
soil15N_test = soil15N_input[tsplit,:]
soil15N_trval = soil15N_input[~tsplit,:]
w_test = w_input[tsplit]
w_trval = w_input[~tsplit]

# Rescale input data by mean and stdev (not d15N)
for n in np.arange(1,soil15N_input.shape[1]):
    tmp = soil15N_trval[:,n].copy()
    # tmp = (tmp - min(soil15N_input[:,n]))/(max(soil15N_input[:,n]) - min(soil15N_input[:,n])) old scaling
    tmp = (tmp-np.nanmean(soil15N_input[:,n]))/np.nanstd(soil15N_input[:,n])
    soil15N_trval[:,n] = tmp
    tmp = soil15N_test[:,n].copy()
    tmp = (tmp-np.nanmean(soil15N_input[:,n]))/np.nanstd(soil15N_input[:,n])
    soil15N_test[:,n] = tmp
# Check input matrix (gridded inputs) somewhat matches point inputs before scaling
a = np.nanmean(inputmatrix,axis=0)
a_s = np.nanmean(inputmatrix,axis=0)
b = np.nanmean(soil15N_input[:,1:],axis=0)
b_s = np.nanstd(soil15N_input[:,1:],axis=0)
diffs_std = (a-b)/b_s
print(diffs_std)
print(a_s/b_s)
# Scale input matrix (gridded inputs)
inputmatrix_scale = inputmatrix.copy()
for n in np.arange(1,soil15N_input.shape[1]):
    tmp = inputmatrix[:,n-1].copy()
    tmp = (tmp-np.nanmean(soil15N_input[:,n]))/np.nanstd(soil15N_input[:,n])
    tmp[tmp < 0.8*min(soil15N_trval[:,n])] = 0.8*min(soil15N_trval[:,n]) # Cut to close to the same min and max and the point inputs to avoid big extraps...
    tmp[tmp > 1.2*max(soil15N_trval[:,n])] = 1.2*max(soil15N_trval[:,n])
    inputmatrix_scale[:,n-1] = tmp

# Bootstrapping to check uncertainty
n_bootstraps = 250
bootstrap_d15N_grid_calc = np.zeros((len(gridtovec[0]),n_bootstraps))+np.nan
bootstrap_rmse = np.zeros((n_bootstraps))*np.nan
bootstrap_R2 = bootstrap_rmse.copy()
bootstrap_mean = bootstrap_rmse.copy()
bootstrap_std = bootstrap_rmse
for bootstrap in np.arange(0,n_bootstraps): # go through bootstrapping runs
    print("bootstrap number",str(bootstrap))
    # add uncertainty to the d15N
    soil15N_thisit = soil15N_trval[:,0].copy()
    unc = (d15N_reperr - min(soil15N_input[:,0]))/(max(soil15N_input[:,0]) - min(soil15N_input[:,0]))
    soil15N_thisit = soil15N_thisit + np.random.normal(0,1,len(soil15N_thisit))*(unc)
    # random numbers length of dataset to split into training and validation data 
    split = np.random.uniform(0,1,len(soil15N_thisit))>0.75 # True = val data, False = training data
    # get y data
    y_train = soil15N_thisit[~split].copy().reshape(-1,1)
    y_val = soil15N_thisit[split].copy().reshape(-1,1)
    # get weights
    w_train = w_trval[~split]
    w_test = w_trval[split]
    # get x data, scale by columns
    x_train = soil15N_trval[~split,1:].copy()#.reshape(-1,1)
    x_val = soil15N_trval[split,1:].copy()#.reshape(-1,1)
    # compile model
    model = Sequential([
        Dense(units=16, input_shape=(len(xcols),), activation='relu'), # first hidden layer, 1 input variable specified
        Dense(units=16, activation='tanh',kernel_regularizer=regularizers.l2(0.01)), # second hidden layer
        Dense(units=16, activation='linear',kernel_regularizer=regularizers.l2(0.01)), # second hidden layer
        BatchNormalization(axis=1),
        Dense(units=1, activation='relu') # output layer 
        ])
    model.compile(optimizer=Adam(learning_rate=0.0005), loss='mean_squared_error', metrics = ['accuracy'], weighted_metrics=['accuracy'])
    # markus recommends mean absolute error, least median squares and least trimmed squares as more robust to outliers
    # define early stopping and saving of best model
    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=30)
    from datetime import date; date = date.today().strftime("%Y%m%d")
    mc = ModelCheckpoint("data/d15Nmodels/d15N_ANN_"+date+"_bs"+str(bootstrap)+".h5", monitor='val_loss', mode='min', save_best_only=True, verbose=0)
    # train the model
    bs = 100 # set the batch size (used in predictions also)
    history = model.fit(x=x_train, y=y_train, validation_split=0.2, batch_size=bs, epochs=80, shuffle=True, verbose=0, callbacks=[es,mc], sample_weight=w_train)
    # test the model and save the stats
    smodel = load_model("data/d15Nmodels/d15N_ANN_"+date+"_bs"+str(bootstrap)+".h5")
    predictions = smodel.predict(x=x_val, batch_size=bs, verbose=2)
    bootstrap_rmse[bootstrap] = ( np.nansum(w_test*((predictions-y_val)**2)[:,0])/np.nansum(w_test) )**0.5 # weighted RMSE
    #bootstrap_rmse[bootstrap] = (np.nanmean((predictions-y_test)**2))**0.5 # unweighted RMSE
    print("w-rmse:",str(bootstrap_rmse[bootstrap]))
    lm_model = LinearRegression().fit(predictions,y_val,sample_weight=w_test)
    bootstrap_R2[bootstrap] = lm_model.score(predictions,y_val,sample_weight=w_test)
    print("r2:",str(bootstrap_R2[bootstrap]))
    # check if this is the best bootstrap model so far; rename if it is
    if bootstrap_rmse[bootstrap] == np.nanmin(bootstrap_rmse):
        if bootstrap_R2[bootstrap] >= np.nanmax(bootstrap_R2)-0.02:
            current_best = bootstrap
            # os.rename("d15Nmodels/d15N_ANN_"+date+"_bs"+str(bootstrap)+".h5","d15Nmodels/d15N_ANN_"+date+"_bsBEST.h5") 
    # run prediction for the grid, also convert back to d15N scale
    predforgrid = smodel.predict(x=inputmatrix_scale, batch_size=bs, verbose=2)
    bootstrap_d15N_grid_calc[:,bootstrap] = predforgrid[:,0]
    
current_best = np.where(bootstrap_R2==np.max(bootstrap_R2))[0][0]
current_best = np.where(bootstrap_rmse==np.min(bootstrap_rmse))[0][0]

# Copy and rename the best bootstrap run
from shutil import copyfile
copyfile("data/d15Nmodels/d15N_ANN_"+date+"_bs"+str(current_best)+".h5","data/d15Nmodels/d15N_ANN_"+date+"_bsBEST.h5")

# show the training progress for the last bootstrap
plt.figure()
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='validate')
plt.legend()
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/trainingprogress.png") 
plt.savefig("figs/isoscape_submodel/trainingprogress.pdf") 
fig.show() 

# Plot the current best BS
tmp = MAT_grid.copy()*np.nan
tmp[gridtovec] = bootstrap_d15N_grid_calc[:,current_best]
plot_map(LON,LAT,tmp,"d15N, current best model",filename="figs/isoscape_submodel/map_bestBS")

# Std and mean of the bootstrap results
meangrid = MAT_grid.copy()*np.nan
meangrid[gridtovec] = np.nanmean(bootstrap_d15N_grid_calc,axis=1)
plot_map(LON,LAT,meangrid,"d15N, mean of bootstraps",filename="figs/isoscape_submodel/map_mean-of-BS")
stdgrid = MAT_grid.copy()*np.nan
stdgrid[gridtovec] = np.nanstd(bootstrap_d15N_grid_calc,axis=1)
plot_map(LON,LAT,stdgrid,"d15N, stdev of bootstraps",filename="figs/isoscape_submodel/map_std-of-BS")

#%% load and use the best model (can just load from here if not running the model new)

# Adapt to use the "current best" instead next time!
smodel = load_model("data/d15Nmodels/d15N_ANN_20231024_bsBEST.h5") # this is the current best model! 
# Note: The runs without added uncertainty in d15N (L391) were used to define the best bootstrap
# The runs with added uncertainty were used to defined the uncertainty; around 0.4 permil higher than without

# test the best model
y_pred = smodel.predict(x=soil15N_test[:,1:], batch_size=bs, verbose=2)
y_true = soil15N_test[:,0]
fig, ax = plt.subplots(2,1,figsize=(8,8))
ax[0].scatter(y_pred,y_true,color="royalblue")
ax[0].plot(y_pred,y_pred,color="navy")
ax[0].set_xlabel("modelled d15N")
ax[0].set_ylabel("observed d15N")
ax[1].scatter(soil15N_test[:,2],y_true,color="red")
ax[1].scatter(soil15N_test[:,2],y_pred,color="royalblue")
ax[1].set_ylabel("MAP")
ax[1].set_ylabel("d15N")
ax[1].legend((["true","pred"]),fontsize=5)
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/test_predictions.png") 
plt.savefig("figs/isoscape_submodel/test_predictions.pdf") 
fig.show() 

# stats
wrmse = ( np.nansum(w_test*((y_pred-y_true)**2)[:,0])/np.nansum(w_test) )**0.5 # weighted RMSE
rmse = (np.nanmean((y_pred-y_true)**2))**0.5 # unweighted RMSE
print("w-rmse:",str(wrmse)," and rmse:",str(rmse))
lm_model = LinearRegression().fit(y_pred,y_true,sample_weight=w_test)
print('intercept:', lm_model.intercept_)
print('slope:', lm_model.coef_[0])
print('R2:', lm_model.score(y_pred,y_true,sample_weight=w_test))

### final model (12.2.21) = -0.41, 1.05, 0.44, 2.57 for int, slope, R2, RMSE for training+val data
### final model (4.8.21) = -0.08, 1.04, 0.41, 2.60 for int, slope, R2, RMSE for test data
### final model (4.8.21) = -0.08, 1.04, 0.41, 2.60 for int, slope, R2, RMSE for test data

# Calculate global gridded d15N with best model
d15N_grid_calc_ML = MAT_grid.copy()*np.nan
predforgrid = predforgrid = smodel.predict(x=inputmatrix_scale, batch_size=bs, verbose=2)
d15N_grid_calc_ML[gridtovec] = predforgrid[:,0]
#d15N_grid_calc_ML = d15N_grid_calc_ML.clip(-5,15) # clip if needed
plot_map(LON,LAT,d15N_grid_calc_ML,"d15N grid, final model",cmap="viridis",filename="figs/isoscape_submodel/final_d15N_isoscape")

# plot without function to set both points and contour same colour range
fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_title("soil d15N, all data (permil)")
cont = plt.contourf(LON,LAT,d15N_grid_calc_ML,cmap="viridis")
points = plt.scatter(soil15N[:,1],soil15N[:,0],c=soil15N[:,5],s=15,edgecolor="black",linewidth=0.5,cmap="viridis")
cbar = plt.colorbar(cont,fraction=0.016, pad=0.04) # add colour bar
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/final_d15N_isoscape_wPoints.png") 
plt.savefig("figs/isoscape_submodel/final_d15N_isoscape_wPoints.pdf") 
fig.show() 

#%% Look at sensitivity of the model to different parameters

# 1. set all other parameters to means to find the range induced by a particular param
fig, ax = plt.subplots(4,3,figsize=(8,8))
full_x = np.concatenate((soil15N_trval[:,1:],soil15N_test[:,1:]), axis=0)
full_y = np.concatenate((soil15N_trval[:,0],soil15N_test[:,0]), axis=0)
full_x_means = full_x.copy()
inducedrange = np.zeros((len(xcols)+1,4)) # for each param save the min,max,std,%var with all other params at mean
for i in range(0,full_x.shape[1]): # x dataset with all values replaced with colmeans (basically 0 due to scaling)
    full_x_means[:,i] = np.nanmean(full_x[:,i])
for n in range(0,len(xcols)+1):
    i = round((n-1)/3)
    j = n % 3
    tmp = full_x.copy()
    if (n != 0): 
        tmp = full_x_means.copy()
        tmp[:,n-1] = full_x[:,n-1]
    full_y_mod = smodel.predict(x=tmp, batch_size=bs, verbose=2)
    ax[i,j].plot(full_y_mod,full_y,"cx")
    ax[i,j].plot(full_y_mod,full_y_mod,"b-")
    if n==0:
        ax[i,j].set_xlabel("Modelled d15N")
        ax[i,j].set_ylabel("Observed d15N")
    else:
        ax[i,j].set_xlabel(xcols[n-1])
    # get stats
    inducedrange[n,0] = np.min(full_y_mod)
    inducedrange[n,1] = np.max(full_y_mod)
    inducedrange[n,2] = np.nanstd(full_y_mod)
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/sensitivity_otherparamsatmean.png") 
plt.savefig("figs/isoscape_submodel/sensitivity_otherparamsatmean.pdf") 
fig.show() 
inducedrange[1:,3] = inducedrange[1:,2]/sum(inducedrange[1:,2])*100 # min,max,std,%var with all other params at mean
inducedrange = pd.DataFrame(inducedrange)
inducedrange.columns=["min","max","std","%var"]
inducedrange.index = ["fullmodel"]+list(xcols)
inducedrange

#%% Look at sensitivity of the model to different parameters

# 2. Shuffle a particular param to see how well model can fit without it
fig, ax = plt.subplots(4,3,figsize=(8,8))
full_x = np.concatenate((soil15N_trval[:,1:],soil15N_test[:,1:]), axis=0)
full_y = np.concatenate((soil15N_trval[:,0],soil15N_test[:,0]), axis=0)
sensitivity = np.zeros((len(xcols)+1,5)) # for each param save the intercept, slope, R2, RMSE
for n in range(0,len(xcols)+1):
    i = round((n-1)/3)
    j = n % 3
    tmp = full_x.copy()
    if (n != 0): 
        shuffler = random.uniform(0,1,len(full_x)).argsort()
        tmp[:,n-1] = full_x[shuffler,n-1]
    full_y_mod = smodel.predict(x=tmp, batch_size=bs, verbose=2)
    ax[i,j].plot(full_y_mod,full_y,"cx")
    ax[i,j].plot(full_y_mod,full_y_mod,"b-")
    if n==0:
        ax[i,j].set_xlabel("Modelled d15N")
        ax[i,j].set_ylabel("Observed d15N")
    else:
        ax[i,j].set_xlabel(xcols[n-1])
    # get stats
    lm_model = LinearRegression().fit(full_y_mod,full_y)
    sensitivity[n,0] = lm_model.intercept_
    sensitivity[n,1] = lm_model.coef_[0]
    sensitivity[n,2] = lm_model.score(full_y_mod,full_y)
    sensitivity[n,3] = (mean((full_y_mod-full_y)**2))**0.5
fig.tight_layout()
plt.savefig("figs/isoscape_submodel/sensitivity_paramshuffled.png") 
plt.savefig("figs/isoscape_submodel/sensitivity_paramshuffled.pdf") 
fig.show() 
sensitivity[1:,4] = sensitivity[1:,3]-sensitivity[0,3]
sensitivity = pd.DataFrame(sensitivity)
sensitivity.columns=["intercept","slope","R2","RMSE","RMSE_reduction"]
sensitivity.index = ["fullmodel"]+list(xcols)
sensitivity

#%% Get WFPS data also (need this for the N2O:N2:NO partitioning in the model)
    
WFPS = np.genfromtxt(fname = "data/Bai2012_Files/Other/WFPS.txt",delimiter=",",skip_header=1) 
# X,X,lon,lat,...,WFPSmean
WFPS_grid = griddata((WFPS[:,2],WFPS[:,3]),WFPS[:,-1], (LON,LAT), method='linear') # lin gives best results!
WFPS_grid = (WFPS_grid*ocean_mask).clip(0,1)
plot_map(LON,LAT,WFPS_grid,"WFPS",filename="figs/input_data/WFPS")

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/climate_input_data.nc','w','NETCDF4'); # using netCDF3 for output format 
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
soild15N = ncout.createVariable('soild15N','f4',('lat','lon'))
soild15N.setncattr('units','permil')
soild15N[:,:] = d15N_grid_calc_ML
soild15N_BS = ncout.createVariable('soild15N_BS','f4',('lat','lon'))
soild15N_BS.setncattr('units','permil')
soild15N_BS[:,:] = meangrid
soild15N_BSUnc = ncout.createVariable('soild15N_BSUnc','f4',('lat','lon'))
soild15N_BSUnc.setncattr('units','permil')
soild15N_BSUnc[:,:] = stdgrid
ncout.close()

# test the file!
f = nc4.Dataset('data/climate_input_data.nc','r')
plot_map(LON,LAT,f.variables["soilC"][:,:],"globalC_grid")
plot_map(LON,LAT,f.variables["MAPrecip"][:,:],"MAP")
plot_map(LON,LAT,f.variables["MATemp"][:,:],"MAT")
plot_map(LON,LAT,f.variables["soilpH"][:,:],"pH")
plot_map(LON,LAT,f.variables["soilN"][:,:],"N")
plot_map(LON,LAT,f.variables["soild15N"][:,:],"d15N")
plot_map(LON,LAT,f.variables["soild15N_BS"][:,:],"d15N, Bootstrap")
plot_map(LON,LAT,f.variables["soild15N_BSUnc"][:,:],"d15N, BS. Uncertainty")
plot_map(LON,LAT,f.variables["WFPS"][:,:],"WFPS")
