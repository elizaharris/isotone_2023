#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Set up to run the IsoTONE model (data imports, modules, etc.)

# Author: Eliza Harris
# Created on Tue Mar  3 09:09:23 2020

# Tasks:
# Import all model input data, some checking and plotting
# Import the observational data for model constraint
# Set the model parameters, some of which will be optimized in the MCMC

#%% Preamble 

# basic packages 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import xarray as xr
from datetime import datetime as dt
import cartopy.crs as ccrs
import os
import netCDF4 as nc4
import time
import datetime

# Get the configs
import configs as configs

# my functions for this model
from src.utils import plot_map
from src.utils import climzone_means
import src.parameterisations_v3 as para # partition nitrif and denitrif
plt.close("all")

# Set up time range
years = np.arange(configs.res["ystart"],configs.res["yend"])

# set resolution and grid for model (standard res, may be overwritten later for arctic)
resolution = configs.res["ystart"]
lat_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lon_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

# find the area of grid cells in m2
lats = np.deg2rad(np.arange(configs.res["latstart"]-0.5*configs.res["latres"], configs.res["latend"]+0.5*configs.res["latres"], configs.res["latres"])) # Calculations needs to be in radians
r_sq = 6371000**2 # earth's radius in m, squared
n_lons = int(360./configs.res["latres"]) # how many longitudes 
area_grid = r_sq*np.ones(n_lons)[:, None]*np.deg2rad(configs.res["latres"])*(np.sin(lats[1:]) - np.sin(lats[:-1])) # area in m2
area_grid = np.transpose(area_grid)


#%% Import and plot the input data
### Note: if the model grid range or resolution is changed, these inputs need to be rerun!
print("Importing input data...")
    
# Import the data processed when created the isoscape
inputs = nc4.Dataset('data/climate_input_data.nc','r')
C_grid = inputs.variables["soilC"][:,:]; plot_map(LON,LAT,C_grid,"soil C (kg m-2)")
MAP_grid = inputs.variables["MAPrecip"][:,:]; plot_map(LON,LAT,MAP_grid,"MAP (mm)")
MAT_grid = inputs.variables["MATemp"][:,:]; plot_map(LON,LAT,MAP_grid,"MAT (degC)",vminmax=(0,3000))
pH_grid = inputs.variables["soilpH"][:,:]; plot_map(LON,LAT,pH_grid,"pH")
N_grid = inputs.variables["soilN"][:,:]; plot_map(LON,LAT,N_grid,"N (g m-2)")
AI_grid = inputs.variables["AridIndex"][:,:]; plot_map(LON,LAT,AI_grid,"Aridity Index")
BulkD_grid = inputs.variables["BulkD"][:,:]; plot_map(LON,LAT,BulkD_grid,"Bulk Density (g cm-3)")
WFPS_grid = inputs.variables["WFPS"][:,:]*100; plot_map(LON,LAT,WFPS_grid,"WFPS (%)")
d15N_grid_new = inputs.variables["soild15N"][:,:]; plot_map(LON,LAT,d15N_grid_new,"d15N")
d15Nerr_grid_new = inputs.variables["soild15N_BSUnc"][:,:]; plot_map(LON,LAT,d15Nerr_grid_new,"d15N uncertainty")

# For d15N, use an old dataset for now...
inputs = nc4.Dataset('data/climate_input_data_old.nc','r') 
d15N_grid = d15N_grid_new.copy()*np.nan
d15N_grid[0:290,0:720] = inputs.variables["soild15N"][:,:]
plot_map(LON,LAT,d15N_grid,"d15N")
d15Nerr_grid = d15Nerr_grid_new.copy()*np.nan
d15Nerr_grid[0:290,0:720] = inputs.variables["soild15N_BSUnc"][:,:]
plot_map(LON,LAT,d15Nerr_grid,"d15N uncertainty")

plt.close("all")

# Import estimates of the fraction of N volatilized as NH3, from Bai 2012, processed in src/regrid_fNH3_fromBai
f = nc4.Dataset('data/fNH3_data.nc','r')
fNH3_grid = f.variables["fNH3"][:,:]; plot_map(LON,LAT,fNH3_grid,"fNH3")

# Import EDGAR gridded emission estimates, processed in src/regrid_EDGAR
# Note: In Arctic version, separate files created for each cat! This is to better manage computations with hgih res
f = nc4.Dataset('data/largeData/EDGAR_data_extrap.nc','r')
AGS_grid = f.variables["AGS"][:,:,:]; plot_map(LON,LAT,AGS_grid[150,:,:],"EDGAR agriculture (g N2O-N/m2/year): 2000")
TRO_grid = f.variables["TRO"][:,:,:]; plot_map(LON,LAT,TRO_grid[150,:,:],"EDGAR transport (g N2O-N/m2/year): 2000")
IDE_grid = f.variables["IDE"][:,:,:]; plot_map(LON,LAT,IDE_grid[150,:,:],"EDGAR indirect (g N2O-N/m2/year): 2000")
WWT_grid = f.variables["WWT"][:,:,:]; plot_map(LON,LAT,WWT_grid[150,:,:],"EDGAR wastewater (g N2O-N/m2/year): 2000")
CHE_grid = f.variables["IDE"][:,:,:]; plot_map(LON,LAT,CHE_grid[150,:,:],"EDGAR chem industry (g N2O-N/m2/year): 2000")
ENE_grid = f.variables["WWT"][:,:,:]; plot_map(LON,LAT,ENE_grid[150,:,:],"EDGAR energy (g N2O-N/m2/year): 2000")

# Import N inputs (fixation, deposition, fertilisation), processed in src/regrid_N_inputs
Ninputs = nc4.Dataset('data/largeData/N_input_data.nc','r')
datarange = Ninputs.variables["time"][:].data

# Import T anomalies (annual, Hadley centre), processed in src/regrid_T_anomalies
Tanom = nc4.Dataset('data/largeData/Tanom_data.nc','r')
Tnorm = Tanom.variables["Tnorm"][:,:,:].data
plot_map(LON,LAT,Tnorm[150,:,:],"Normalised temp anomaly: "+str(Tanom.variables["time"][150]),filename="figs/input_data/T_anomalies_norm")
mean_T_rise = Tanom.variables["mean_T_rise"][:].data

plt.close("all")

#%% Observational data for model comparison
print("Importing obs data...")

# Atmospheric time series data (N2O ppb and isotopes)
# Combined from various sources, see data/atmos_data/combine_atmos_data.py for details
N2O_atmos = pd.read_csv("data/atmos_data/SummAtmosData.csv", sep=',') 

# flux data from Chris Dorich, divided into climate regions
if 1: # Sometimes don't run this; for Arctic should have Arctic fluxes, and matching this detailed grid is too slow
    N2O_fluxes = pd.read_csv("data/N2O_ChrisDorich/globaln2o_sites_filled.csv")
    # cut to: MAT, MAP, BulkD, C, N, pH, lat, lon, flux, flux_se, EF, fert rate for easy comparison
    # fix NAs and string values also
    def fix_nas(data,NA_value): # function to make different nans to true nans
        for n in range(0,len(NA_value)):
            tmp = np.where(new == NA_value[n])
            data[tmp] = np.nan
        return(data)
    def fix_strings(data,chars):
        for n in range(0,len(chars)):
            tmp = np.where(np.char.find(data.astype("str"),chars[n])!=-1)
            data[tmp] = np.nan
        return(data)
    N2O_flux_params = (("Temp_C","Precip_mm","Soil_BD_gcm3","Soil_SOC_perc","Soil_SON_perc","pH","Lat","Long",'N2O_kgN2O-N_ha', 'N2O_se_kgN2O-N_ha',"EF",'N_App_Rate_kgN_ha'))
    N2O_fluxes_short = np.array(N2O_fluxes[N2O_flux_params[0]])
    for n in range(1,len(N2O_flux_params)):
        new = np.array(N2O_fluxes[N2O_flux_params[n]])
        new = fix_nas(new,("*","na","nan"))
        new = fix_strings(new,(">","-"))
        new = new.astype(np.float64)
        N2O_fluxes_short = np.vstack((N2O_fluxes_short,new))
    N2O_fluxes_short = N2O_fluxes_short.transpose()
    tmp = np.where(~np.isnan(N2O_fluxes_short[:,8])) # Flux is NA here...
    # Match the anc data to each flux point
    flux_ancdata = np.zeros((N2O_fluxes.shape[0],8))*np.nan # MAT, MAP, BulkD, C, N, pH, AI, d15N
    for n in range(0,N2O_fluxes.shape[0]):
        if (~np.isnan(N2O_fluxes["Lat"][n])) & (~np.isnan(N2O_fluxes["Long"][n])):
            r = np.where(np.nanmin(abs(LAT[:,0] - N2O_fluxes["Lat"][n])) == abs(LAT[:,0] - N2O_fluxes["Lat"][n])) # match flux lat long to gridcells
            c = np.where(np.nanmin(abs(LON[0,:] - N2O_fluxes["Long"][n])) == abs(LON[0,:] - N2O_fluxes["Long"][n]))
            if len(r[0])>1: r=r[0][0] # if two grid cells match lat/lon equally, take first
            if len(c[0])>1: c=c[0][0]
            flux_ancdata[n,0] = MAT_grid[r,c]
            flux_ancdata[n,1] = MAP_grid[r,c]
            flux_ancdata[n,2] = BulkD_grid[r,c]
            flux_ancdata[n,3] = C_grid[r,c]/10 # change to %
            flux_ancdata[n,4] = N_grid[r,c]/10 # change to %
            flux_ancdata[n,5] = pH_grid[r,c]
            flux_ancdata[n,6] = AI_grid[r,c]
            flux_ancdata[n,7] = d15N_grid[r,c]
    N2O_fluxes_short = N2O_fluxes_short[tmp[0],:]
    flux_ancdata = flux_ancdata[tmp[0],:]
    # calc mean EFs by climate zone      
    ### If this changes, edit analogous function for model data in 3_Final_RunMCMC!
    N2O_fluxes_zones = climzone_means(var1_grid = MAT_grid, var2_grid = MAP_grid, datavar1 = flux_ancdata[:,0],
                                    datavar2 = flux_ancdata[:,1],data = N2O_fluxes_short[:,10],LON=LON,LAT=LAT,bins=4,plotfigs="Y")
    # final output
    print(str(N2O_fluxes_short.shape[0])+" flux measurements from "+str(len(set(N2O_fluxes["Location"][tmp[0]])))+" locations")

plt.close("all")

#%% Collect input variables that may be optimised
print("Collecting input variables...")
 
# Dates for model output, to be able to compare to observations in df N2O_atmos (these are where values are averaged to in combine_atmos_data)
n_windows = len(configs.mod_obs_time["ystart"])
for n in np.arange(0,n_windows):
    tmp = np.arange(configs.mod_obs_time["ystart"][n],configs.mod_obs_time["yend"][n],configs.mod_obs_time["yres"][n])
    if n == 0:
        tstarts = tmp
    else:
        tstarts = np.hstack((tstarts,tmp))

# 1. For soil model
# Input vars: Fractionation and partitioning 
# Uncertainty windows in form of (value, low, high) for uniform windows or (value, sd) for gaussian windows
fracex = 0.4
d15N_inat = -0.5 # inputs d15N; lower = more gas loss 
E_NH3 = -17.9 # d15N frac for amm vol losses; lower = less gas loss 
E_L	= [-1,0,5] # d15N frac for leaching from Bai et al.; higher = less gas loss; model breaks down with much deviation
E_nit = [-56.6,7.3] # d15N frac for "exiting" nitrate pool for nitrification; lower = loss gas loss
E_denit	= [-31.3,6.1] # d15N frac for "exiting" nitrate pool for denitification (NO3- to NO2-); lower = loss gas loss
E_denit_N2O = [-14.9,6.7] # d15N frac factor for N2O production during denitrification (NO2- to N2O); minor impact on d15N emitted
E_red = [-6.6,2.7] # d15N frac factor for N2O reduction; minor impact on d15N emitted
E_SP_nit = [29.9,2.9] # SP endmember for nitrification
E_SP_denit = [-1.6,3.0] # SP endmember for denitrification
E_SP_red = [-5,-8,-2] # SP frac factor for N2O reduction
f_soilpoolloss = 0.5 # fraction of soil pool lost to atm in each timestep; no impact on results, just a model feature
scale_fitNO = np.array((1.,1.,1.,1.)) # scaling factors for the four parameters a,b,c,d of the sigmoid fit for NO/(N2O+NO) where: a + d = RHS plateau; b = mid point of rise; c = slope (larger = steeper); d = LHS plateau
scale_fitN2 = np.array((1.,1.,1.,1.)) # scaling factors for N2/(N2O+N2) sigmoid fit
scale_fitnit = np.array((1.,1.,1.,1.)) # scaling factors for nit/(nit+denit) sigmoid fit
#gaspart = para.gp_gaspart(20,scale_fitNO = scale_fitNO,scale_fitN2 = scale_fitN2,plotYN="Y")
N2Opart = para.gp_nitdenit(20,scale_fitnit = scale_fitnit,plotYN="Y")

# 2. For emissions
fertEFred = 0.4 #1 # reduction in EF for N2O emitted from fertiliser relative to natural conditions 
scaleDep = 1 # scale for deposition N inputs
scaleFix = 1 # scale for fixation N inputs
temp_sens = 1.1 # increase in emissions for 1 degree of warming
# wwt source
d15_WWT = [-11.6,12.7] # harris et al, 2017, SI synthesis
SP_WWT = [10.5,5.7] # harris et al, 2017, SI synthesis
# tro source
d15_TRO = [-7.2,1.2] # harris et al, 2017, SI synthesis
SP_TRO = [10.0,4.3] # harris et al, 2017, SI synthesis
# che source
d15_CHE = [-8.3,10.6] # harris et al, 2017, SI synthesis
SP_CHE = [3.3,5.5] # harris et al, 2017, SI synthesis
# ene source
d15_ENE = [3.9,2.9] # harris et al, 2017, SI synthesis
SP_ENE = [17.6,0.5] # harris et al, 2017, SI synthesis

# 3. For atm model
# atmospheric vars
MW_air = 28.9647 # molecular weight of air in g/mol
MW_N2O = 44 # molecular weight of N2O in g/mol
T_S = np.array([5.37e17,1.26e17]) # Troposphere-stratosphere exchange in kg yr-1; from Schilt 2014
T_Sm = T_S*1000/MW_air # mol/year
m_trop = 0.85*1.77e20; # moles of gas in the troposphere
m_strat = 0.15*1.77e20; # moles of gas in the stratosphere 
F_ocean = [5.1,1.8] # Ocean N2O flux in Tg N y-1 (from Tian 2020)
ltPD = [116,9] # N2O lifetime in years for PD
ltPD_ltPI = [1.06,0.02]
ltPI = [ltPD[0]*ltPD_ltPI[0],ltPD[1]*(ltPD_ltPI[1]+1)] # N2O lifetime in years for PI (from Prather et al; 118 in Sowers et al.)
# preanth trop
c_prea = [270,7.5] # preanthropogenic N2O concentration (ppb) - 260-270
d15_prea = [8.9,2.0] # del values for the preanth troposphere, from Toyoda et al. 2013 - average from Bernard, Rockmann, Sowers
d18_prea = [46.1,2.0] 
SP_prea = [19.05,2.0]
# ocean source
d15_ocean = [5.1,1.9] # del values for the ocean source, from Schilt et al.
d18_ocean = [44.8,3.6] # del values for the ocean source, from Schilt et al.
SP_ocean = [15.8,7.1] # SP values from Snider dataset ("marine")

plt.close("all")

