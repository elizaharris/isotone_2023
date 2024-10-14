#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:50:35 2020

This file imports and regrids all N inputs (dep, fix, fert) and saves the results for
1850 - 2050 to save modelling time.

@author: elizaharris


"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import xarray as xr
import cartopy.crs as ccrs
from datetime import datetime as dt
from scipy.interpolate import griddata
import os

parentdir = os.getcwd()
os.chdir(parentdir+"/isotone_arcticBranch") # Change to the isotone repo as working directory
from src.utils import plot_map, nc_to_list
import configs as configs

# set resolution and grid for model
resolution = configs.res["ystart"]
lat_out = np.arange(configs.res["latstart"], configs.res["latend"], configs.res["latres"]) # min max spacing
lon_out = np.arange(configs.res["lonstart"], configs.res["lonend"], configs.res["lonres"])
(LON,LAT) = np.meshgrid(lon_out,lat_out)

# Set years for this regridding and create space for output
years = np.arange(1850,2050)
fix_griddata = np.zeros((len(years),LAT.shape[0],LON.shape[1])) 
dep_griddata = np.zeros((len(years),LAT.shape[0],LON.shape[1]))
fert_griddata = np.zeros((len(years),LAT.shape[0],LON.shape[1]))

startyear_prev = 0
datayear_prev = 0
datayear_d_prev = 0
for n in range(0,years.shape[0]) :
    print(years[n])
    
    ### Fixation data from CABLE from Ying Ping Wang, CSIRO (A1 scenario)
    fix_datafreq = 10 # years per netcdf file
    fix_datarange = np.arange(1901,2100,fix_datafreq) # years that each file starts
    # find the year to use for fixation data - nearest covered by the dataset
    if (years[n]<min(fix_datarange)) : datayear = min(fix_datarange)
    elif (years[n]>max(fix_datarange)) : datayear = max(fix_datarange)
    else: datayear = years[n]
    # find the name of the dataset file and import
    startyear = math.floor((datayear-1)/fix_datafreq)*10+1
    if (startyear != startyear_prev) :
        filename = "../isotone-rawdata/YPWang_N_fixation_CABLE_v2/CABLE_bnf_A1_annual_2D_"+str(startyear)+"-"+str(startyear+fix_datafreq-1)+".nc"    
        fixation = xr.open_dataset(filename)
    if (datayear != datayear_prev) :
        tmp = fixation.bnf.sel(year=(datayear-startyear)).values.squeeze() 
        # plot_map(fixation.lon, fixation.lat, tmp, "fixation: "+str(years[n])) # plot original if needed
        # regrid the data 
        fixation_list = nc_to_list(fixation.lon, fixation.lat, tmp)
        fixation_grid = griddata((fixation_list[:,0], fixation_list[:,1]), fixation_list[:,2], (LON,LAT), method='linear')
    fix_griddata[n,:,:] = fixation_grid

    ### Deposition data from CABLE from Ying Ping
    dep_datafreq = 1 # years per netcdf file
    dep_datarange = np.arange(1860,2051,dep_datafreq) # years that each file starts
    # find the year to use for fixation data - nearest covered by the dataset
    if (years[n]<min(dep_datarange)) : datayear_d = min(dep_datarange)
    elif (years[n]>max(dep_datarange)) : datayear_d = max(dep_datarange)
    else: datayear_d = years[n]
    # find the name of the dataset file and import
    filename = "../isotone-rawdata/YPWang_N_deposition/ndep_"+str(datayear_d)+"_0.5x0.5.nc"
    if (datayear_d != datayear_d_prev) :
        deposition = xr.open_dataset(filename)
        tmp = deposition.N_total_deposition.values.squeeze() 
        # plot_map(deposition.longitude, deposition.latitude, tmp, "N deposition (orig, g N/m2/year)")
        # regrid the data 
        deposition_list = nc_to_list(deposition.longitude, deposition.latitude, tmp)
        deposition_grid = griddata((deposition_list[:,0], deposition_list[:,1]), deposition_list[:,2], (LON,LAT), method='linear')
    dep_griddata[n,:,:] = deposition_grid    

    ### Fertiliser - from https://luh.umd.edu/data.shtml (see docs in folder for more info)
    opened = 0
    if (years[n]<=2015) : # historical fert dataset v2h
        if (opened==0) : # gives fert intensity mean per cropland area
            filename = "../isotone-rawdata/management/management.nc"    
            fert = xr.open_dataset(filename,decode_times=False) 
            opened = 1 + opened
        if (opened<2) : # gives cropland fraction in a gridcell
            filename = "../isotone-rawdata/management/states.nc"    
            states = xr.open_dataset(filename,decode_times=False) 
            opened = 1 + opened
        fert1 = fert.fertl_c3ann.sel(time=years[n]-850).values.squeeze()*states.c3ann.sel(time=years[n]-850).values.squeeze()
        fert2 = fert.fertl_c4ann.sel(time=years[n]-850).values.squeeze()*states.c4ann.sel(time=years[n]-850).values.squeeze()
        fert3 = fert.fertl_c3per.sel(time=years[n]-850).values.squeeze()*states.c3per.sel(time=years[n]-850).values.squeeze() 
        fert4 = fert.fertl_c4per.sel(time=years[n]-850).values.squeeze()*states.c4per.sel(time=years[n]-850).values.squeeze() 
        fert5 = fert.fertl_c3nfx.sel(time=years[n]-850).values.squeeze()*states.c3nfx.sel(time=years[n]-850).values.squeeze() 
        fert_total = fert1+fert2+fert3+fert4+fert5 # kg N ha-1 y-1
        fert_total = fert_total*1000/10000
        fert_list = nc_to_list(fert.lon, fert.lat, fert_total)
        fert_grid = griddata((fert_list[:,0], fert_list[:,1]), fert_list[:,2], (LON,LAT), method='linear')
    opened1 = 0
    if (years[n]>2015) : # predicted fert dataset - note the 2015 grids for hist/pred match perfectly!!
        if (opened1==0) :
            filename = "../isotone-rawdata/management/multiple-management_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp534-2-1-f_gn_2015-2100.nc"    
            fert = xr.open_dataset(filename,decode_times=False) 
            opened1 = 1 + opened1
        if (opened1<2) : # gives cropland fraction in a gridcell
            filename = "../isotone-rawdata/management/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-GCAM-ssp434-2-1-f_gn_2015-2100.nc"    
            states = xr.open_dataset(filename,decode_times=False) 
            opened1 = 1 + opened1
        fert1 = fert.fertl_c3ann.sel(time=years[n]-2015).values.squeeze()*states.c3ann.sel(time=years[n]-2015).values.squeeze()
        fert2 = fert.fertl_c4ann.sel(time=years[n]-2015).values.squeeze()*states.c4ann.sel(time=years[n]-2015).values.squeeze()
        fert3 = fert.fertl_c3per.sel(time=years[n]-2015).values.squeeze()*states.c3per.sel(time=years[n]-2015).values.squeeze() 
        fert4 = fert.fertl_c4per.sel(time=years[n]-2015).values.squeeze()*states.c4per.sel(time=years[n]-2015).values.squeeze() 
        fert5 = fert.fertl_c3nfx.sel(time=years[n]-2015).values.squeeze()*states.c3nfx.sel(time=years[n]-2015).values.squeeze() 
        fert_total = fert1+fert2+fert3+fert4+fert5 # kg N ha-1 y-1
        fert_total = fert_total*1000/10000
        fert_list = nc_to_list(fert.lon, fert.lat, fert_total)
        fert_grid = griddata((fert_list[:,0], fert_list[:,1]), fert_list[:,2], (LON,LAT), method='linear')
    fert_griddata[n,:,:] = fert_grid    
  
    startyear_prev = startyear
    datayear_prev = datayear
    datayear_d_prev = datayear_d

# Save plots for some years
n = 70 # About 1920
plot_map(LON,LAT,fix_griddata[n,:,:],"N fixation (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))
plot_map(LON,LAT,dep_griddata[n,:,:],"N deposition (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))
plot_map(LON,LAT,fert_griddata[n,:,:],"N fertilization (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))
n = 170 # About 2020
plot_map(LON,LAT,fix_griddata[n,:,:],"N fixation (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))
plot_map(LON,LAT,dep_griddata[n,:,:],"N deposition (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))
plot_map(LON,LAT,fert_griddata[n,:,:],"N fertilization (g N/m2/year): "+str(years[n]),filename="figs/input_data/fix_"+str(years[n]))

#%% save the data as a netcdf file

import netCDF4 as nc4

ncout = nc4.Dataset('data/largeData/N_input_data.nc','w','NETCDF4'); # using netCDF3 for output format 
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
fertilisation = ncout.createVariable('fertilisation','f4',('time','lat','lon'))
fertilisation.setncattr('units','g N/m2/year')
fertilisation.setncattr('source','Land use harmonization datasets for past and future combined (https://luh.umd.edu/data.shtml)')
fertilisation[:,:,:] = fert_griddata[:,:,:]
deposition = ncout.createVariable('deposition','f4',('time','lat','lon'))
deposition.setncattr('units','g N/m2/year')
deposition.setncattr('source','CABLE model - data provided by Ying Ping Wang, CSIRO, February 2020')
deposition[:,:,:] = dep_griddata[:,:,:]
fixation = ncout.createVariable('fixation','f4',('time','lat','lon'))
fixation.setncattr('units','g N/m2/year')
fixation.setncattr('source','CABLE model - data provided by Ying Ping Wang, CSIRO, February 2020')
fixation[:,:,:] = fix_griddata[:,:,:]
ncout.close()

# test the file!
f = nc4.Dataset('data/largeData/N_input_data.nc','r')
plot_map(LON,LAT,f.variables["fertilisation"][n,:,:],"N fertilization (g N/m2/year): "+str(f.variables["time"][n]))
plot_map(LON,LAT,f.variables["deposition"][n,:,:],"N deposition (g N/m2/year): "+str(f.variables["time"][n]))
plot_map(LON,LAT,f.variables["fixation"][n,:,:],"N fixation (g N/m2/year): "+str(f.variables["time"][n]))
