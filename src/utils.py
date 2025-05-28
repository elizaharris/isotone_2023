#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 09:19:12 2021

@author: elizaharris
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import rasterio as rio
import cartopy.crs as ccrs
import scipy.stats 

# Uniform probability: 0 if within range, 1 if no
def prob_u(x,lo,hi): return float( (x < hi).all() and (x > lo).all() )

# Gaussian probability function
def prob_g(x,x0,sigma): 
    tmp = (x-x0)/sigma # Normalise each value
    res = scipy.stats.norm(0, 1).pdf(tmp) # Find probability density of each value in Gaussian
    return np.nanmean(res)

# Export model results if desired
def export_model_results(fullmod,suffix,run_name): 
    fname = "outputs/"+run_name+'/isotone_output'+suffix+'.nc'
    if os.path.exists(fname):
        os.remove(fname)
    ncout = nc4.Dataset(fname,'w','NETCDF4'); # using netCDF3 for output format 
    ncout.createDimension('obs_dim',fullmod_post[0].shape[0])
    ncout.createDimension('lon',preamble.LON.shape[1])
    ncout.createDimension('lat',preamble.LAT.shape[0])
    ncout.createDimension('datarng_len',len(datarng[0]))
    ncout.createDimension('datarng_dim',2)
    ncout.createDimension('time',len(preamble.years))
    ncout.createDimension('soilmodres_dim',fullmod_post[3].shape[1])
    ncout.createDimension('N_processes_vec_dim',fullmod_post[4].shape[1])
    ncout.createDimension('N_summary_dim',fullmod_post[5].shape[1])
    ncout.createDimension('N_summary_full_dim',fullmod_post[6].shape[1])
    lonvar = ncout.createVariable('lon','f4',('lon'))
    lonvar[:] = preamble.LON[0,:]
    lonvar.setncattr('units','degrees east')
    latvar = ncout.createVariable('lat','f4',('lat'))
    latvar[:] = preamble.LAT[:,0]
    latvar.setncattr('units','degrees north')
    timevar = ncout.createVariable('time','f4',('time'))
    timevar[:]= preamble.years
    x = np.zeros((len(datarng[0]),2)) # Convert shape of datarng to array
    x[:,0] = datarng[0]; x[:,1] = datarng[1]
    datarngvar = ncout.createVariable('datarng','f4',('datarng_len','datarng_dim'))
    datarngvar[:,:] = x[:,:]
    datarngvar.setncattr('description','positions on lat-lon grid with non-nan values')
    model_res = ncout.createVariable('model_res','f4',('obs_dim'))
    model_res.setncattr('units','Varying; see linked observations')
    model_res.setncattr('source','fullmod[0]')
    model_res[:] = fullmod[0]
    soilmodres = ncout.createVariable('soilmodres','f4',('datarng_len','soilmodres_dim'))
    soilmodres.setncattr('units','Varying; see model')
    soilmodres.setncattr('source','fullmod[3]')
    soilmodres[:,:] = fullmod[3]
    N_processes_vec = ncout.createVariable('N_processes_vec','f4',('datarng_len','N_processes_vec_dim'))
    N_processes_vec.setncattr('units','Varying; see model')
    N_processes_vec.setncattr('source','fullmod[4]')
    N_processes_vec[:,:] = fullmod[4]
    N_summary = ncout.createVariable('N_summary','f4',('time','N_summary_dim'))
    N_summary.setncattr('units','Varying; see model')
    N_summary.setncattr('source','fullmod[5]')
    N_summary[:,:] = fullmod[5]
    N_summary_full = ncout.createVariable('N_summary_full','f4',('time','N_summary_full_dim'))
    N_summary_full.setncattr('units','Varying; see model')
    N_summary_full.setncattr('source','fullmod[6]')
    N_summary_full[:,:] = fullmod[6]
    ncout.close()
    # test the file!
    f = nc4.Dataset(fname,'r')
    a_nc = np.nanmean(f.variables['N_summary_full'][:,:])
    a = np.nanmean(fullmod[6][:,:])
    b_nc = np.nanmean(f.variables['N_processes_vec'][:,:])
    b = np.nanmean(fullmod[4][:,:])
    if (round(a,3)-round(a_nc,3)<0.0001) & (round(b,3)-round(b_nc,3)<0.0001):
        res = "Success"
    else:
        res = "Something went wrong!"
    return res

# define function for plotting    
def plot_map(longi,lati,gridval,title="title",vminmax=(np.nan,np.nan),cmap="viridis",filename="figs/testfig",show=0) :
    fig = plt.figure(figsize=(12,6))
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
    if show==1:
        fig.show() 

# function for mean by "climate zones" defined by 2 variables
def climzone_means(var1_grid, var2_grid, datavar1,datavar2,data,LON,LAT,bins=4,plotfigs="N"):  
    if (bins % 2) != 0: 
        print("bins must be even number!! ")
    else: 
        # create bins based on mean +- standard deviations and calc bin means for data
        #var1_bins_orig = np.linspace(np.nanmean(var1_grid)-bins/2*np.nanstd(var1_grid),np.nanmean(var1_grid)+bins/2*np.nanstd(var1_grid),num=bins+1)
        #var2_bins_orig = np.linspace(0,np.nanmean(var2_grid)+bins/2*np.nanstd(var2_grid),num=bins+1)
        # create bins based on percentiles of the data
        var1_bins_orig = np.percentile(var1_grid[~np.isnan(var1_grid)],(0,25,50,75,100))
        var2_bins_orig = np.percentile(var2_grid[~np.isnan(var2_grid)],(0,25,50,75,100))
        var1_bins = var1_bins_orig.copy(); var1_bins[0] = -1000; var1_bins[-1] = +1000 # make sure edge bins should span full range
        var2_bins = var2_bins_orig.copy(); var2_bins[0] = -1000; var2_bins[-1] = +1000 # edge bins should span full range
        res = np.zeros((bins**2,7))*np.nan # var1 bin #, var1 bin mid, var2 bin #, var2 bin mid, data mean, data std, data n
        l = 0
        for n in range(0,bins):
            for i in range(0,bins):
                res[l,0] = n
                res[l,1] = (var1_bins_orig[n]+var1_bins_orig[n+1])/2
                res[l,2] = i
                res[l,3] = (var2_bins_orig[i]+var2_bins_orig[i+1])/2        
                tmp = np.where((datavar1>=var1_bins[n]) & (datavar1<=var1_bins[n+1]) & (datavar2>=var2_bins[i]) & (datavar2<=var2_bins[i+1]))
                res[l,6] = len(tmp[0])
                if len(tmp[0])>1:
                    res[l,4] = np.nanmean(data[tmp])
                    res[l,5] = np.nanstd(data[tmp])
                l = l+1
        # plot the results
        if plotfigs == "Y":
            fig, ax = plt.subplots(1,1)
            for i in range(0,bins):
                tmp = (res[:,2] == i) & ~np.isnan(res[:,4])
                ax.plot(res[tmp,1],res[tmp,4],"-",c="C"+str(i))
            tmp = (res[:,0] == 0)
            ax.legend(res[tmp,3].astype(str))  
            for i in range(0,bins):
                tmp = (res[:,2] == i) & ~np.isnan(res[:,4])
                ax.scatter(res[tmp,1],res[tmp,4],c="C"+str(i))
                ax.plot(res[tmp,1],res[tmp,4]-res[tmp,5],":",c="C"+str(i))
                ax.plot(res[tmp,1],res[tmp,4]+res[tmp,5],":",c="C"+str(i))
            for i, txt in enumerate(res[:,6]):
                ax.annotate(int(txt), (res[i,1]+0.01,res[i,4]+0.01))
            ax.set_xlabel("var1")
            ax.set_ylabel("N2O EF")
        # plot the zones if needed
        if plotfigs == "Y":
            climzones = var1_grid.copy()*np.nan
            l = 0
            for n in range(0,bins):
                for i in range(0,bins):
                    tmp = np.where((var1_grid>=var1_bins[n]) & (var1_grid<=var1_bins[n+1]) & (var2_grid>=var2_bins[i]) & (var2_grid<=var2_bins[i+1]))
                    climzones[tmp] = l       
                    l = l+1
            plot_map(LON,LAT,climzones,"Climate zones - n(bins) = "+str(bins))
        return(res)
    
# function to transform to a list of lon, lat, values and average every X points
def rast_to_list(rastdata,Xred,takemean="Y",dataclip="N",replace_w_nan="N") :
    arrdata = rastdata.read(1)
    if dataclip != "N" : arrdata = arrdata.clip(dataclip[0],dataclip[1])
    if replace_w_nan != "N" : arrdata[arrdata==replace_w_nan] = np.nan
    output = np.zeros((round(arrdata.shape[0]/Xred)*round(arrdata.shape[1]/Xred),3))
    z = 0
    for n in range(0,round(arrdata.shape[0]/Xred)) :
        if n % 100 == 0: 
            print("n = "+str(n)+"(total = "+str(arrdata.shape[0]/Xred)+")")
        for i in range(0,round(arrdata.shape[1]/Xred)) :
            #if i % 100 == 0: 
                #print("n = "+str(n)+"; i = "+str(i)+" (total = "+str(arrdata.shape[0]/Xred)+","+str(arrdata.shape[1]/Xred)+")")
            nmid = int(n*Xred+Xred/2)
            imid = int(i*Xred+Xred/2)
            tmp = rio.transform.TransformMethodsMixin.xy(rastdata,nmid,imid)
            output[z,0] = tmp[0]
            output[z,1] = tmp[1]
            if takemean == "Y" :  # take mean of all values within xred window
                nvals = np.repeat(range(int(nmid-Xred/2),int(nmid+Xred/2)),Xred) # repeats each element X times
                ivals = np.tile(range(int(imid-Xred/2),int(imid+Xred/2)),Xred) # repeats whole range X times
            if takemean == "N" :  # take only the single points at each Xred
                nvals = nmid 
                ivals = imid 
            output[z,2] = np.nanmean(arrdata[nvals,ivals])
            z = z+1
    return output

# function to transform to a list of lon, lat, values and average every X points
def nc_to_list(loni,lati,data) :
    output = np.zeros((loni.shape[0]*lati.shape[0],3))
    output[:,0] = np.tile(loni,lati.shape[0])
    output[:,1] = np.repeat(lati,loni.shape[0])
    output[:,2] = data.flatten()
    return output

# function to transform to a list of lon, lat, values and average every X points
def grid_to_list(loni,lati,data) :
    output = np.zeros((loni.shape[0]*lati.shape[0],3))
    output[:,0] = np.tile(loni,lati.shape[0])
    output[:,1] = np.repeat(lati,loni.shape[0])
    output[:,2] = data.flatten()
    return output

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

# Sometimes the accept_MC file has accept-stepsize as a single column, other times as two columns; this function corrects this
def fix_acceptMC(accept_MC_res): 
    if len(accept_MC_res.shape) == 1: 
        new = np.zeros((int(accept_MC_res.shape[0]/2),2))
        new[:,0] = accept_MC_res[0:int(accept_MC_res.shape[0]/2)]
        new[:,1] = accept_MC_res[int(accept_MC_res.shape[0]/2):]
    else : new = accept_MC_res
    return(new)

# find combined mean and sd from the run means and stds, accounting for different numbers of accepted values in each run
def combined_mean(data,std,n) :
    data_dim = len(data.shape)
    Sx = data*0
    for i in range(0,len(n)) :
        if data_dim == 3: Sx[i,:,:] = data[i,:,:]*n[i]
        if data_dim == 2: Sx[i,:] = data[i,:]*n[i]
    tx = np.nansum(Sx,axis=0) # total sum of all values from all runs
    m = tx/np.nansum(n) # combined mean
    Sxsd = data*0   
    for i in range(0,len(n)) :
        if data_dim == 3: Sxsd[i,:,:] = std[i,:,:]*std[i,:,:]*(n[i]-1) + (Sx[i,:,:]*Sx[i,:,:]/n[i])
        if data_dim == 2: Sxsd[i,:] = std[i,:]*std[i,:]*(n[i]-1) + Sx[i,:]*Sx[i,:]/n[i]
    txx = np.nansum(Sxsd,axis=0)
    tn = np.nansum(n)
    s = abs((txx-(tx*tx)/tn)/(tn-1))**0.5 # combined SD
    nans = np.isnan(data).all(axis=0)
    m[nans] = np.nan
    s[nans] = np.nan
    return(m,s) 

