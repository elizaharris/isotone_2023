#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Bring atmospheric datasets to the same scale and smooth for model input

"""
Created on Tue Mar  3 10:39:11 2020

@author: elizaharris
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import configs as configs

# Function to turn string date to approx. decimal
def date_as_decimal(date,sep) :
    [y,m,d] = str(date).split(sep)
    if int(y)<1000: y = int(y)+2000
    return int(y)+(int(m)-1)/12+(int(d)-1)/30/12

# Get data from Longfei Yu (Empa), from Cape Grim, JFJ and EGRIP
N2O_longfei = pd.read_csv("data/atmos_data/N2Otrends_longfei_sort_secondEGRIP_04082021.csv", sep=',')
date_dec = np.zeros((N2O_longfei.shape[0],1))
for n in range(0,N2O_longfei.shape[0]) :
    date = N2O_longfei["Date"][n]
    date_dec[n] = date_as_decimal(date,"-")
N2O_longfei["dt"] = date_dec
N2O_longfei["dataset"] = "Longfei"

# GAW data (from JFJ, provided by Martin Steinbacher (Empa))
N2O_GAW = pd.read_csv("data/atmos_data/JFJ_GAWdata_daily_v2.csv", sep=',')
date_dec = np.zeros((N2O_GAW.shape[0],1))
for n in range(0,N2O_GAW.shape[0]) :
    date = N2O_GAW["Date"][n]
    date_dec[n] = date_as_decimal(date,"-")
N2O_GAW["dt"] = date_dec
N2O_GAW = N2O_GAW.loc[N2O_GAW["N2O_ppb"]>-100] # remove nans
N2O_GAW["dataset"] = "GAW"

# Cape Grim conc data (CSIRO)
N2O_CG1 = pd.read_csv("data/atmos_data/CG_old_data.csv", sep=',')
date_dec = N2O_CG1["year"]+(N2O_CG1["month"]-1)/12+(N2O_CG1["day"]-1)/31/12
N2O_CG1["dt"] = date_dec
N2O_CG2 = pd.read_csv("data/atmos_data/CG_mid_data_short.csv", sep=',')
date_dec = N2O_CG2["year"]+(N2O_CG2["month"]-1)/12+(N2O_CG2["day"]-1)/31/12
N2O_CG2["dt"] = date_dec
N2O_CG3 = pd.read_csv("data/atmos_data/CG_new_data_short.csv", sep=',')
date_dec = N2O_CG3["year"]+(N2O_CG3["month"]-1)/12+(N2O_CG3["day"]-1)/31/12
N2O_CG3["dt"] = date_dec
# correct old- and mid-data to new data (different instruments, protocols...)
def correct_scale(scale_data_dt,scale_data,corr_data_dt,corr_data) :
    t_min = max(np.nanmin(scale_data_dt),np.nanmin(corr_data_dt)) # find the overlapping time window
    t_max = min(np.nanmax(scale_data_dt),np.nanmax(corr_data_dt))
    t_step = (t_max-t_min)/50
    t_window = np.arange(t_min,t_max+t_step,step = t_step)
    scale_data_i = np.interp(t_window,scale_data_dt,scale_data) # interp both to steps within the overlap window
    corr_data_i = np.interp(t_window,corr_data_dt,corr_data)
    offset = np.nanmean(corr_data_i - scale_data_i)
    corr_data_out = corr_data - offset
    return(corr_data_out,offset)
N2O_CG2_raw = N2O_CG2.copy()
(N2O_CG2["value"],offset) = correct_scale(scale_data_dt = N2O_CG3["dt"],scale_data = N2O_CG3["value"],corr_data_dt = N2O_CG2_raw["dt"],corr_data = N2O_CG2_raw["value"])
N2O_CG1_raw = N2O_CG1.copy()
(N2O_CG1["value"],offset) = correct_scale(scale_data_dt = N2O_CG2["dt"],scale_data = N2O_CG2["value"],corr_data_dt = N2O_CG1_raw["dt"],corr_data = N2O_CG1_raw["value"])
# combine the datasets
N2O_CG = pd.concat([N2O_CG1,N2O_CG2,N2O_CG3],join='inner',sort=False)
N2O_CG["N2O_ppb"] = N2O_CG["value"]
N2O_CG = N2O_CG.loc[N2O_CG["N2O_ppb"]>-100] # remove nans
N2O_CG["sd_N2O_ppb"] = N2O_CG["value_unc"]
N2O_CG["dataset"] = "CG_online"
N2O_CG["Site"] = "CapeGrim"

# Machida data (conc only)
N2O_Mach = pd.read_csv("data/atmos_data/machida1995.csv", sep=',')
N2O_Mach["dataset"] = "Machida"

# Park et al. data: from 10.1038/NGEO1421 (Trends and seasonal cycles in the isotopic composition of nitrous oxide since 1940)
N2O_Park1 = pd.read_csv("data/atmos_data/Park2011_FirnData_edit.csv", sep=',')
N2O_Park1["dataset"] = "Park_Firn"
N2O_Park2 = pd.read_csv("data/atmos_data/Park2011_ArchiveData_edit.csv", sep=',')
date_dec = np.zeros((N2O_Park2.shape[0],1))
for n in range(0,N2O_Park2.shape[0]) :
    date = N2O_Park2["Date"][n]
    [y,m,d] = [str(date)[0:4],str(date)[4:6],str(date)[6:8]]    
    date_dec[n] = int(y)+(int(m)-1)/12+(int(d)-1)/31
N2O_Park2["dt"] = date_dec
N2O_Park2["dataset"] = "Park_Archive"
N2O_Park2 = N2O_Park2.loc[N2O_Park2["SP_permil"]<1000] # remove nans
N2O_Park = pd.concat([N2O_Park1,N2O_Park2],join='outer',sort=False)
# rescale d15N to Longfei's data to account for isotopic scale offsets
N2O_Park1_raw = N2O_Park1.copy()
(N2O_Park1["d15Nbulk_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["d15Nbulk_permil"],corr_data_dt = N2O_Park1_raw["dt"],corr_data = N2O_Park1_raw["d15Nbulk_permil"])
(N2O_Park1["SP_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["SP_permil"],corr_data_dt = N2O_Park1_raw["dt"],corr_data = N2O_Park1_raw["SP_permil"])
N2O_Park2_raw = N2O_Park2.copy()
(N2O_Park2["d15Nbulk_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["d15Nbulk_permil"],corr_data_dt = N2O_Park2_raw["dt"],corr_data = N2O_Park2_raw["d15Nbulk_permil"])
(N2O_Park2["SP_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["SP_permil"],corr_data_dt = N2O_Park2_raw["dt"],corr_data = N2O_Park2_raw["SP_permil"])

# Sowers data
N2O_Sowers = pd.read_csv("data/atmos_data/sowersdata.csv", sep=',')
N2O_Sowers["dataset"] = "Sowers_2002"
# rescale d15N to Longfei's data
N2O_Sowers_raw = N2O_Sowers.copy()
# (N2O_Sowers["d15Nbulk_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["d15Nbulk_permil"],corr_data_dt = N2O_Sowers_raw["dt"],corr_data = N2O_Sowers_raw["d15Nbulk_permil"])
# There is no significant difference in the overlap period

# Ice core data (Law Dome), CSIRO
N2O_LawDome = pd.read_csv("data/atmos_data/Law_Dome_Ice_Core_2000-Year_CO2__CH4__N2O_and_d13C-CO2/Law_Dome_N2O.csv", sep=',')
N2O_LawDome["dataset"] = "LawDome"
N2O_LawDome["dt"] = N2O_LawDome["N2O Age (year AD)"]
N2O_LawDome = N2O_LawDome.iloc[N2O_LawDome["N2O Age (year AD)"].argsort()]
# rescale d15N to Longfei's data
N2O_LawDome_raw = N2O_LawDome.copy()
(N2O_LawDome["N2O_ppb"],offset) = correct_scale(scale_data_dt = N2O_CG["year"],scale_data = N2O_CG["N2O_ppb"],corr_data_dt = N2O_LawDome["dt"],corr_data = N2O_LawDome["N2O (ppb)"])

# Data from Ghosh et al. 2023 JGR-A (Supp Info)
### Note: This is a compilation, not all measured in this paper! See paper Fig 1 for which ice core was from which paper.
N2O_Ghosh = pd.read_csv("data/atmos_data/Ghosh_2023.csv", sep=',')
N2O_Ghosh["dataset"] = "Ghosh"
N2O_Ghosh["dt"] = N2O_Ghosh["Age"]
N2O_Ghosh = N2O_Ghosh[N2O_Ghosh["dt"].isna()==False]
N2O_Ghosh = N2O_Ghosh.iloc[N2O_Ghosh["dt"].argsort()]
# rescale d15N to Longfei's data
N2O_Ghosh_raw = N2O_Ghosh.copy()
(N2O_Ghosh["N2O_ppb"],offset) = correct_scale(scale_data_dt = N2O_CG["year"],scale_data = N2O_CG["N2O_ppb"],corr_data_dt = N2O_Ghosh_raw["dt"],corr_data = N2O_Ghosh_raw["N2O_ppb"])
(N2O_Ghosh["d15Nbulk_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["d15Nbulk_permil"],corr_data_dt = N2O_Ghosh_raw["dt"],corr_data = N2O_Ghosh_raw["d15Nbulk_permil"])
(N2O_Ghosh["SP_permil"],offset) = correct_scale(scale_data_dt = N2O_longfei["dt"],scale_data = N2O_longfei["SP_permil"],corr_data_dt = N2O_Ghosh_raw["dt"],corr_data = N2O_Ghosh_raw["SP_permil"])

# join the datasets
N2O_atmos = pd.concat([N2O_Park1,N2O_Park2,N2O_GAW,N2O_longfei,N2O_Sowers,N2O_Mach,N2O_CG,N2O_LawDome,N2O_Ghosh],join='outer',sort=False)

# Function to find a weighted average
# Note: higher weights = more contribution!
def weighted_avg_and_se(values, weights):
    v = (~np.isnan(values))
    if sum(v)>1:
        values = values[v]
        weights = weights[v]
        w = np.isnan(weights)
        if sum(~w)==0: weights[w] = 2 # Replace nan weights
        else: weights[w] = np.max(weights)
        average = np.average(values, weights=weights)
        # Fast and numerically precise:
        variance = np.average((values-average)**2, weights=weights)
        se = math.sqrt(abs(variance))/((sum(v)-1)**0.5)
        sd = math.sqrt(abs(variance))
        return (average, sd, sum(v))
    if sum(v)==1:
        average = values[v].values[0]
        se = weights[v].values[0]
        sd = weights[v].values[0]
        return (average, sd, 1)

# Create the start years for each time window
n_windows = len(configs.mod_obs_time["ystart"])
for n in np.arange(0,n_windows):
    tmp = np.arange(configs.mod_obs_time["ystart"][n],configs.mod_obs_time["yend"][n],configs.mod_obs_time["yres"][n])
    if n == 0:
        tstarts = tmp
    else:
        tstarts = np.hstack((tstarts,tmp))

# Weighted averages within time windows to smooth
N2O_atmos_s = np.zeros((len(tstarts),13))+np.nan
N2O_atmos_s[:,0] = tstarts
for n in range(0,(len(tstarts)-1)): 
    tmp = ((N2O_atmos["dt"]>=tstarts[n]) & (N2O_atmos["dt"]<tstarts[n+1]))
    N2O_atmos_s[n,(1,5,9)] = weighted_avg_and_se(N2O_atmos["N2O_ppb"][tmp],1/N2O_atmos["sd_N2O_ppb"][tmp])
    N2O_atmos_s[n,(2,6,10)] = weighted_avg_and_se(N2O_atmos["d15Nbulk_permil"][tmp],1/N2O_atmos["sd_d15Nbulk_permil"][tmp])
    N2O_atmos_s[n,(4,8,12)] = weighted_avg_and_se(N2O_atmos["d18O_permil"][tmp],1/N2O_atmos["sd_d18O_permil"][tmp])
    N2O_atmos_s[n,(3,7,11)] = weighted_avg_and_se(N2O_atmos["SP_permil"][tmp],1/N2O_atmos["sd_SP_permil"][tmp])
N2O_atmos_df = pd.DataFrame(N2O_atmos_s)
N2O_atmos_df.columns = ['dt','N2O_ppb','d15Nbulk_permil','SP_permil','d18O_permil',
                       'sd_N2O_ppb','sd_d15Nbulk_permil', 'sd_SP_permil','sd_d18O_permil',
                       'n_N2O_ppb','n_d15Nbulk_permil', 'n_SP_permil','n_d18O_permil']
# Fill sd columns with max sd
for c in ['sd_N2O_ppb','sd_d15Nbulk_permil', 'sd_SP_permil','sd_d18O_permil']:
    r = N2O_atmos_df[c].isna()
    N2O_atmos_df[c][r] = np.nanmax(N2O_atmos_df[c])

# plot the data
datasets = ("Longfei","Park_Archive","Park_Firn","GAW","Machida","Sowers_2002","CG_online","LawDome","Ghosh")
syms = ("o","x","o","o","d","x","o","x","x")
fig, ax = plt.subplots(3,1)
# Plot each of the datasets
for n in range(0,len(datasets)):
    a = ax[0].plot(N2O_atmos["dt"][N2O_atmos["dataset"]==datasets[n]],N2O_atmos["N2O_ppb"][N2O_atmos["dataset"]==datasets[n]],syms[n])
    b = ax[1].plot(N2O_atmos["dt"][N2O_atmos["dataset"]==datasets[n]],N2O_atmos["d15Nbulk_permil"][N2O_atmos["dataset"]==datasets[n]],syms[n])
    c = ax[2].plot(N2O_atmos["dt"][N2O_atmos["dataset"]==datasets[n]],N2O_atmos["SP_permil"][N2O_atmos["dataset"]==datasets[n]],syms[n])
# Plot the average for each of the isotopes
for n,p in enumerate(['N2O_ppb','d15Nbulk_permil','SP_permil']):
    ax[n].plot(N2O_atmos_df["dt"],N2O_atmos_df[p],"o",c="orange")
    #ax[n].plot(N2O_atmos_df["dt"],N2O_atmos_df[p]-N2O_atmos_df["sd_"+p],":",c="orange")
    #ax[n].plot(N2O_atmos_df["dt"],N2O_atmos_df[p]+N2O_atmos_df["sd_"+p],":",c="orange")
ax[0].set_ylabel("N2O trop (ppb)")
ax[1].set_ylabel("d15N (permil)")
ax[2].set_ylabel("SP (permil)")
for n in np.arange(3): 
    ax[n].set_xlim(1730,2030)
ax[0].legend(datasets,fontsize=5)
fig.tight_layout()
filename = "figs/atmosdata_overview"
plt.savefig(filename+".png") 
plt.savefig(filename+".pdf") 
fig.show() 

# trends for the past X years
from sklearn.linear_model import LinearRegression
trendyears = (1980,2020)
column = 'SP_permil'
a = np.where((N2O_atmos_df["dt"]>=trendyears[0]) & (N2O_atmos_df["dt"]<=trendyears[1]) & (~np.isnan(N2O_atmos_df[column])))
X = np.array(N2O_atmos_df["dt"].iloc[a]).reshape(-1, 1)
y = N2O_atmos_df[column].iloc[a]
model = LinearRegression().fit(X, y)
print(model.coef_)
# Same result but with p values...
import statsmodels.api as sma
X2  = sma.add_constant(X)
est = sma.OLS(y, X2)
est2 = est.fit()
print(est2.summary())
# And for d15N bulk
column = 'd15Nbulk_permil'
a = np.where((N2O_atmos_df["dt"]>=trendyears[0]) & (N2O_atmos_df["dt"]<=trendyears[1]) & (~np.isnan(N2O_atmos_df[column])))
X = np.array(N2O_atmos_df["dt"].iloc[a]).reshape(-1, 1)
y = N2O_atmos_df[column].iloc[a]
X2  = sma.add_constant(X)
est = sma.OLS(y, X2)
est2 = est.fit()
print(est2.summary())

# write to a summary file
N2O_final = N2O_atmos[["dt","dataset",'N2O_ppb','d15Nbulk_permil','d18O_permil','SP_permil',
                       'sd_N2O_ppb','sd_d15Nbulk_permil','sd_d18O_permil', 'sd_SP_permil','CO_ppb', 'CH4_ppb', 'CO2_ppm']]
N2O_final.to_csv('data/atmos_data/AllAtmosData_wGhosh.csv', index = False)
N2O_atmos_df.to_csv('data/atmos_data/SummAtmosData_wGhosh.csv', index = False)
