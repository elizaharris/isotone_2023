

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Set up and run the MCMC

# Author: Eliza Harris
# Created on Tue Jun 16 13:17:17 2020

import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as math
import time
import datetime
from importlib import reload
import netCDF4 as nc4

parentdir = os.getcwd()
os.chdir(parentdir+"/isotone_2025") # Change to the isotone repo as working directory

# Import the full preamble (Pt1_Final_Preamble.py, to load all required data!)
import src.Pt1_Final_Preamble as preamble
from src.Pt1_Final_Preamble import configs # The stuff used often also imported separately
import src.Pt2_Final_Model_v4 as model_script
from src.Pt2_Final_Model_v4 import model
import src.utils as utils
from src.utils import prob_u, prob_g, export_model_results

#%% Set up for the MCMC
# Run this section even if only looking at results with script 4!

run_name = configs.run_name

# define the parameter list : x = [ value, unc/lo, hi, dist ] where dist = 0 = uniform and dist = 1 = gaussian
params = np.zeros((12,4))
params[0,:] = (265, preamble.c_prea[1], preamble.c_prea[1], 1)#(c_prea[0], c_prea[1], c_prea[1], 1) # preanth N2O conc
params[1,:] = (preamble.scale_fitN2[1],0.7,2,0) # N2 emission scaling factor
params[2,:] = (preamble.fracex,0.3,1.001,0) # frac expression factor
params[3,:] = (preamble.fertEFred,0,1.001,0) # fert EF red
params[4,:] = (preamble.temp_sens,0.04,0.04,1)
params[5,:] = (preamble.ltPD_ltPI[0],preamble.ltPD_ltPI[1],preamble.ltPD_ltPI[1],1)
params[6,:] = (121,preamble.ltPD[1],preamble.ltPD[1],1) #(ltPD[0],ltPD[1],ltPD[1],1)
params[7,:] = (preamble.d15_ocean[0],preamble.d15_ocean[1],preamble.d15_ocean[1],1)
params[8,:] = (preamble.SP_ocean[0],preamble.SP_ocean[1],preamble.SP_ocean[1],1)
params[9,:] = (preamble.d15_prea[0],preamble.d15_prea[1],preamble.d15_prea[1],1) 
params[10,:] = (preamble.scale_fitNO[1],0.7,1.3,0) # NO emission scaling factor: restictions keep it in the range of obs
params[11,:] = (preamble.SP_prea[0],preamble.SP_prea[1],preamble.SP_prea[1],1) 

# define the parameter list using optimised values from original IsoTONE study
params_post = np.zeros((12,4))
params_post[0,:] = (276, 2.2, 2.2, 1)#(c_prea[0], c_prea[1], c_prea[1], 1) # preanth N2O conc
params_post[1,:] = (0.9, 0.7, 1.1, 0) # N2 emission scaling factor
params_post[2,:] = (0.55, 0.5, 0.6, 0) # frac expression factor
params_post[3,:] = (0.3, 0.23, 0.37, 0) # fert EF red
params_post[4,:] = (1.1,0.03,0.03,1) # temp sensitivity
params_post[5,:] = (preamble.ltPD_ltPI[0],preamble.ltPD_ltPI[1],preamble.ltPD_ltPI[1],1) # lifetime PD to PI; stays same after MCMC
params_post[6,:] = (131, 6, 6, 1) #(ltPD[0],ltPD[1],ltPD[1],1)
params_post[7,:] = (5.3, 1.6, 1.6, 1) # d15N ocean
params_post[8,:] = (14.2, 9.6, 9.6, 1) # SP ocean
params_post[9,:] = (11.2, 0.3, 0.3, 1) # d15N prea
params_post[10,:] = (1.2, 1.1, 1.3, 0) # NO emission scaling factor: restictions keep it in the range of obs
params_post[11,:] = (19.8, 0.3, 0.3, 1) # SP prea

# Define the observations to be used for param optimisation: obs = [ value, unc, group number ] (obs with the same group number are not varied independently in the MCMC)   
obs_conc = np.hstack((np.array(preamble.N2O_atmos[["N2O_ppb","sd_N2O_ppb"]]),np.zeros((preamble.N2O_atmos[["N2O_ppb"]].shape[0],1))+1))
obs_d15 = np.hstack((np.array(preamble.N2O_atmos[["d15Nbulk_permil","sd_d15Nbulk_permil"]]),np.zeros((preamble.N2O_atmos[["N2O_ppb"]].shape[0],1))+2))
obs_SP = np.hstack((np.array(preamble.N2O_atmos[["SP_permil","sd_SP_permil"]]),np.zeros((preamble.N2O_atmos[["N2O_ppb"]].shape[0],1))+3))
obs = np.vstack((obs_conc,obs_d15,obs_SP))
obs_3last = np.cumsum(np.array((obs_conc.shape[0],obs_d15.shape[0],obs_SP.shape[0],)))-2 # row refs for the third-last obs for probabilities

datarng = np.where(~np.isnan(preamble.d15N_grid) & ~np.isinf(preamble.d15N_grid))

#%% model check  
import time
start = time.time()
#tmp = model(x= params[:,0],d15Nrandomise=0.05) # check with orig params
tmp = model(x= params_post[:,0], preamble=preamble, d15Nrandomise=0.05) # check with orig params
end = time.time()
print("Model run time:",(end-start), "s")

# Rough plots
res = tmp[0]
res_conc = res[0:int(len(res)/3)]
res_d15 = res[int(len(res)/3):int(2*len(res)/3)]
res_SP = res[int(2*len(res)/3):int(len(res))]
fig, ax = plt.subplots(3,1)
ax[0].plot(preamble.tstarts,obs_conc[:,0],"bo")
ax[0].plot(preamble.tstarts,res_conc,"cx")
ax[0].legend((["obs","mod"]))
ax[0].set_ylabel("N2O trop (ppb)")
ax[1].plot(preamble.tstarts,obs_d15[:,0],"bo")
ax[1].plot(preamble.tstarts,res_d15,"cx")
ax[1].set_ylabel("d15N (permil)")
ax[2].plot(preamble.tstarts,obs_SP[:,0],"bo")
ax[2].plot(preamble.tstarts,res_SP,"cx")
ax[2].set_ylabel("SP (permil)")
plt.show()

k_G = np.zeros(preamble.d15N_grid.shape)+np.nan; k_G[datarng] = tmp[3][:,2]
f_N2O = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2O[datarng] = tmp[4][:,0]
utils.plot_map(preamble.LON,preamble.LAT,k_G*f_N2O*100,"fraction of Nr lost to gas production",show=1)
utils.plot_map(preamble.LON,preamble.LAT,k_G*f_N2O*100,"fraction of Nr lost to N2O",show=1)
