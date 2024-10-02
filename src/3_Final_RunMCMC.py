

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Set up and run the MCMC

# Author: Eliza Harris
# Created on Tue Jun 16 13:17:17 2020

### 1_Final_Preamble.py must first be run to load all required data!
### 2_Final_Model_v4.py must first be run to load the model function

import glob
from src.utils import prob_u, prob_g

#%% Set up for the MCMC
# Run this section even if only looking at results with script 4!

run_name = configs.run_name

# Choose the step sizes (this is a fixed steplength - proportion of standard deviation)
step_lengths = [0.5,0.25,0.75]

# Define the observations to be used for param optimisation: obs = [ value, unc, group number ] (obs with the same group number are not varied independently in the MCMC)   
obs_conc = np.hstack((np.array(N2O_atmos[["N2O_ppb","sd_N2O_ppb"]]),np.zeros((N2O_atmos[["N2O_ppb"]].shape[0],1))+1))
obs_d15 = np.hstack((np.array(N2O_atmos[["d15Nbulk_permil","sd_d15Nbulk_permil"]]),np.zeros((N2O_atmos[["N2O_ppb"]].shape[0],1))+2))
obs_SP = np.hstack((np.array(N2O_atmos[["SP_permil","sd_SP_permil"]]),np.zeros((N2O_atmos[["N2O_ppb"]].shape[0],1))+3))
obs_flux = np.vstack((N2O_fluxes_zones[:,4],N2O_fluxes_zones[:,5],(np.zeros((N2O_fluxes_zones[:,4].shape[0],1))+4).transpose())).transpose()
obs = np.vstack((obs_conc,obs_d15,obs_SP,obs_flux))
obs_3last = np.cumsum(np.array((obs_conc.shape[0],obs_d15.shape[0],obs_SP.shape[0],obs_flux.shape[0])))-2 # row refs for the third-last obs for probabilities

# Define model uncertainty
mod_uncert = obs[:,1].copy()
mod_uncert_base = [0.5,0.1,0.1,0.5] # Uncertainty for "2020" data for N2O, d15N, SP, flux; gets larger with time in past
mod_uncert_scale = (abs(tstarts-np.nanmax(tstarts))/(np.nanmax(tstarts)-np.nanmin(tstarts))*2+1) # scale error by 1 (2020) to 3 (1840) to account for higher uncertainty in older data
for n in range(0,3):   
    mod_uncert[obs[:,2]==(n+1)] = mod_uncert_scale*mod_uncert_base[n] # Model uncertainty increases backward in time 
mod_uncert[obs[:,2]==4] = mod_uncert_base[3] # No time dimension on fluxes

#%% Set up a folder for the MCMC results
# Dont't run this section if only looking at results with script 4!

if not os.path.exists('outputs/MCMCresults/'+run_name): # Check if the directory of results already exists and make if needed
   os.makedirs('outputs/MCMCresults/'+run_name)
else: 
    success = False
    while success == False:
        print("A directory with this name already exists... ")
        choice = "x"
        while ((choice=="a")+(choice=="o")==0):
            choice = input('Enter choice: \n o = overwrite and delete old results \n a = add results to directory \n') 
        if choice=="o":
            choice2 = input('Overwrite: Are you sure? Y to confirm. \n') 
            if choice2!="Y": continue
            filenames = glob.glob('outputs/MCMCresults/'+run_name+'/*') # Empty the directory
            res = [ os.remove(f) for f in filenames ]
            print("Choice o: Directory emptied.")
            success = True
        if choice=="a":
            print("Results will be added to: outputs/MCMCresults/"+run_name)
            success = True

#%% Start the MCMC   
# Don't run this section if only looking at results with script 4!

runs = 0 # Start counting of runs
n_iterations = 500 # Define the number of model iterations per run (each run uses a single step size)
n_runs = 9 # Define the number of runs (subsequent runs move through the different step sizes)
if n_runs % len(step_lengths) != 0:
    print("Number of runs does not divide evenly by step size!!")  

for runs in range(4,n_runs):
    step_length = step_lengths[runs % len(step_lengths)] 
    
    # define the parameter list to be optimised: x = [ value, unc/lo, hi, dist ] where dist = 0 = uniform and dist = 1 = gaussian
    params = np.zeros((12,4))
    params[0,:] = (265, c_prea[1], c_prea[1], 1)#(c_prea[0], c_prea[1], c_prea[1], 1) # preanth N2O conc
    params[1,:] = (scale_fitN2[1],0.7,2,0) # N2 emission scaling factor
    params[2,:] = (fracex,0.3,1.001,0) # frac expression factor
    params[3,:] = (fertEFred,0,1.001,0) # fert EF red
    params[4,:] = (temp_sens,0.04,0.04,1)
    params[5,:] = (ltPD_ltPI[0],ltPD_ltPI[1],ltPD_ltPI[1],1)
    params[6,:] = (121,ltPD[1],ltPD[1],1) #(ltPD[0],ltPD[1],ltPD[1],1)
    params[7,:] = (d15_ocean[0],d15_ocean[1],d15_ocean[1],1)
    params[8,:] = (SP_ocean[0],SP_ocean[1],SP_ocean[1],1)
    params[9,:] = (d15_prea[0],d15_prea[1],d15_prea[1],1) 
    params[10,:] = (scale_fitNO[1],0.7,1.3,0) # NO emission scaling factor: restictions keep it in the range of obs
    params[11,:] = (SP_prea[0],SP_prea[1],SP_prea[1],1) 
    params_orig = params.copy()
    # For troubleshooting, some values to try...: [285,1.7,0.9,1.2,1.1,132,131,5.4,7.1,10]
    # params[:,0] = [280,1.7,0.4,0.3,1.1,1.06,121,5.4,7.1,10,1,SP_prea[0]]

    # select if the most recent model final priors should be used? Use N if major changes have been made    
    UseMostRecent = "Y" 
    if UseMostRecent == "Y" :
        filenames = glob.glob('outputs/MCMCresults/'+run_name+'/ModelResults_MCMC_*.nc')
        if len(filenames)>0: 
            order = np.zeros((len(filenames))) # make sure file names ordered by date
            for n in range(0,len(filenames)):
                tmp = filenames[n].replace('.', '_').split("_")
                order[n] = int(tmp[-3]+tmp[-2])
            r = np.where(order==np.max(order))[0][0]
            ncres = nc4.Dataset(filenames[r],'r')
            params[:,0] = ncres.variables["last_prior"][:]  
            if params[0,0] > 9e36: params[:,0] = ncres.variables["prior_res"][:][:,0]
        else:
            print("No previous results existing; using base priors")
    
    # Create space for solutions and add first priors
    prior_solutions = np.zeros((len(params[:,0]),n_iterations))+np.nan
    prior_solutions[:,0] = params[:,0]
    tmp = model(x=prior_solutions[:,0]) # First run of the model with prior inputs
    datarng = np.where(~np.isnan(d15N_grid) & ~np.isinf(d15N_grid))
    N_summary = tmp[5] # Get the full N summary results for comparison with obs
    Ntot_prev = np.array((N_summary[N_summary[:,0] == 1860,1],N_summary[N_summary[:,0] == 2010,1])) # As well as the total emissions in the PI and 2010
    modres_prev = expand_modres(tmp) # Put model results in same format as obs

    # Create space for the walked obs in each iteration also
    obs_solutions = np.zeros((len(obs[:,0]),n_iterations))+np.nan
    obs_solutions[:,0] = obs[:,0]
    
    # Create space for the other model results
    accept_MC = np.zeros((n_iterations,1)) # States whether the model in this iteration was accepted or not
    accept_MC[0] = 1
    atmos_opt = np.zeros((n_iterations,2)) # Space for the parameters optimised in the atm model (strat trop flux and ocean emissions)
    atmos_opt[0,:] = [tmp[1][0],tmp[2][0]]
    soil_results = np.zeros((n_iterations,tmp[3].shape[0],tmp[3].shape[1])) # Space for other model outputs that we want to save
    soil_results[0,:,:] = tmp[3]
    Nproc_results = np.zeros((n_iterations,tmp[4].shape[0],tmp[4].shape[1]))
    Nproc_results[0,:,:] = tmp[4]
    
    # run the MCMC
    i_prev = 0 
    for i in range(1,n_iterations):
        # Get the input parameters from the previous step, walk them according to the step length
        prior_i = prior_solutions[:,i_prev].copy()
        prior_i_stepscale = params[:,1].copy() # For Gaussian params, the walk distance is the stdev * the step length
        prior_i_stepscale[params[:,3]==0] = abs(params[params[:,3]==0,2]-params[params[:,3]==0,1])*0.25 # For uniform params, the walk distance is 25% of the range * the step length
        prior_i += step_length * prior_i_stepscale * np.random.uniform(-1,1,len(params[:,0])) # "Walk": vary all the parameters independently
        prior_solutions[:,i] = prior_i.copy() # save the results for this iteration
        
        # Get the obs values from the previous step, walk them according to the step length
        obs_i = obs_solutions[:,i_prev].copy()
        tmp = np.random.uniform(-1,1,int(max(obs[:,2])))
        obs_i += step_length * obs[:,1] * tmp[(obs[:,2]-1).astype(int)] # Vary different observation groups independently from each other; all are Gaussian so use stdev
        obs_solutions[:,i] = obs_i.copy()

        ## now test if inputs and obs pass the acceptance test for priors
        # first check uniform priors
        tmp = np.where(params[:,3] == 0)
        if np.random.uniform(0,1,1) > prob_u(prior_i[tmp],params_orig[tmp,1],params_orig[tmp,2]) / prob_u(prior_solutions[tmp,i_prev],params_orig[tmp,1],params_orig[tmp,2]) : 
            print("don't accept 1") # 
            continue # don't accept
        # then check gaussian priors
        tmp = np.where(params[:,3] == 1)[0]
        if np.random.uniform(0,1,1) > prob_g(prior_i[tmp],params_orig[tmp,0],params_orig[tmp,1]) / prob_g(prior_solutions[tmp,i_prev],params_orig[tmp,0],params_orig[tmp,1]) : 
            print("don't accept 2") # 
            continue # don't accept
        # now do it for the data (all normally dist.)  (use only one from each group, as all in each group vary the same)  
        if np.random.uniform(0,1,1) > prob_g(obs_i[obs_3last],obs[obs_3last,0],obs[obs_3last,1]) / prob_g(obs_solutions[obs_3last,i_prev],obs[obs_3last,0],obs[obs_3last,1]) : 
            print("don't accept 3") # 
            continue # don't accept

        # Now run the model and check the acceptance, with randomisation selected for d15N
        ### FROM HERE: Why is the model so overly sensitive to param 2, frac ex??? Look in the soil model...
        tmp = model(prior_i,d15Nrandomise=0.05)
        print(i)
        N_summary = tmp[5]
        # For troubleshooting to get a stable model: Check total flux in 1860 and 2010 approx right acc. Tian 2019: in TgN y-1: 1860 6.3pm1.1, 2010 10pm2.2... much higher in Tian 2020??
        Ntot = np.array((N_summary[N_summary[:,0] == 1860,1],N_summary[N_summary[:,0] == 2010,1]))
        #if interval_random() > prob_g(Ntot[0],6.3,1.1*1) / prob_g(Ntot_prev[0],6.3,1.1*1) : continue # don't accept
        #if interval_random() > prob_g(Ntot[1],10,2.2*1) / prob_g(Ntot_prev[1],10,2.2*1) : continue # don't accept
        Ntot_prev = Ntot
        # then check the rest of the observations
        modres = expand_modres(tmp)
        if isinstance(tmp[1], float): # For some reason, this is nearly always a float in an array, but sometimes just a float...
            atmos_opt[i,:] = (tmp[1],tmp[2])
        else:
            atmos_opt[i,:] = (tmp[1][0],tmp[2][0])
        soil_results[i,:,:] = tmp[3]
        Nproc_results[i,:,:] = tmp[4]
        tmp_c = abs((modres-obs_i)/mod_uncert) # current model: deviation from observations relative to the model uncertainty
        tmp_p = abs((modres_prev-obs_solutions[:,i_prev])/mod_uncert) # prev model: deviation from observations relative to the model uncertainty
        mod_summ_c = np.zeros(int(max(obs[:,2]))); 
        mod_summ_p = np.zeros(int(max(obs[:,2]))) # summarise model by observation type: create space
        for n in range(0,int(max(obs[:,2]))) : # find average deviation relative to model uncertainty for each observation type, not including EFs (too uncertain)
            mod_summ_c[n] = np.nanmean(tmp_c[ obs[:,2]==(n+1) ])
            mod_summ_p[n] = np.nanmean(tmp_p[ obs[:,2]==(n+1) ])
        current_mod_prob = prob_g(mod_summ_c,mod_summ_c*0,mod_summ_c*0+1) # gaussian likelihood that this model res are equal to 0 deviation (std 1 since we normalised)
        prev_mod_prob = prob_g(mod_summ_p,mod_summ_p*0,mod_summ_p*0+1)
        # note: this probability is not valid in absolute terms but gives relative probability between previous and current model correctly, and delivers results with a useful number of decimals for the metropolis hastings to work well
        if (prev_mod_prob==0): 
            prev_mod_prob=1e-300 # make sure we are not trying to divide by zero is the probability is below the number of decimals that can be represented in a float
        if np.random.uniform(0,1,1) > current_mod_prob/prev_mod_prob : 
            print("don't accept mod") # 
            continue # don't accept
        # it's passed all tests, append it to the solutions and make it the new starting point
        i_prev = i
        modres_prev = modres
        accept_MC[i] = 1
        go = sum(accept_MC==1)
        print("found a solution! n="+str(int(go))+" - rate="+str(int(100*go/(i+1) ))+"%" )
        print(modres)
        
    ### save MCMC results
    
    # take means of results to reduce data size  
    go = np.where(accept_MC==1)[0]
    soil_res_mean = np.nanmean(soil_results[go,:,:],axis=0)
    soil_res_sd = np.nanstd(soil_results[go,:,:],axis=0)
    Nproc_mean = np.nanmean(Nproc_results[go,:,:],axis=0)
    Nproc_sd = np.nanstd(Nproc_results[go,:,:],axis=0)
    prior_solutions_meansd = np.vstack((np.nanmean(prior_solutions[:,go],axis=1),np.nanstd(prior_solutions[:,go],axis=1))).transpose()
    atmos_opt_meansd = np.vstack((np.nanmean(atmos_opt[go,:],axis=0),np.nanstd(atmos_opt[go,:],axis=0))).transpose()
    obs_solutions_meansd = np.vstack((np.nanmean(obs_solutions[:,go],axis=1),np.nanstd(obs_solutions[:,go],axis=1))).transpose()
        
    name = ("MCMC_"+datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    ncout = nc4.Dataset('outputs/MCMCresults/'+run_name+"/ModelResults_"+name+".nc",'w','NETCDF4'); # using netCDF3 for output format 
    ncout.createDimension('points',1)
    ncout.createDimension('rows',2)
    ncout.createDimension('soilparams',soil_results.shape[2])
    ncout.createDimension('Nprocparams',Nproc_results.shape[2])
    ncout.createDimension('datarng',soil_results.shape[1])
    ncout.createDimension('n_params',12)
    ncout.createDimension('n_obs',obs_solutions_meansd.shape[0])
    n_its = ncout.createVariable('n_its','f4',('points'))
    n_its[:] = i+1
    n_accepted = ncout.createVariable('n_accepted','f4',('points'))
    n_accepted[:] = len(go)
    soil_res_m = ncout.createVariable('soil_res_m','f4',('datarng','soilparams'))
    soil_res_m.setncattr('cols','D1 = datarng (non-na indices of d15N in global grid), D2 = k_L, k_NH3, k_G, soil_NO3_steady, d15N, SP')
    soil_res_m[:,:] = soil_res_mean[:,:]
    soil_res_s = ncout.createVariable('soil_res_s','f4',('datarng','soilparams'))
    soil_res_s[:,:] = soil_res_sd[:,:]
    Nproc_res_m = ncout.createVariable('Nproc_res_m','f4',('datarng','Nprocparams'))
    Nproc_res_m.setncattr('cols','D1 = datarng (non-na indices of d15N in global grid), D2 = N2O_all, N2_all, NO_all, denit_N2O, nit_N2O, red_N2O')
    Nproc_res_m[:,:] = Nproc_mean[:,:]
    Nproc_res_s = ncout.createVariable('Nproc_res_s','f4',('datarng','Nprocparams'))
    Nproc_res_s[:,:] = Nproc_sd[:,:]
    prior_res = ncout.createVariable('prior_res','f4',('n_params',"rows"))
    prior_res[:,:] = prior_solutions_meansd[:,:]
    atmos_opt_res = ncout.createVariable('atmos_opt_res','f4',('rows',"rows"))
    atmos_opt_res[:,:] = atmos_opt_meansd[:,:]
    obs_res = ncout.createVariable('obs_res','f4',('n_obs',"rows"))
    obs_res[:,:] = obs_solutions_meansd[:,:]
    last_prior = ncout.createVariable('last_prior','f4',('n_params'))
    last_prior[:] = prior_solutions[:,go[-1]]
    print("last prior is "+str(last_prior[:]))
    # full prior res as text also to check coverage of parameter space
    np.savetxt('outputs/MCMCresults/'+run_name+"/AcceptMC_"+name+".txt", np.vstack((accept_MC.astype("int"),accept_MC.astype("int")*0+step_length)).transpose())
    np.savetxt('outputs/MCMCresults/'+run_name+"/PriorSolutions_"+name+".txt", prior_solutions)
    np.savetxt('outputs/MCMCresults/'+run_name+"/AtmosSolutions_"+name+".txt", atmos_opt)
    
#%% model check  
tmp = model(x= a,d15Nrandomise=0.05) # check with orig params
tmp = model(x= params[:,0],d15Nrandomise=0.05) # check with orig params
tmp = model(x= prior_solutions[:,go[-1]],d15Nrandomise=0.05) # check with last prior
res = tmp[0]
res_conc = res[0:int(len(res)/3)]
res_d15 = res[int(len(res)/3):int(2*len(res)/3)]
res_SP = res[int(2*len(res)/3):int(len(res))]
fig, ax = plt.subplots(3,1)
ax[0].plot(tstarts,obs_conc[:,0],"bo")
ax[0].plot(tstarts,res_conc,"cx")
ax[0].legend((["obs","mod"]))
ax[0].set_ylabel("N2O trop (ppb)")
ax[1].plot(tstarts,obs_d15[:,0],"bo")
ax[1].plot(tstarts,res_d15,"cx")
ax[1].set_ylabel("d15N (permil)")
ax[2].plot(tstarts,obs_SP[:,0],"bo")
ax[2].plot(tstarts,res_SP,"cx")
ax[2].set_ylabel("SP (permil)")
plt.show()

