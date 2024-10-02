
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Import and plot detailed MCMC results (parameter space, stability, etc.)

# Author: Eliza Harris
# Created on Tue Apr  7 10:19:28 2020

### 1_Final_Preamble.py must first be run to load all required data!
### 2_Final_Model_v4.py must first be run to load the model function
### First section of 3_Final_RunMCMC.py must first be run to load the MCMC params

import glob
import scipy.stats as stats
import seaborn as sns
from src.utils import combined_mean, fix_acceptMC
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#%% Import MCMC results 

# MC results are in the netcdf:
# accept_MC_res = 0 for not accept, 1 for accept
# obs_solutions: observations posterior: D1 = N2O, d15N, SP end on end, D2 = n_its
# prior_solutions: D1 = the params, D2 = n_its
# atmos_opt: D1 = n_its, D2 = T_Sm_new, F_ocean_new (optimised within atm model)
# soil_results: D1 = n_its, D2 = datarng (non-na indices of d15N in global grid), D3 = k_L, k_NH3, k_G, soil_NO3_steady, d15N, SP )))
# Nproc_results: D1 = n_its, D2 = datarng (non-na indices of d15N in global grid), D3 = N2O_all, N2_all, NO_all, denit_N2O, nit_N2O, red_N2O

run_name = configs.run_name

# find all files (by datetime, same in all)
dates = glob.glob('outputs/MCMCresults/'+run_name+'/AcceptMC_MCMC_*')
dates_as_num = np.zeros(len(dates))
for n in range(0,len(dates)) :
    tmp = dates[n].split("_")
    dates[n] = (tmp[3]+"_"+tmp[4].split(".")[0])
    dates_as_num[n] = int(tmp[3]+tmp[4].split(".")[0])
# make sure dates are in order
tmp = dates_as_num.argsort()
dates_o = [ dates[i] for i in tmp ]

# Sometimes the accept_MC file has accept-stepsize as a single column, other times as two columns; correct this...
def fix_acceptMC(accept_MC_res): 
    if len(accept_MC_res.shape) == 1: 
        new = np.zeros((int(accept_MC_res.shape[0]/2),2))
        new[:,0] = accept_MC_res[0:int(accept_MC_res.shape[0]/2)]
        new[:,1] = accept_MC_res[int(accept_MC_res.shape[0]/2):]
    else : new = accept_MC_res
    return(new)

# Get the full prior solutions 
for n in range(0,len(dates_o)) :
    name = ("MCMC_"+dates_o[n])
    if n==0: # read files in
        accept_MC_res = fix_acceptMC( np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/AcceptMC_"+name+".txt") )
        stepsize = accept_MC_res[:,1].copy()
        accept_MC_res = accept_MC_res[:,0]
        prior_solutions = np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/PriorSolutions_"+name+".txt")
        atmos_opt = np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/AtmosSolutions_"+name+".txt")
    else : # cocatenate
        tmp = fix_acceptMC( np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/AcceptMC_"+name+".txt") )
        accept_MC_res = np.append(accept_MC_res,tmp[:,0])
        stepsize = np.append(stepsize,tmp[:,1])
        tmp = np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/PriorSolutions_"+name+".txt")
        prior_solutions = np.hstack((prior_solutions,tmp))
        tmp = np.loadtxt(fname = "outputs/MCMCresults/"+run_name+"/AtmosSolutions_"+name+".txt")
        atmos_opt = np.concatenate((atmos_opt,tmp),0)

### then get and combine the obs data
name = ("MCMC_"+dates_o[0])
ncres = nc4.Dataset("outputs/MCMCresults/"+run_name+"/ModelResults_"+name+".nc",'r') # Get the first nc in order to set up space for all the results
datarng = np.where(~np.isnan(d15N_grid) & ~np.isinf(d15N_grid))
n_acc = np.zeros(len(dates)) # number accepted
n_its = np.zeros(len(dates)) # number of iterations
obs_res_msd = np.zeros((len(dates),ncres.dimensions["n_obs"].size,2))
soilmodres_m = np.zeros((len(dates),ncres.dimensions["datarng"].size,ncres.dimensions["soilparams"].size))
Nprocesses_m = np.zeros((len(dates),ncres.dimensions["datarng"].size,ncres.dimensions["Nprocparams"].size))
soilmodres_s = np.zeros((len(dates),ncres.dimensions["datarng"].size,ncres.dimensions["soilparams"].size))
Nprocesses_s = np.zeros((len(dates),ncres.dimensions["datarng"].size,ncres.dimensions["Nprocparams"].size))
for n in range(0,len(dates)) :
    name = ("MCMC_"+dates_o[n])
    ncres = nc4.Dataset("outputs/MCMCresults/"+run_name+"/ModelResults_"+name+".nc",'r')
    n_acc[n] = ncres.variables["n_accepted"][:]
    n_its[n] = ncres.variables["n_its"][:]
    obs_res_msd[n,:,:] = ncres.variables["obs_res"][:,:]
    soilmodres_m[n,:,:] = ncres.variables["soil_res_m"][:,:]
    Nprocesses_m[n,:,:] = ncres.variables["Nproc_res_m"][:,:]
    soilmodres_s[n,:,:] = ncres.variables["soil_res_s"][:,:]
    Nprocesses_s[n,:,:] = ncres.variables["Nproc_res_s"][:,:]

# Find the mean of the posteriors, accounting for the different number of accepted values in each run
post_obs, post_obs_sd = combined_mean(data = obs_res_msd[:,:,0], std = obs_res_msd[:,:,1], n = n_acc) # Each slice is the mean/std from each "run" (dim = nruns * nobs) but the combined mean/std needs to account for differing numbers of accepted values in each run
post_kG_m, post_kG_sd = combined_mean(data = soilmodres_m[:,:,2], std = soilmodres_s[:,:,2], n = n_acc)
post_fN2O_m, post_fN2O_sd = combined_mean(data = Nprocesses_m[:,:,0], std = Nprocesses_s[:,:,0], n = n_acc)

#%% MCMC diagnostics

# all tested values
full_params =  pd.DataFrame(np.hstack((prior_solutions.transpose(),atmos_opt)))
full_params.columns = ["c_prea","N2 fit scale","frac express","fert EF_red","temp sens","lifetime PD/PI","lifetime PD","d15N ocean","SP ocean","d15N prea","NO fit scale","SP prea","TtoS mass", "F ocean"]
full_params[full_params.columns[-2]] = full_params[full_params.columns[-2]]/1000*MW_air
# accepted posterior values
go = np.where(accept_MC_res==1)[0]
post = np.nanmean(full_params.iloc[go],axis=0)
post_sd = np.nanstd(full_params.iloc[go],axis=0)
# comllect the prior values and ranges
full_prior = np.zeros((full_params.shape[1],4)) 
full_prior[0:params.shape[0],:] = params_orig.copy()  
full_prior[-2,:] = [T_S[0],T_S[1],T_S[1],1]
full_prior[-1,:] = [F_ocean[0],F_ocean[1],F_ocean[1],1]
# select to plot only every nth point of the full data (figure too large otherwise)
subset = np.arange(0,full_params.shape[0],4)

### check the coverage of the parameter space 
fig, ax = plt.subplots(2,3,figsize=(8,6))
i1 = 0
i2 = 0
col_pairs = [[0,2],[3,4],[5,6],[7,9],[8,11],[1,10]] # Chose the params to plot against each other
for n in range(0,6) :
    cols = full_params.columns[col_pairs[n]] # Columns as names
    cols_n = col_pairs[n] # Columns as numbers
    # plot all points coloured by stepsize
    cp = ax[i1,i2].scatter(full_params[cols[0]].iloc[subset],full_params[cols[1]].iloc[subset],c=stepsize[subset],marker=".",s=6,alpha=0.6) 
    # plot accepted points in grey
    ax[i1,i2].scatter(full_params[cols[0]].iloc[go],full_params[cols[1]].iloc[go],c="silver",marker=".",s=3) 
    # Plot lines to show range of y-axis variable : Post in red, prior in blue
    ax[i1,i2].plot([post[cols_n[0]]]*3,[ post[cols_n[1]]-post_sd[cols_n[1]],post[cols_n[1]],post[cols_n[1]]+post_sd[cols_n[1]] ],"r")
    if full_prior[cols_n[1],3] == 1: # gaussian error: solid lines
        ax[i1,i2].plot([full_prior[cols_n[0],0]]*3,[ full_prior[cols_n[1],0]-full_prior[cols_n[1],1],full_prior[cols_n[1],0],full_prior[cols_n[1],0]+full_prior[cols_n[1],1] ],"b")
    else: 
        ax[i1,i2].plot([full_prior[cols_n[0],0]]*3,[ full_prior[cols_n[1],1],full_prior[cols_n[1],0],full_prior[cols_n[1],2] ],"b--") # uniform: dashed
    # Plot lines to show range of x-axis variable : Post in red, prior in blue
    ax[i1,i2].plot([ post[cols_n[0]]-post_sd[cols_n[0]],post[cols_n[0]],post[cols_n[0]]+post_sd[cols_n[0]] ],[post[cols_n[1]]]*3,"r")
    if full_prior[cols_n[0],3] == 1: # gaussian error
        ax[i1,i2].plot([ full_prior[cols_n[0],0]-full_prior[cols_n[0],1],full_prior[cols_n[0],0],full_prior[cols_n[0],0]+full_prior[cols_n[0],1] ],[full_prior[cols_n[1],0]]*3,"b")
    else: 
        ax[i1,i2].plot([ full_prior[cols_n[0],1],full_prior[cols_n[0],0],full_prior[cols_n[0],2] ],[full_prior[cols_n[1],0]]*3,"b--") # uniform
    if ((i1==0) & (i2==0)): ax[i1,i2].legend((["all","accepted","post mean","prior mean"]),fontsize=4)
    # if ((i1==0) & (i2==0)): plt.colorbar(cp)
    ax[i1,i2].set_xlabel(cols[0])
    ax[i1,i2].set_ylabel(cols[1])
    i2 = i2+1
    if i2 == 3: i2 = 0; i1 = i1+1; 
fig.tight_layout()
plt.savefig("figs/MCMC_paramspace.png") 
plt.savefig("figs/MCMC_paramspace.pdf") 
fig.show() 
    
#%% Compare results from different step sizes and look at significance
    
nsteps = len(np.unique(stepsize))
step_means = np.zeros((nsteps+1,full_params.shape[1]))*np.nan # means for 0.3, 0.5, 0.75, all, for each post. parameter
step_stds = np.zeros((nsteps+1,full_params.shape[1]))*np.nan
step_ps = np.zeros((nsteps+1,full_params.shape[1]))*np.nan
for n in range(0,nsteps+1):
    if n<nsteps: 
        total = np.where((stepsize == np.unique(stepsize)[n]))
        tmp = np.where((stepsize == np.unique(stepsize)[n]) & (accept_MC_res==1)) # results from this stepsize
        tmp2 = np.where((stepsize != np.unique(stepsize)[n]) & (accept_MC_res==1)) # results from the other step sizes
        t,p = stats.ttest_ind(np.array(full_params)[tmp,:],np.array(full_params)[tmp2,:],axis=1) # ttest: are results from this step different to all other accepted results?
        step_ps[n,:] = p # p of whether this group of stepsize results are sig diff to the rest of the accepted results
        print("stepsize = "+str(np.unique(stepsize)[n])+": runs = "+str(len(total[0]))+", accepted = "+str(len(tmp[0]))+", rate = "+str(len(tmp[0])/len(total[0])*100))
    if n==nsteps+1: 
        tmp = np.where((accept_MC_res==1))
        print("all steps: runs = "+str(len(accept_MC_res))+", accepted = "+str(len(tmp[0]))+", rate = "+str(len(tmp[0])/len(accept_MC_res)*100))
    # Save the means and stds for each step size
    step_means[n,:] = np.nanmean(full_params.iloc[tmp],axis=0)
    step_stds[n,:] = np.nanstd(full_params.iloc[tmp],axis=0)
    
# test the mean and stdev for all values against the prior
for n in range(0,full_params.shape[1]):
    if full_prior[n,3] == 1: # gaussian priors
        t,p = stats.ttest_ind_from_stats(step_means[3,n],step_stds[3,n],len(go),full_prior[n,0],full_prior[n,1],len(go))
        step_ps[nsteps,n] = p # p of whether this group of stepsize results are sig diff to the rest of the accepted results
        print(p) ### note: all gaussian parameters are significantly different to priors
step_ps_summ = step_ps.copy()*0
step_ps_summ[np.where(step_ps<0.01)] = 1
step_ps_summ = pd.DataFrame(step_ps_summ,columns=full_params.columns,index=np.append(np.unique(stepsize),"all_to_prior")).transpose()
print(step_ps_summ)
print("1 = significantly different")

# Plot the means for each step size
labels = np.append(np.unique(stepsize),"all")
w = 0.7 # width of the bars
fig, ax = plt.subplots(5,3,figsize=(8,6))
x,y = 0,0 # Subplot index
for n in range(0,full_params.shape[1]):
    ax[x,y].bar(range(0,4),step_means[:,n],yerr=step_stds[:,n],width=w,linewidth=0,color="r",tick_label=labels)
    low = min(step_means[:,n]-step_stds[:,n])
    high = max(step_means[:,n]+step_stds[:,n])
    ax[x,y].set_ylim([(low-1*(high-low)), (high+1*(high-low))])
    ax[x,y].set_title(full_params.columns[n],fontsize=8)
    ax[x,y].set_xlim([-0.5,3.5])
    ax[x,y].plot((-0.5,3.5),(0,0),"k") # plot 0 line
    ax[x,y].plot((-0.5,3.5),(full_prior[n,0],full_prior[n,0]),"b")
    if full_prior[n,3] == 0:  # uniform error in prior
        ax[x,y].plot((-0.5,3.5),(full_prior[n,1],full_prior[n,1]),"b--")
        ax[x,y].plot((-0.5,3.5),(full_prior[n,2],full_prior[n,2]),"b--")
    if full_prior[n,3] == 1:  # gaussian error in prior
        ax[x,y].plot((-0.5,3.5),(full_prior[n,0]-full_prior[n,1],full_prior[n,0]-full_prior[n,1]),"b:")
        ax[x,y].plot((-0.5,3.5),(full_prior[n,0]+full_prior[n,1],full_prior[n,0]+full_prior[n,1]),"b:")
    x = x+1
    if x==5: x=0; y = y+1
fig.tight_layout()
plt.savefig("figs/MCMC_stepsize-comparison.png") 
plt.savefig("figs/MCMC_stepsize-comparison.pdf") 
fig.show()
    
#%% Look at correlations between the different parameters
    
# select only accepted params and reorder params for this figure
accepted = np.zeros((len(go),full_params.shape[1]))
a_names = ["c_prea","lifetime PD/PI","lifetime PD","d15N ocean","SP ocean","d15N prea","SP prea","frac express","fert EF_red","temp sens","N2 fit scale","NO fit scale","TtoS mass","F ocean"]
for n, v in enumerate(a_names):
    accepted[:,n] = full_params[v].iloc[go]    

# look at correlations
correlations_pearson = np.zeros((len(a_names),len(a_names),2)) # corr, p
correlations_spearman = np.zeros((len(a_names),len(a_names),2)) # corr, p
for n in range(0,len(a_names)):
    for i in range(0,len(a_names)):
        r,p = stats.pearsonr(accepted[:,n],accepted[:,i])
        if n!=i: correlations_pearson[n,i,0:2] = [r,p]
        r,p = stats.spearmanr(accepted[:,n],accepted[:,i])
        if n!=i: correlations_spearman[n,i,0:2] = [r,p]
corrs_p = correlations_pearson[:,:,0].copy() # collect just the correlation coefficients
corrs_p[correlations_pearson[:,:,1]>0.01]=0 # set the "insigificant" correlation coefficients to 0 for the plot
corrs_p_df = pd.DataFrame(corrs_p, index=a_names,columns=a_names)
corrs_s = correlations_spearman[:,:,0].copy()
corrs_s[correlations_spearman[:,:,1]>0.01]=0
corrs_s_df = pd.DataFrame(corrs_s, index=a_names,columns=a_names)

# plot the results
fig, ax = plt.subplots(1,1,figsize=(6,6))
sns.heatmap(corrs_p_df, vmin=-1, vmax=1, center=0,cmap=sns.diverging_palette(220, 20, n=200),square=True)
plt.title("Pearson", fontsize =16)
plt.tight_layout()
plt.savefig("figs/MCMC_corrs_pearson.png") 
plt.savefig("figs/MCMC_corrs_pearson.pdf") 
fig.show()
fig, ax = plt.subplots(1,1,figsize=(6,6))
sns.heatmap(corrs_s_df, vmin=-1, vmax=1, center=0,cmap=sns.diverging_palette(220, 20, n=200),square=True)
plt.title("Spearman", fontsize =16)
plt.tight_layout()
plt.savefig("figs/MCMC_corrs_spearman.png") 
plt.savefig("figs/MCMC_corrs_spearman.pdf") 
fig.show()

#%% PCA to look at multiparam correlations

# standardise the data scaling
full_params_norm = StandardScaler().fit_transform(full_params.iloc[go])
# perform PCA
pca = PCA(n_components=10)
principalComponents = pca.fit_transform(full_params_norm)
colnames = ["PC"+str(n+1) for n in range(0,10)]
principalDf = pd.DataFrame(data = principalComponents, columns = colnames)
exvar = pca.explained_variance_ratio_ # tells how much each PC explains (21, 17, 11, 10%...): cumsum(exvar) for cumulative sum: 8 PCs for 90% of variance
components = pca.components_
# visualise
fig, ax = plt.subplots(2,1,figsize=(8,6))
for n in (0,1):
    ax[n].set_xlabel('Principal Component 1: '+str(exvar[0].round(3)*100)+"%", fontsize = 8)
    ax[n].set_ylabel('Principal Component '+str(n+2)+": "+str(exvar[n+1].round(3)*100)+"%", fontsize = 8)
    #ax.scatter(principalDf.loc[:, 'PC1'], principalDf.loc[:, 'PC2'], s = 50) plot all data if needed
    for i, txt in enumerate(full_params.columns):
        ax[n].annotate(txt, (components[0,i]+0.01, components[n+1,i]+0.01),fontsize=5)
        ax[n].plot((0,components[0,i]), (0,components[n+1,i]),"k")
    ax[n].scatter(components[0,:],components[n+1,:])
plt.tight_layout()
plt.savefig("figs/MCMC_PCA.png") 
plt.savefig("figs/MCMC_PCA.pdf") 
fig.show()

#%% Figures
        
print(str(round(sum(n_acc)/sum(n_its)*100,1))+"% of solutions accepted; ",str(sum(n_acc))+" solutions")
go = np.where(accept_MC_res==1)[0]
post = np.nanmean(prior_solutions[:,go],axis=1)
post_sd = np.nanstd(prior_solutions[:,go],axis=1)

# Final N gas parameterisations
para.gp_finalfig(scale_fitNO = (1,post[10],1,1),scale_fitN2 = (1,post[1],1,1))

# post model
tmp = model(post)
modres_post = expand_modres(tmp)
k_G = np.zeros(d15N_grid.shape)+np.nan; k_G[datarng] = tmp[3][:,2]
f_N2O = np.zeros(d15N_grid.shape)+np.nan; f_N2O[datarng] = tmp[4][:,0]
EF = k_G*f_N2O*100
EF_zones_post = climzone_means(var1_grid = MAT_grid, var2_grid = BulkD_grid, datavar1 = MAT_grid,
                              datavar2 = BulkD_grid,data = EF,bins=4,plotfigs="N",LON=LON,LAT=LAT)
# prior model
tmp = model(x=params_orig[:,0])
modres_prior = expand_modres(tmp)

# Plot conc and isotopes with time
labels = ("N2O trop (ppb)","d15N (permil)","SP (permil)","EF")
fig, ax = plt.subplots(4,2)
for n in range(0,4): # For N2O, d15N, SP, EF
    tmp = obs[:,2]==n+1
    if n<3:
        ax[n,0].plot(tstarts,obs[tmp,0],"bo")
        ax[n,0].plot(tstarts,post_obs[tmp],"cx")
        ax[n,0].plot(tstarts,modres_prior[tmp],"ro")
        ax[n,0].plot(tstarts,modres_post[tmp],"mx")
    if n==3:
        ax[n,0].errorbar(np.arange(1,17),obs[tmp,0],yerr=obs[tmp,1],color="blue")
        ax[n,0].plot(np.arange(1,17),obs[tmp,0],"bo")
        ax[n,0].errorbar(np.arange(1,17),post_obs[tmp],yerr=post_obs_sd[tmp],color="cyan")
        ax[n,0].plot(np.arange(1,17),post_obs[tmp],"cx")
        ax[n,0].plot(np.arange(1,17),modres_prior[tmp],"ro")
        ax[n,0].plot(np.arange(1,17),modres_post[tmp],"mx")
    if n==0: ax[n,0].legend((["obs_prior","obs_post","mod_prior","mod_post"]))
    ax[n,0].set_ylabel(labels[n])
    # plot the 1:1
    ymin = np.nanmin(np.append(obs[tmp,0],modres_prior[tmp]))
    ymax = np.nanmax(np.append(obs[tmp,0],modres_prior[tmp]))
    ax[n,1].plot(obs[tmp,0],modres_prior[tmp],"go")
    ax[n,1].plot(post_obs[tmp],modres_post[tmp],"yx")
    ax[n,1].plot((ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10),(ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10),"k-")
    ax[n,1].set_xlim(ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10)
    ax[n,1].set_ylim(ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10)
    ax[n,1].set_xlabel("Obs: "+labels[n])
    ax[n,1].set_ylabel("Mod: "+labels[n])
    if n==0: ax[n,1].legend((["prior","post"]))
    # agreement stats
    print("RMSE prior, "+labels[n]+" = "+str((np.nanmean((obs[tmp,0] - modres_prior[tmp])**2))**0.5))
    print("RMSE post, "+labels[n]+" = "+str((np.nanmean((post_obs[tmp] - modres_post[tmp])**2))**0.5))
fig.show()

# Print the posterior params in order 
table_order = [0,2,3,4,5,6,7,9,8,11,1,10]
for n in table_order:
    print(f"{full_params.columns[n]+': ' :<17} {str(round(post[n],4)) :<10} +/-  {str(round(post_sd[n],4))}")

