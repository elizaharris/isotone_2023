#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 07:50:20 2020

Run after "FullModel.." and "Read_MCMC.." to generate final figures

@author: elizaharris
"""

#%% Run the final model with full results
post = params_post[:,0]
fullmod_post = model(post,preamble,fullres="Y")
prior = params[:,0]
fullmod_prior = model(prior,preamble,fullres="Y")

# unpack the posterior results
soilmodres = fullmod_post[3]
N_processes_vec = fullmod_post[4]
k_L_post = np.zeros(preamble.d15N_grid.shape)+np.nan; k_L_post[datarng] = soilmodres[:,0]
k_G_post = np.zeros(preamble.d15N_grid.shape)+np.nan; k_G_post[datarng] = soilmodres[:,2]
f_denit_post = np.zeros(preamble.d15N_grid.shape)+np.nan; f_denit_post[datarng] = N_processes_vec[:,3]
f_N2O_post = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2O_post[datarng] = N_processes_vec[:,0]
f_N2_post = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2_post[datarng] = N_processes_vec[:,1]
f_NO_post = np.zeros(preamble.d15N_grid.shape)+np.nan; f_NO_post[datarng] = N_processes_vec[:,2]
d15N_modoutput_post = np.zeros(preamble.d15N_grid.shape)+np.nan; d15N_modoutput_post[datarng] = soilmodres[:,3]
N2O_d15N_post = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_d15N_post[datarng] = soilmodres[:,4]
N2O_SP_post = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_SP_post[datarng] = soilmodres[:,5]

# unpack the prior results
soilmodres = fullmod_prior[3]
N_processes_vec = fullmod_prior[4]
k_L_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; k_L_prior[datarng] = soilmodres[:,0]
k_G_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; k_G_prior[datarng] = soilmodres[:,2]
f_denit_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; f_denit_prior[datarng] = N_processes_vec[:,3]
f_N2O_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2O_prior[datarng] = N_processes_vec[:,0]
d15N_modoutput_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; d15N_modoutput_prior[datarng] = soilmodres[:,3]
N2O_d15N_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_d15N_prior[datarng] = soilmodres[:,4]
N2O_SP_prior = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_SP_prior[datarng] = soilmodres[:,5]

# N input mod results
N_summary = fullmod_post[5]
# 0-6: year, total emissions (Tg N-N2O a-1), emissions from fert, dep, fix, mean d15N, mean SP, EDGAR total, mean d15N, mean SP for soil only
N_summary_full = fullmod_post[6]
# 0-6 : year, N input fert, dep, fix (Tg N a-1), N2O denit only fert, dep, fix (Tg N-N2O a-1), 
# 7-12: N2 fert dep fix, NO fert dep fix, 
# 13-18: NH3 fert dep fix, leach fert dep fix

#%% basics
utils.plot_map(preamble.LON,preamble.LAT,k_G_post*100,"fraction of Nr lost to microbial gas production, post",filename="figs/202502_NewWFPS/f_gas_map",show=1)
utils.plot_map(preamble.LON,preamble.LAT,k_G_post*f_N2O_post*100,"fraction of Nr lost to N2O, post",filename="figs/202502_NewWFPS/f_n2o_map",show=1)
utils.plot_map(preamble.LON,preamble.LAT,k_G_prior*100,"fraction of Nr lost to microbial gas production, post",filename="figs/202502_NewWFPS/f_gas_map_prior",show=1)
utils.plot_map(preamble.LON,preamble.LAT,k_G_prior*f_N2O_prior*100,"fraction of Nr lost to N2O, post",filename="figs/202502_NewWFPS/f_n2o_map_prior",show=1)

# global EF
factor_area = (preamble.area_grid/1000/1000/1000/1000)
print("Mean Arctic EF = "+str(round(np.nansum(k_G_post*f_N2O_post*factor_area)/np.nansum(factor_area)*100,2))+"%")
print("Mean Arctic EF = "+str(round(np.nansum(k_G_prior*f_N2O_prior*factor_area)/np.nansum(factor_area)*100,2))+"% (prior)")
# global EF for inputs
y = N_summary[:,0] == 2022
print("Effective Arctic global EF = "+str(N_summary[y,1]/np.sum(N_summary_full[y,1:4])*100)+"%")
print("Effective Arctic global EF = "+str(fullmod_prior[5][y,1]/np.sum(fullmod_prior[6][y,1:4])*100)+"% (prior)")
# Total permafrost emissions
print("Terrestrial N2O emissions in TgN y-1: 1800 for clim. sens., 1.03 for NH Permafrost from Voigt et al.")
print("My results: 1860 = "+str(N_summary[N_summary[:,0] == 1860,1])+"; 2010 = "+str(N_summary[N_summary[:,0] == 2010,1]))
print("My results: 2020 = "+str(N_summary[N_summary[:,0] == 2020,1]))
print("My results: 2022 = "+str(N_summary[N_summary[:,0] == 2022,1]))

# Emissions against time
fig = plt.figure(figsize=(12,6))
plt.plot(N_summary[:,0],N_summary[:,1],"o",label="total")
plt.plot(N_summary[:,0],N_summary[:,2],"-",label="fert")
plt.plot(N_summary[:,0],N_summary[:,3],":",label="dep")
plt.plot(N_summary[:,0],N_summary[:,4],"-.",label="fix")
plt.ylabel("N2O emissions (Tg N2O-N a-1)")
plt.legend()
fig.tight_layout()
plt.savefig("figs/202502_NewWFPS/total_cat_emissions.png") 
plt.savefig("figs/202502_NewWFPS/total_cat_emissions.pdf") 
fig.show() 

# Get map of emissions for a certain year
def getmap(y,preamble,post,k_G_post,f_N2O_post):
    y = 2022
    # get T sens
    tindex = int(np.where(abs(preamble.Tanom.variables["time"][:]-y) == np.min(abs(preamble.Tanom.variables["time"][:]-y)))[0])
    tanomaly = preamble.Tnorm[tindex,:,:]
    tanomaly_sens = (post[4]-1)*tanomaly + 1
    # get N input data (g m-2 y-1)
    index = np.where(preamble.Ninputs.variables["time"][:].data == y)[0] # select the year
    fix = ( preamble.Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens )[0,:,:]
    dep = ( preamble.Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens  )[0,:,:]
    fert = ( preamble.Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] )[0,:,:]
    total = fix + dep + fert
    n2o_total = k_G_post*f_N2O_post*total
    return(n2o_total, fix, dep, fert)
emissions_2022, _, _, _ = getmap(2022,preamble,post,k_G_post,f_N2O_post)
emissions_1800, _, _, _ = getmap(1800,preamble,post,k_G_post,f_N2O_post)

utils.plot_map(preamble.LON,preamble.LAT,emissions_2022,"N2O (g-N m-2 a-1)",filename="figs/202502_NewWFPS/n2o_map",show=1)
utils.plot_map(preamble.LON,preamble.LAT,emissions_2022-emissions_1800,"N2O anthrop (g-N m-2 a-1)",filename="figs/202502_NewWFPS/n2o_map_anthrop",show=1)
vminmax = (-8,-1)
utils.plot_map(preamble.LON,preamble.LAT,np.log(emissions_2022).clip(vminmax[0],vminmax[1]),"N2O (log, g-N m-2 a-1)",vminmax=vminmax,filename="figs/202502_NewWFPS/n2o_map_log",show=1)
utils.plot_map(preamble.LON,preamble.LAT,np.log(emissions_2022-emissions_1800).clip(vminmax[0],vminmax[1]),"N2O, anthrop (log, g-N m-2 a-1)",vminmax=vminmax,filename="figs/202502_NewWFPS/n2o_map_anthrop_log",show=1)

# Emissions with latitude and time
from matplotlib import colormaps
years_to_plot = np.arange(1850,2030,10)
colours_to_plot = colormaps['cool'].resampled(len(years_to_plot))
emissions_grid = k_G_post*f_N2O_post*total
fig = plt.figure(figsize=(12,6))
for n,y in enumerate(years_to_plot):
    # get T sens
    tindex = int(np.where(abs(preamble.Tanom.variables["time"][:]-y) == np.min(abs(preamble.Tanom.variables["time"][:]-y)))[0])
    tanomaly = preamble.Tnorm[tindex,:,:]
    tanomaly_sens = (post[4]-1)*tanomaly + 1
    # get data (g m-2 y-1)
    index = np.where(preamble.Ninputs.variables["time"][:].data == y)[0] # select the year
    fix = ( preamble.Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens )[0,:,:]
    dep = ( preamble.Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens  )[0,:,:]
    fert = ( preamble.Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] )[0,:,:]
    total = fix + dep + fert
    # Find emissions by lat and plot
    emissions = k_G_post*f_N2O_post*total*preamble.area_grid/1000/1000 # tonnes (10^6 grams) y-1
    emissions_by_lat = np.nansum(emissions,axis=1)
    colour_location = (y-np.nanmin(years_to_plot))/(np.nanmax(years_to_plot)-np.nanmin(years_to_plot)) # chosen year normalised to 0-1
    plt.plot(preamble.LAT[:,0],emissions_by_lat,"-",color=colours_to_plot(colour_location),label=y)
plt.ylabel("N2O emissions (10^6 g N2O-N a-1)")
plt.legend()
fig.tight_layout()
fig.show() 
plt.savefig("figs/202502_NewWFPS/emissions_by_lat.png") 
plt.savefig("figs/202502_NewWFPS/emissions_by_lat.pdf") 

# f_n2o against d15N_grid
x = preamble.d15N_grid[datarng]
y = (k_G_post*f_N2O_post)[datarng]
df = pd.DataFrame({"d15N":x,"fN2O":y})
df["d15N_round"] = np.round(df["d15N"],1)
d15N_fbyd = np.sort(np.unique(df["d15N_round"]))
mean_fbyd = df.groupby("d15N_round")["fN2O"].mean()
sd_fbyd = df.groupby("d15N_round")["fN2O"].std()
fig = plt.figure(figsize=(12,6))
plt.plot(x,y,"o",alpha=0.02)
plt.plot(d15N_fbyd,mean_fbyd,"m-")
plt.plot(d15N_fbyd,mean_fbyd-sd_fbyd,"m:")
plt.plot(d15N_fbyd,mean_fbyd+sd_fbyd,"m:")
plt.ylabel("d15N (permil)")
plt.ylabel("f_N2O")
fig.tight_layout()
fig.show() 
plt.savefig("figs/202502_NewWFPS/d15N_fN2O.png") 
plt.savefig("figs/202502_NewWFPS/d15N_fN2O.pdf") 

#%% Compare to the Voigt flux dataset

voigt = pd.read_csv("data/N2OFluxes_Voigt-etal_2020_dk20250212.csv")

voigt["flux_model_year"] = np.nan # space for flux 
voigt["flux_model_1800"] = np.nan # space for flux in 1800 (preindustrial)
voigt["flux_model_2020"] = np.nan # space for flux in 2020 (current)
voigt["flux_model_lon_mismatch"] = np.nan # save the closeness of lat and lon match
voigt["flux_model_lat_mismatch"] = np.nan # 
# Get fluxes
lats = preamble.LAT[:,0]
lons = preamble.LON[0,:]
emissions_1800, _, _, _ = getmap(1800,preamble,post,k_G_post,f_N2O_post)
emissions_2020, _, _, _ = getmap(2020,preamble,post,k_G_post,f_N2O_post)
for n in np.arange(voigt.shape[0]):
    print(n)
    # Find the matching cell
    voigt.loc[n,"flux_model_lat_mismatch"] = np.nanmin(abs(lats-voigt["latitude"].iloc[n]))
    i = np.where( abs(lats-voigt["latitude"].iloc[n]) == np.nanmin(abs(lats-voigt["latitude"].iloc[n])) )[0]
    voigt.loc[n,"flux_model_lon_mismatch"] = np.nanmin(abs(lons-voigt["longitude"].iloc[n]))
    j = np.where( abs(lons-voigt["longitude"].iloc[n]) == np.nanmin(abs(lons-voigt["longitude"].iloc[n])) )[0]
    # 1800 and 2020 fluxes
    voigt.loc[n,"flux_model_1800"] = emissions_1800[i,j]
    voigt.loc[n,"flux_model_2020"] = emissions_2020[i,j]
    # matching year flux
    y = np.nanmean([voigt['year_start'].iloc[n],voigt['year_end'].iloc[n]])
    emissions_y, _, _, _ = getmap(2020,preamble,post,k_G_post,f_N2O_post)
    voigt.loc[n,"flux_model_year"] = emissions_y[i,j]

# Plot and compare
lat_res = np.nanmean(np.diff(lats))
lon_res = np.nanmean(np.diff(lons))
r = np.where( (voigt["flux_model_lat_mismatch"]<2*lat_res) & (voigt["flux_model_lon_mismatch"]<2*lon_res) )[0]

plt.plot(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r]/10**6,voigt["flux_model_year"].iloc[r],"x")
plt.show()

plt.plot(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r]/10**6,voigt["flux_model_year"].iloc[r],"x")
plt.xlim([0,0.3])
plt.show()

r_veg = [ rr for rr in r if voigt['veg_type'].iloc[rr]=="Vegetated" ]
r_pveg = [ rr for rr in r if voigt['veg_type'].iloc[rr]=="Partly vegetated" ]
r_bare = [ rr for rr in r if voigt['veg_type'].iloc[rr]=="Bare" ]
plt.plot(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r_veg]/10**6,voigt["flux_model_year"].iloc[r_veg],"go",label="Vegetated")
plt.plot(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r_pveg]/10**6,voigt["flux_model_year"].iloc[r_pveg],"bo",label="Partly vegetated")
plt.plot(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r_bare]/10**6,voigt["flux_model_year"].iloc[r_bare],"ro",label="Bare")
#plt.xlim([0,0.3])
plt.show()

plt.plot(np.log(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r]/10**6),np.log(voigt["flux_model_year"].iloc[r]),"x")
plt.show()

from sklearn.linear_model import LinearRegression
X = np.array(voigt["annual_N2O_flux_ug_m2_100days-x-2"].iloc[r]/10**6)
y = np.array(voigt["flux_model_year"].iloc[r])
rr = [ not np.isnan(yy) for yy in y ]
reg = LinearRegression().fit(X[rr].reshape(-1, 1), y[rr])
reg.score(X[rr].reshape(-1, 1), y[rr])

#### FROM global version; incorporate some of these later...

#%% Estimate uncertainty in posterior model by running multiple times using randomly selected accepted solutions

n_its = 100
rand = [int(x) for x in np.random.uniform(low=0, high=len(go), size=n_its)] # randomly select solutions
soilmodres_all = np.zeros((soilmodres.shape[0],soilmodres.shape[1],n_its))+np.nan
N_processes_vec_all = np.zeros((N_processes_vec.shape[0],N_processes_vec.shape[1],n_its))+np.nan
N_summary_all = np.zeros((N_summary.shape[0],N_summary.shape[1],n_its))+np.nan
N_summary_full_all = np.zeros((N_summary_full.shape[0],N_summary_full.shape[1],n_its))+np.nan
for its in np.arange(n_its):
    if its%5==0: print(its)
    post_i = prior_solutions[:,go[rand][its]]
    fullmod_i = model(post_i,fullres="Y")

    # unpack the posterior results
    soilmodres_all[:,:,its] = fullmod_i[3]
    N_processes_vec_all[:,:,its] = fullmod_i[4]
    N_summary_all[:,:,its] = fullmod_i[5]
    N_summary_full_all[:,:,its] = fullmod_i[6]

soilmodres_sd = np.nanstd(soilmodres_all,axis=2)
N_processes_vec_sd = np.nanstd(N_processes_vec_all,axis=2)
N_summary_sd = np.nanstd(N_summary_all,axis=2)
N_summary_full_sd = np.nanstd(N_summary_full_all,axis=2)

k_L_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; k_L_sd[datarng] = soilmodres_sd[:,0]
k_G_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; k_G_sd[datarng] = soilmodres_sd[:,2]
f_denit_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; f_denit_sd[datarng] = N_processes_vec_sd[:,3]
f_N2O_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2O_sd[datarng] = N_processes_vec_sd[:,0]
f_N2_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; f_N2_sd[datarng] = N_processes_vec_sd[:,1]
f_NO_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; f_NO_sd[datarng] = N_processes_vec_sd[:,2]
d15N_modoutput_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; d15N_modoutput_sd[datarng] = soilmodres_sd[:,3]
N2O_d15N_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_d15N_sd[datarng] = soilmodres_sd[:,4]
N2O_SP_sd = np.zeros(preamble.d15N_grid.shape)+np.nan; N2O_SP_sd[datarng] = soilmodres_sd[:,5]

#%% basics

plot_map(LON,LAT,k_G_post*f_N2O_post*(1-f_denit_post)*100,"fraction of Nr lost to nitrification N2O")
plot_map(LON,LAT,k_G_post*f_N2O_post*100,"fraction of Nr lost to N2O, post")

# global EF
factor_area = (area_grid/1000/1000/1000/1000)
print("Mean global EF = "+str(round(np.nansum(k_G_post*f_N2O_post*factor_area)/np.nansum(factor_area)*100,2))+"%")
err = ((k_G_sd/k_G_post)**2 + (f_N2O_sd/f_N2O_post)**2)**0.5
err[np.isinf(err)] = np.nan
print("Err in mean global EF = "+str(round(np.nansum(err*k_G_post*f_N2O_post*factor_area)/np.nansum(factor_area)*100,2))+"%")
print("Mean global EF = "+str(round(np.nansum(k_G_prior*f_N2O_prior*factor_area)/np.nansum(factor_area)*100,2))+"% (prior)")
# global EF for inputs
y = N_summary[:,0] == 1940
print("Effective global EF = "+str(N_summary[y,1]/np.sum(N_summary_full[y,1:4])*100)+"%")
err = N_summary_sd[y,1]/N_summary[y,1]
print("Err in effective global EF = "+str(err*N_summary[y,1]/np.sum(N_summary_full[y,1:4])*100)+"%")
print("Effective global EF = "+str(fullmod_prior[5][y,1]/np.sum(fullmod_prior[6][y,1:4])*100)+"% (prior)")
# comparison:
print("Terrestrial N2O emissions in TgN y-1: 1800 for clim. sens., 1860 6.3pm1.1, 2010 10pm2.2 in Tian2019")
print("My results: 1860 = "+str(N_summary[N_summary[:,0] == 1860,1])+"; 2010 = "+str(N_summary[N_summary[:,0] == 2010,1]))
print("My results: 2020 = "+str(N_summary[N_summary[:,0] == 2020,1]))
# NO emissions
print("My results, NO, 2010 = "+str(sum(N_summary_full[N_summary[:,0] == 2010,10:13])))

# EFs for climate zones
EF_zones_post = utils.climzone_means(var1_grid = MAT_grid, var2_grid = BulkD_grid, datavar1 = MAT_grid,
                              datavar2 = BulkD_grid,data = k_G_post*f_N2O_post*100,bins=4,plotfigs="Y",LON=LON,LAT=LAT)

#%% Plot k_G and f_N2O with uncertainties

kG_post_fromallres = np.zeros(preamble.d15N_grid.shape)+np.nan; kG_post_fromallres[datarng] = post_kG_m
kG_post_sd_fromallres = np.zeros(preamble.d15N_grid.shape)+np.nan; kG_post_sd_fromallres[datarng] = post_kG_sd
fN2O_post_fromallres = np.zeros(preamble.d15N_grid.shape)+np.nan; fN2O_post_fromallres[datarng] = post_fN2O_m
fN2O_post_sd_fromallres = np.zeros(preamble.d15N_grid.shape)+np.nan; fN2O_post_sd_fromallres[datarng] = post_fN2O_sd
# note: these fields are almost identical to the ones got from the mean post, then the model run!
tmp = ((post_kG_sd/post_kG_m)**2 + (post_fN2O_sd/post_fN2O_m)**2)**0.5 * post_kG_m * post_fN2O_m
kGfN2O_post_sd_fromallres = np.zeros(preamble.d15N_grid.shape)+np.nan; kGfN2O_post_sd_fromallres[datarng] = tmp

plot_map(LON,LAT,kG_post_fromallres*100,"fraction of Nr lost to gas production")
plot_map(LON,LAT,kG_post_sd_fromallres*100,"fraction of Nr lost to gas production, uncertainty")
plot_map(LON,LAT,kG_post_fromallres*fN2O_post_fromallres*100,"fraction of Nr lost to N2O")
plot_map(LON,LAT,kGfN2O_post_sd_fromallres*100,"fraction of Nr lost to N2O, uncertainty")

fig, ax = plt.subplots(1,1)
a=ax.scatter(MAP_grid[datarng],(kG_post_fromallres*100)[datarng],c=WFPS_grid[datarng],marker=".",cmap="coolwarm")
ax.set_xlim(0,5000)
cbar = plt.colorbar(a,fraction=0.016, pad=0.04)

#%% Plot conc and isotopes with time for post and prior

# cyan and magenta show obs error ranges!
modres_post = fullmod_post[0]
modres_prior = fullmod_prior[0]
labels = ("N2O trop (ppb)","d15N (permil)","SP (permil)")
fig, ax = plt.subplots(3,2)
for n in range(0,3):
    istart = np.min(np.where(obs[:,2]==n+1))
    iend = np.max(np.where(obs[:,2]==n+1))+1
    ax[n,0].plot(tstarts,obs[istart:iend,0],"bo")
    ax[n,0].plot(tstarts,post_obs[istart:iend],"cx")
    ax[n,0].plot(tstarts,modres_prior[istart:iend],"ro")
    ax[n,0].plot(tstarts,modres_post[istart:iend],"mx")
    ax[n,0].plot(tstarts,obs[istart:iend,0]+2*obs[istart:iend,1],"b:")
    ax[n,0].plot(tstarts,obs[istart:iend,0]-2*obs[istart:iend,1],"b:")
    if n==0: ax[n,0].legend((["obs_prior","obs_post","mod_prior","mod_post"]))
    ax[n,0].set_ylabel(labels[n])
    ax[n,1].plot(obs[istart:iend,0],modres_prior[istart:iend],"go")
    ax[n,1].plot(post_obs[istart:iend],modres_post[istart:iend],"yx")
    ax[n,1].plot(post_obs[istart:iend]+2*obs[istart:iend,1],modres_post[istart:iend],"m.")
    ax[n,1].plot(post_obs[istart:iend]-2*obs[istart:iend,1],modres_post[istart:iend],"c.")
    ax[n,1].plot(modres_post[istart:iend],modres_post[istart:iend],"k-")
    ax[n,1].set_xlabel("Obs: "+labels[n])
    ax[n,1].set_ylabel("Mod: "+labels[n])
    if n==0: ax[n,1].legend((["prior","post"]))
fig.show()

#%% plot total and anthrop flux for a given year

tmp = 2019
# get T sens
tindex = int(np.where(abs(Tanom.variables["time"][:]-tmp) == np.min(abs(Tanom.variables["time"][:]-tmp)))[0])
tanomaly = Tnorm[tindex,:,:]
tanomaly_sens = (post[4]-1)*tanomaly + 1
# get data (g m-2 y-1)
index = np.where(Ninputs.variables["time"][:].data == tmp)[0] # select the year
fix = ( Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens )[0,:,:]
dep = ( Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens  )[0,:,:]
fert = ( Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] )[0,:,:]

tmp = 1800 # base year also
# get T sens
tindex = int(np.where(abs(Tanom.variables["time"][:]-tmp) == np.min(abs(Tanom.variables["time"][:]-tmp)))[0])
tanomaly = Tnorm[tindex,:,:]
tanomaly_sens = (post[4]-1)*tanomaly + 1
# get data (g m-2 y-1)
index = int(np.where(abs(Ninputs.variables["time"][:]-tmp) == np.min(abs(Ninputs.variables["time"][:]-tmp)))[0])
fix1800 = Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens 
dep1800 = Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens
fert1800 = Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] 

vminmax = (-8,-1)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*fix - k_G_post*f_N2O_post*fix1800).clip(vminmax[0],vminmax[1]),"N2O anth, fix (g-N m-2 a-1)",vminmax=vminmax)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*dep - k_G_post*f_N2O_post*dep1800).clip(vminmax[0],vminmax[1]),"N2O anth, dep (g-N m-2 a-1)",vminmax=vminmax)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*fert - k_G_post*f_N2O_post*fert1800).clip(vminmax[0],vminmax[1]),"N2O anth, fert (g-N m-2 a-1)",vminmax=vminmax)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*fix).clip(vminmax[0],vminmax[1]),"N2O total, fix (g-N m-2 a-1",vminmax=vminmax)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*dep).clip(vminmax[0],vminmax[1]),"N2O total, dep (g-N m-2 a-1",vminmax=vminmax)
plot_map(LON,LAT,np.log(k_G_post*f_N2O_post*fert).clip(vminmax[0],vminmax[1]),"N2O total, fert (g-N m-2 a-1",vminmax=vminmax)

# Compare to Rona's emissions
isotonemissions = k_G_post*f_N2O_post*(fix+dep+fert*post[3])
isotonemissions_sd = (( (k_G_sd/k_G_post)**2 + (f_N2O_sd/f_N2O_post)**2)**0.5)  *(fix+dep+fert*post[3])
tmp = np.isnan(isotonemissions)
inversemissions = np.array(nc4.Dataset('Emissions_Rona.nc','r').variables['annualN2O'])
inversemissions[tmp] = np.nan
inversemissions_sd = np.array(nc4.Dataset('Emissions_Rona.nc','r').variables['sdN2O'])
inversemissions_sd[tmp] = np.nan
if 0: # log plots
    vminmax = (-8,-1)
    plot_map(LON,LAT,np.log(isotonemissions).clip(vminmax[0],vminmax[1]),"2019 emissions, IsoTONE (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,np.log(inversemissions).clip(vminmax[0],vminmax[1]),"2019 emissions, inversion (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,np.log(isotonemissions_sd).clip(vminmax[0],vminmax[1]),"Err: 2019 emissions, IsoTONE (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,np.log(inversemissions_sd).clip(vminmax[0],vminmax[1]),"Err: 2019 emissions, inversion (g-N m-2 a-1)",vminmax=vminmax)
if 0: # lin plots
    vminmax = (0,1)
    plot_map(LON,LAT,(isotonemissions).clip(vminmax[0],vminmax[1]),"2019 emissions, IsoTONE (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,(inversemissions).clip(vminmax[0],vminmax[1]),"2019 emissions, inversion (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,(isotonemissions_sd).clip(vminmax[0],vminmax[1]),"Err: 2019 emissions, IsoTONE (g-N m-2 a-1)",vminmax=vminmax)
    plot_map(LON,LAT,(inversemissions_sd).clip(vminmax[0],vminmax[1]),"Err: 2019 emissions, inversion (g-N m-2 a-1)",vminmax=vminmax)
# Plot the difference also
diff = isotonemissions - inversemissions
diff_sd = ( (isotonemissions_sd)**2 + (inversemissions_sd)**2 )**0.5
plot_map(LON,LAT,diff,"2019 emissions, difference (g-N m-2 a-1)",cmap="bwr")
tmp = diff.copy(); tmp[abs(diff) <= diff_sd] = np.nan
plot_map(LON,LAT,tmp,"2019 emissions, difference (g-N m-2 a-1)",cmap="bwr")

# stats of comparison
tmp = ~np.isnan(isotonemissions)
gpy = stats.linregress(isotonemissions[tmp],inversemissions[tmp])
fig, ax = subplots(1)
ax.plot(isotonemissions[tmp],inversemissions[tmp],marker="o")
ax.plot(isotonemissions[tmp],isotonemissions[tmp]*gpy[0] + gpy[1],"r-")
# Significant correlation but far from 1:1
r,p = stats.spearmanr(isotonemissions[tmp],inversemissions[tmp])

# Total emissions from inversion
total = np.nansum(inversemissions[datarng]*factor_area)
total_sd = (np.nansum((inversemissions_sd[datarng]*factor_area)**2))**0.5

#%% emission factor comparison with Cui database (agri) and Dorich database (natural)

# Get the Cui data, extract columns, rename
EFs_Cui = pd.read_csv("input_data/Cui2021_EFdata/Cui2021_EFdata.csv") # Longitude, Latitude, E, E0, EF
for n in np.arange(EFs_Cui.shape[0]):
    if not(isinstance(EFs_Cui["E"].iloc[n],str)): continue
    if len(EFs_Cui["E"].iloc[n]) == 6:
        tmp = EFs_Cui["E"].iloc[n][1:5]
    else: tmp = EFs_Cui["E"].iloc[n]
    EFs_Cui["E"].iloc[n] = np.float(tmp)
    if not(isinstance(EFs_Cui["EF"].iloc[n],str)): continue
    if len(EFs_Cui["EF"].iloc[n]) == 6:
        tmp = EFs_Cui["EF"].iloc[n][1:5]
    else: tmp = EFs_Cui["EF"].iloc[n]
    EFs_Cui["EF"].iloc[n] = np.float(tmp)
EF_summary = EFs_Cui[["Latitude","Longitude","Crop type","E","EF","EF_sd"]]
EF_summary.columns = ("Lat","Lon","Crop","N2O_kgN2O-N_ha","EF","EF_sd")
EF_summary["Database"] = "Cui"

# Get the Dorich data, extract columns, rename
EFs_Dorich = pd.read_csv("input_data/N2O_ChrisDorich/globaln2o_sites_filled.csv")
EFs_Dorich = EFs_Dorich[['Lat', 'Long','Man_Crop','N2O_kgN2O-N_ha','EF']]
EFs_Dorich["EF_sd"] = [ max(EF*0.2,0.2) for EF in EFs_Dorich["EF"] ] # no sd given so assume 20%
EFs_Dorich.columns = ("Lat","Lon","Crop","N2O_kgN2O-N_ha","EF","EF_sd")
EFs_Dorich["Database"] = "Dorich"

EF_summary = EF_summary.append(EFs_Dorich,ignore_index=True)
EF_summary["EF"] = pd.to_numeric(EF_summary["EF"])
EF_summary["N2O_kgN2O-N_ha"] = pd.to_numeric(EF_summary["N2O_kgN2O-N_ha"])

# Get the 2020 effective EFs from the model (accounting for fert EF red)
emissions = k_G_post * f_N2O_post * (fix + dep) + k_G_post * f_N2O_post * post[3] * fert # g N2O-N m-2 a-1
EF = emissions / (fix + dep + fert)
EF_sd = ( ( (k_G_sd/k_G_post)**2 + (f_N2O_sd/f_N2O_post)**2 )**0.5 ) * EF

# Match a grid cell value to each Cui/Dorich data point
EF_summary[["model_fgas","model_fN2O","model_EF","model_EF_eff","model_N2O"]] = np.nan # Add columns for comparable data from model
for n in np.arange(EF_summary.shape[0]):
    if np.isnan(EF_summary["Lon"].iloc[n]): continue
    x = np.where(abs(LON[0,:] - EF_summary["Lon"].iloc[n]) == min(abs(LON[0,:] - EF_summary["Lon"].iloc[n])))[0]
    y = np.where(abs(LAT[:,0] - EF_summary["Lat"].iloc[n]) == min(abs(LAT[:,0] - EF_summary["Lat"].iloc[n])))[0]
    if len(x)>1: x = x[0]
    if len(y)>1: y = y[0]
    EF_summary["model_fgas"].iloc[n] = k_G_post[y,x]
    EF_summary["model_fN2O"].iloc[n] = f_N2O_post[y,x]
    EF_summary["model_EF"].iloc[n] = k_G_post[y,x]*f_N2O_post[y,x]*100
    EF_summary["model_EF_eff"].iloc[n] = EF[y,x]*100
    EF_summary["model_N2O"].iloc[n] = emissions[y,x]/1000*10000
    
plt.errorbar(EF_summary["model_EF_eff"],EF_summary["EF"],EF_summary["EF_sd"],fmt='o')
plt.scatter(EF_summary["model_N2O"],EF_summary["N2O_kgN2O-N_ha"],marker="o")
  
# Find the mean and std EF for each gridcell for Cui data
tmp = np.where(~np.isnan(EF))
EF_summary_gridmatch = np.zeros((len(tmp[0]),4)) + np.nan # EF, EFsd, N2O, N2Osd
for n in np.arange(len(tmp[0])):
    match1 = (EF_summary["Lon"] >= LON[tmp[0][n],tmp[1][n]]-0.25) & (EF_summary["Lon"] <= LON[tmp[0][n],tmp[1][n]]+0.25)
    match2 = (EF_summary["Lat"] >= LAT[tmp[0][n],tmp[1][n]]-0.25) & (EF_summary["Lat"] <= LAT[tmp[0][n],tmp[1][n]]+0.25)
    match = np.where(match1 & match2)[0]
    if len(match)==1:
        EF_summary_gridmatch[n,0] = EF_summary["EF"].iloc[match]
        EF_summary_gridmatch[n,1] = EF_summary["EF_sd"].iloc[match]
        EF_summary_gridmatch[n,2] = EF_summary["N2O_kgN2O-N_ha"].iloc[match]
    if len(match)>1:
        EF_summary_gridmatch[n,0] = np.nanmean(pd.to_numeric(EF_summary["EF"].iloc[match]))
        EF_summary_gridmatch[n,1] = np.nanstd(pd.to_numeric(EF_summary["EF"].iloc[match]))
        EF_summary_gridmatch[n,2] = np.nanmean(pd.to_numeric(EF_summary["N2O_kgN2O-N_ha"].iloc[match]))
        EF_summary_gridmatch[n,3] = np.nanstd(pd.to_numeric(EF_summary["N2O_kgN2O-N_ha"].iloc[match]))
EF_summary_grid = EF.copy()
EF_summary_grid[tmp] = EF_summary_gridmatch[:,0]

# note! Only a VERY small proportion of gridcells have obs data; 316 out of 55478 grid cells with a modelled EF
plot_map(LON,LAT,EF_summary_grid,"fraction of Nr lost to N2O, observations")
plot_map(LON,LAT,EF,"fraction of Nr lost to N2O, post")

# Point to point comparison of Gaussian intersection probabilities
diff = EF_summary_gridmatch[:,0] - EF[tmp]*100
err = (EF_summary_gridmatch[:,1]**2 + (EF_sd[tmp]*100)**2)**0.5
gauss_prod = (diff / err)
gp = np.exp(-1/8*gauss_prod)
plt.plot(diff,"o")
plt.plot(err,"o")
# The difference is much larger than the error in nearly all cases

# Correlation?
nxa = np.where((~np.isnan(EF)) & (~np.isnan(EF_summary_grid))) # linear regression - need to add x and y error!
gpy = stats.linregress(EF[nxa]*100,EF_summary_grid[nxa])
fig, ax = subplots(1)
ax.errorbar(EF[tmp]*100,EF_summary_gridmatch[:,0],EF_summary_gridmatch[:,1],fmt="o")
ax.plot(EF[tmp]*100,EF[tmp]*100*gpy[0] + gpy[1],"r-")
# Significant correlation but far from 1:1
r,p = stats.spearmanr(EF[nxa]*100,EF_summary_grid[nxa])

# Note: This paper shows no significant relationship between application rate and EF (shown by the non-linearity term)
mean(EFs_Cui['ΔEF'])
# Overall, seems like the two correlate (still need to add error to regression) but do not agree - the model
# results are consistently higher.

#%% The anthropogenic flux with time

# 0-6: year, total emissions (Tg N-N2O a-1), emissions from fert, dep, fix, mean d15N, mean SP 
Fanth = N_summary.copy()
Fanth[:,1:] = 0
Fanth[:,1] = N_summary[:,1]-N_summary[0,1] # anth total emissions
Fanth[:,2] = N_summary[:,2]-N_summary[0,2] # anth fert emissions
Fanth[:,3] = N_summary[:,3]-N_summary[0,3] # anth dep emissions
Fanth[:,4] = N_summary[:,4]-N_summary[0,4] # anth fix emissions
frac_anth = Fanth[:,1]/(Fanth[:,1]+N_summary[0,1])
Fanth[:,5] = (N_summary[:,5] - N_summary[0,5]*(1-frac_anth))/frac_anth 
Fanth[frac_anth<0.1,5] = np.nan # at least 10% to calculate d15N
Fanth[:,6] = (N_summary[:,6] - N_summary[0,6]*(1-frac_anth))/frac_anth 
Fanth[frac_anth<0.1,6] = np.nan # at least 10% to calculate SP

Fanth_sd = N_summary.copy()
Fanth_sd[:,1] = ((N_summary_sd[:,1]**2 + N_summary_sd[0,1]**2)/2)**0.5 # anth total emissions
Fanth_sd[:,2] = ((N_summary_sd[:,2]**2 + N_summary_sd[0,2]**2)/2)**0.5# anth fert emissions
Fanth_sd[:,3] = ((N_summary_sd[:,3]**2 + N_summary_sd[0,3]**2)/2)**0.5# anth dep emissions
Fanth_sd[:,4] = ((N_summary_sd[:,4]**2 + N_summary_sd[0,4]**2)/2)**0.5 # anth fix emissions
Fanth_sd[:,5] = N_summary_sd[:,5]/2
Fanth_sd[frac_anth<0.1,5] = np.nan # at least 10% to calculate d15N
Fanth_sd[:,6] = N_summary_sd[:,6]/2
Fanth_sd[frac_anth<0.1,6] = np.nan # at least 10% to calculate SP

fig, ax = plt.subplots(3,1)
ax[0].plot(Fanth[:,0],Fanth[:,1])
ax[0].plot(Fanth[:,0],Fanth[:,2])
ax[0].plot(Fanth[:,0],Fanth[:,3])
ax[0].plot(Fanth[:,0],Fanth[:,4])
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,3])
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,4])
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,1])
ax[0].plot(Fanth[:,0],Fanth[:,1]-Fanth_sd[:,1],":",c="C0")
ax[0].plot(Fanth[:,0],Fanth[:,2]-Fanth_sd[:,1],":",c="C1")
ax[0].plot(Fanth[:,0],Fanth[:,3]-Fanth_sd[:,1],":",c="C2")
ax[0].plot(Fanth[:,0],Fanth[:,4]-Fanth_sd[:,1],":",c="C3")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,3]-N_summary_sd[0,3],":",c="C4")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,4]-N_summary_sd[0,4],":",c="C5")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,1]-N_summary_sd[0,1],":",c="C6")
ax[0].plot(Fanth[:,0],Fanth[:,1]+Fanth_sd[:,1],":",c="C0")
ax[0].plot(Fanth[:,0],Fanth[:,2]+Fanth_sd[:,1],":",c="C1")
ax[0].plot(Fanth[:,0],Fanth[:,3]+Fanth_sd[:,1],":",c="C2")
ax[0].plot(Fanth[:,0],Fanth[:,4]+Fanth_sd[:,1],":",c="C3")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,3]+N_summary_sd[0,3],":",c="C4")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,4]+N_summary_sd[0,4],":",c="C5")
ax[0].plot(Fanth[:,0],Fanth[:,1]*0+N_summary[0,1]+N_summary_sd[0,1],":",c="C6")
ax[0].legend((["total anth (Tg-N a-1)","fert","dep","fix","natural_dep","natural_fix","natural"]))
ax[1].plot(N_summary[:,0],N_summary[:,5])
ax[1].plot(Fanth[:,0],Fanth[:,5])
ax[1].plot(N_summary[:,0],N_summary[:,5]-N_summary_sd[:,5]/2,":",c="C0")
ax[1].plot(Fanth[:,0],Fanth[:,5]-Fanth_sd[:,5],":",c="C1")
ax[1].plot(N_summary[:,0],N_summary[:,5]+N_summary_sd[:,5]/2,":",c="C0")
ax[1].plot(Fanth[:,0],Fanth[:,5]+Fanth_sd[:,5],":",c="C1")
ax[1].legend((["terr total flux","anthrop flux"]))
ax[2].plot(N_summary[:,0],N_summary[:,6])
ax[2].plot(Fanth[:,0],Fanth[:,6])
ax[2].plot(N_summary[:,0],N_summary[:,6]-N_summary_sd[:,6]/2,":",c="C0")
ax[2].plot(Fanth[:,0],Fanth[:,6]-Fanth_sd[:,6],":",c="C1")
ax[2].plot(N_summary[:,0],N_summary[:,6]+N_summary_sd[:,6]/2,":",c="C0")
ax[2].plot(Fanth[:,0],Fanth[:,6]+Fanth_sd[:,6],":",c="C1")
ax[2].legend((["terr total flux","anthrop flux"]))


#%% total N losses to nitrification and denitrification

# map
Nloss_nit = (fix + dep + fert*post[3])*k_G_post*(f_N2O_post*(1-f_denit_post)+f_NO_post)
Nloss_denit = (fix + dep + fert*post[3])*k_G_post*(f_N2O_post*f_denit_post+f_N2_post)
plot_map(LON,LAT,np.log(Nloss_nit).clip(-8,2),"N lost to nitrification (g-N m-2 a-1)")
plot_map(LON,LAT,np.log(Nloss_denit).clip(-8,2),"N lost to denitrification (g-N m-2 a-1)")

# time series
nit_N2O = N_summary[:,1]-np.nansum(N_summary_full[:,4:7],axis=1)-(Fanth[:,1]-np.nansum(Fanth[:,2:5],axis=1))
denit_N2O = np.nansum(N_summary_full[:,4:7],axis=1)
nit = np.nansum(N_summary_full[:,10:13],axis=1)+nit_N2O
denit = np.nansum(N_summary_full[:,7:10],axis=1)+denit_N2O
nit_sd = np.column_stack((N_summary_full_sd[:,10:13],N_summary_sd[:,1]))
nit_sd = (np.nansum(nit_sd**2,axis=1))**0.5 
denit_sd = np.column_stack((N_summary_full_sd[:,7:10],N_summary_full_sd[:,4:7]))
denit_sd = (np.nansum(denit_sd**2,axis=1))**0.5
nit_N2O_sd = N_summary_sd[:,1]
denit_N2O_sd = (np.nansum(N_summary_full_sd[:,4:7]**2,axis=1))**0.5

fig, ax = plt.subplots(3,1)
ax[0].plot(N_summary[:,0],nit)
ax[0].plot(N_summary[:,0],denit)
ax[0].plot(N_summary[:,0],np.nansum(N_summary_full[:,1:4],axis=1)/10)
ax[0].plot(N_summary[:,0],nit-nit_sd,":",c="C0")
ax[0].plot(N_summary[:,0],nit+nit_sd,":",c="C0")
ax[0].plot(N_summary[:,0],denit-denit_sd,":",c="C1")
ax[0].plot(N_summary[:,0],denit+denit_sd,":",c="C1")
ax[0].legend((["nitrification N (Tg-N a-1)","denitrification N (Tg-N a-1)","total N inputs /10"]))
ax[0].set_ylabel("N flux (Tg-N a-1)")
ax[1].plot(N_summary[:,0],nit/np.nansum(N_summary_full[:,1:4],axis=1)*100)
ax[1].plot(N_summary[:,0],denit/np.nansum(N_summary_full[:,1:4],axis=1)*100)
ax[1].plot(N_summary[:,0],(nit-nit_sd)/np.nansum(N_summary_full[:,1:4],axis=1)*100,":",c="C0")
ax[1].plot(N_summary[:,0],(nit+nit_sd)/np.nansum(N_summary_full[:,1:4],axis=1)*100,":",c="C0")
ax[1].plot(N_summary[:,0],(denit-denit_sd)/np.nansum(N_summary_full[:,1:4],axis=1)*100,":",c="C1")
ax[1].plot(N_summary[:,0],(denit+denit_sd)/np.nansum(N_summary_full[:,1:4],axis=1)*100,":",c="C1")
ax[1].set_ylabel("Proportion of total N inputs")
ax[2].plot(N_summary[:,0],nit_N2O/N_summary[:,1]*100)
ax[2].plot(N_summary[:,0],denit_N2O/N_summary[:,1]*100)
ax[2].plot(N_summary[:,0],(nit_N2O-nit_N2O_sd)/N_summary[:,1]*100,":",c="C0")
ax[2].plot(N_summary[:,0],(nit_N2O+nit_N2O_sd)/N_summary[:,1]*100,":",c="C0")
ax[2].plot(N_summary[:,0],(denit_N2O-denit_N2O_sd)/N_summary[:,1]*100,":",c="C1")
ax[2].plot(N_summary[:,0],(denit_N2O+denit_N2O_sd)/N_summary[:,1]*100,":",c="C1")
ax[2].set_ylabel("Proportion of N2O emissions")

# Product ratio of N2O to N2 from denitrification
RN2O = denit_N2O/(denit_N2O + np.nansum(N_summary_full[:,7:10],axis=1))
a = (denit_N2O_sd**2 + np.nansum(N_summary_full_sd[:,7:10]**2,axis=1)**2)**0.5
a = a/(denit_N2O + np.nansum(N_summary_full[:,7:10],axis=1))
RN2O_sd = ((denit_N2O_sd/denit_N2O)**2 + a**2)**0.5 * RN2O

#%% N gas parameterisations

# most of the code for this is in the parameterisations file!

# final params
para.gp_finalfig(scale_fitNO = (1,post[10],1,1),scale_fitN2 = (1,post[1],1,1))

# final only, with error, to add to first
para.gp_finalfig_error(scale_fitNO = (1,post[10],1,1),scale_fitN2 = (1,post[1],1,1),err_NO = (0,post_sd[10],0,0),err_N2 = (0,post_sd[1],0,0))


#%% Plot growth rates time

def growth_rate_calc(dt,data,period) :
    dt = np.array(dt)
    data = np.array(data)
    res = data.copy()*np.nan
    res_dt = data.copy()*np.nan
    for n in range(0,len(dt)):
        tmp = np.where((dt>=(dt[n]-period)) & (dt<=dt[n]))[0] # identify times from dt-period to dt
        if (len(tmp)>1) & (len(tmp)<5): # diffs between points if few points
            gpy = (data[tmp[1:]] - data[tmp[:-1]])/(dt[tmp[1:]] - dt[tmp[:-1]]) # growth per year
            gpy[np.isinf(gpy)] = np.nan
            res[n] = np.nanmean(gpy)
            res_dt[n] = np.nanmean(dt[tmp])
        if (len(tmp)>=5): # linear fit if many points
            gpy = stats.linregress(dt[tmp],data[tmp])
            if gpy[3] < 0.01: res[n] = gpy[0]
            res_dt[n] = np.nanmean(dt[tmp])
    return(res,res_dt)
           
# find emissions growth rate: dt, all, fert, dep, fix
periods = (3,10)
N2Ogrowth_short = np.zeros((len(N_summary[:,0]),5)) # all, fert, dep, fix
N2Ogrowth_short[:,0] = growth_rate_calc(N_summary[:,0],N_summary[:,0],periods[0])[1]
N2Ogrowth_long = np.zeros((len(N_summary[:,0]),5)) # all, fert, dep, fix
N2Ogrowth_long[:,0] = growth_rate_calc(N_summary[:,0],N_summary[:,0],periods[1])[1]
for i in range(0,4): # all, fert, dep, fix
    N2Ogrowth_short[:,i+1] = growth_rate_calc(N_summary[:,0],N_summary[:,i+1],periods[0])[0]
    N2Ogrowth_long[:,i+1] = growth_rate_calc(N_summary[:,0],N_summary[:,i+1],periods[1])[0]

# plot emissions growth rates
xlimit=((1900,2020))
fig, ax = plt.subplots(4,1)
for n in range(0,4):
    if n==0:
        ax[0].plot(N_summary[:,0],N_summary_full[:,1]+N_summary_full[:,2]+N_summary_full[:,3])
    else: 
        ax[0].plot(N_summary[:,0],N_summary_full[:,n])
    ax[1].plot(N_summary[:,0],N_summary[:,n+1],color="C"+str(n))
    ax[2].plot(N2Ogrowth_short[:,0],N2Ogrowth_short[:,n+1],color="C"+str(n),linestyle="--")
    ax[2].plot(N2Ogrowth_long[:,0],N2Ogrowth_long[:,n+1],color="C"+str(n),linestyle="-")     
ax[0].set_ylabel("N inputs (Tg N a-1)")
ax[0].set_xlim(xlimit)
ax[1].legend((["N2O: all","fert","dep","fix"]))
ax[1].set_ylabel("Emissions (Tg N2O-N a-1)")
ax[1].set_xlim(xlimit)
ax[2].legend(([str(periods[0])+"-yr mean",str(periods[1])+"-yr mean"]))
ax[2].set_ylabel("Emissions growth rate (Tg N2O-N a-1 / a-1)")
ax[2].set_xlim(xlimit)

# calculate also atmos growth rate from model 
N2O_ppb = modres_post[obs[0:147,2]==1].transpose() # no growth rate yet...
atmosgrowth = (N2O_ppb[1:] - N2O_ppb[0:-1])/(tstarts[1:]- tstarts[0:-1])
atmosgrowth_chg = (atmosgrowth[1:] - atmosgrowth[0:-1])/(tstarts[2:]- tstarts[0:-2])
ax[3].plot(tstarts[0:-1]+(tstarts[0:-1]- tstarts[0:-1])/2,atmosgrowth[:])
ax[3].plot(tstarts[0:-2]+(tstarts[2:]- tstarts[0:-2])/2,atmosgrowth_chg[:]*50,color="C0",linestyle="--")
ax[3].legend((["model","growth chg x 50"]))
ax[3].set_ylabel("Atmospheric growth rate")
ax[3].set_xlim(xlimit)


#%% atmos growth rate ### doesn't work yet!

N2O_atmos_all = pd.read_csv("AllAtmosData.csv", sep=',') # get full atmos dataset
# order by date
tmp = np.argsort(N2O_atmos_all["dt"])
N2O_atmos_all = N2O_atmos_all.iloc[tmp]
# calc growth rate
N2Ogrowth_atmos = np.zeros((len(N2O_atmos_all["dt"]),4)) # dt_s, short, dt_l, long
a,b = growth_rate_calc(N2O_atmos_all["dt"],N2O_atmos_all["N2O_ppb"],periods[0])

N2Ogrowth_atmos = np.zeros((len(plot_years),4)) # short, long periods growth; growth rate change
for n in range(0,len(plot_years)) :
    # "short" growth rate
    istart = n-speriod
    if istart<0: istart = 0
    iend = n+speriod
    if iend>=N_summary.shape[0]: iend = N_summary.shape[0]-1
    prev = (N2O_atmos_all["dt"]>plot_years[istart]) & (N2O_atmos_all["dt"]<plot_years[n])
    prev = np.nanmean(N2O_atmos_all["N2O_ppb"][prev])
    now = (N2O_atmos_all["dt"]>plot_years[n]) & (N2O_atmos_all["dt"]<plot_years[iend])
    now = np.nanmean(N2O_atmos_all["N2O_ppb"][now])
    N2Ogrowth_atmos[n,0] = (now-prev)/speriod
    # "long" growth rate
    istart = n-lperiod
    if istart<0: istart = 0
    iend = n+lperiod
    if iend>=N_summary.shape[0]: iend = N_summary.shape[0]-1
    prev = (N2O_atmos_all["dt"]>plot_years[istart]) & (N2O_atmos_all["dt"]<plot_years[n])
    prev = np.nanmean(N2O_atmos_all["N2O_ppb"][prev])
    now = (N2O_atmos_all["dt"]>plot_years[n]) & (N2O_atmos_all["dt"]<plot_years[iend])
    now = np.nanmean(N2O_atmos_all["N2O_ppb"][now])
    N2Ogrowth_atmos[n,1] = (now-prev)/lperiod
N2Ogrowth_atmos[1:,2] = (N2Ogrowth_atmos[1:,0] - N2Ogrowth_atmos[0:-1,0])/(plot_years[1:]-plot_years[0:-1])
N2Ogrowth_atmos[1:,3] = (N2Ogrowth_atmos[1:,1] - N2Ogrowth_atmos[0:-1,1])/(plot_years[1:]-plot_years[0:-1])
# calculate also atmos growth rate from model 
N2O_ppb = np.array((modres_post[0:(int(len(obs)/3))],obs[0:(int(len(obs)/3)),0])).transpose()
atmosgrowth = (N2O_ppb[1:,:] - N2O_ppb[0:-1,:])/np.array((tstarts[1:-1]- tstarts[0:-2],tstarts[1:-1]- tstarts[0:-2])).transpose()
atmosgrowth_chg = (atmosgrowth[1:,:] - atmosgrowth[0:-1,:])/np.array((tstarts[2:-1]- tstarts[0:-3],tstarts[2:-1]- tstarts[0:-3])).transpose()

#%% emission factors with time for leach, NH3, N2, NO, N2O for fert, dep, fix, all

Nsumm_fertscale = N_summary_full.copy()
Nsumm_fertscale[:,1] = Nsumm_fertscale[:,1]*post[3] # scale fert inputs to already remove plant uptake, so that EFs are comparable
Nsumm_fertscale_sd = N_summary_full_sd.copy()
Nsumm_fertscale_sd[:,1] = Nsumm_fertscale_sd[:,1]*post[3]
xlimit=((1940,2020))
fig, ax = plt.subplots(5,1)
for n in range(0,4):
    # leaching
    if n == 0: # all inputs
        tmp = np.nansum(Nsumm_fertscale[:,16:19],axis=1)/np.nansum(Nsumm_fertscale[:,1:4],axis=1)*100
        tmp_sd = np.nansum(Nsumm_fertscale_sd[:,16:19]**2,axis=1)**0.5/np.nansum(Nsumm_fertscale[:,16:19],axis=1)*tmp
    if n > 0: # fert, dep, fix
        tmp = Nsumm_fertscale[:,15+n]/Nsumm_fertscale[:,n]*100
        tmp_sd = Nsumm_fertscale_sd[:,15+n]/Nsumm_fertscale[:,15+n]*tmp
        if n==1: tmp_sd = Nsumm_fertscale_sd[:,15+n]/Nsumm_fertscale[:,15+n]*tmp/10
    ax[0].plot(N_summary[:,0],tmp,c="C"+str(n))
    ax[0].plot(N_summary[:,0],tmp-tmp_sd,":",c="C"+str(n)); ax[0].plot(N_summary[:,0],tmp+tmp_sd,":",c="C"+str(n))
    ax[0].set_ylabel("Leaching: % of N")
    ax[0].set_xlim(xlimit)
    # NH3
    if n == 0: # all inputs
        tmp = np.nansum(Nsumm_fertscale[:,13:16],axis=1)/np.nansum(Nsumm_fertscale[:,1:4],axis=1)*100
        tmp_sd = np.nansum(Nsumm_fertscale_sd[:,13:16]**2,axis=1)**0.5/np.nansum(Nsumm_fertscale[:,13:16],axis=1)*tmp
    if n > 0: # fert, dep, fix
        tmp = Nsumm_fertscale[:,12+n]/Nsumm_fertscale[:,n]*100
        tmp_sd = Nsumm_fertscale_sd[:,12+n]/Nsumm_fertscale[:,12+n]*tmp
        if n==1: tmp_sd = Nsumm_fertscale_sd[:,12+n]/Nsumm_fertscale[:,12+n]*tmp/4
    ax[1].plot(N_summary[:,0],tmp,c="C"+str(n))
    ax[1].plot(N_summary[:,0],tmp-tmp_sd,":",c="C"+str(n)); ax[1].plot(N_summary[:,0],tmp+tmp_sd,":",c="C"+str(n))
    ax[1].set_ylabel("NH3: % of N")
    ax[1].set_xlim(xlimit)
    # N2
    if n == 0: # all inputs
        tmp = np.nansum(Nsumm_fertscale[:,7:10],axis=1)/np.nansum(Nsumm_fertscale[:,1:4],axis=1)*100
        tmp_sd = np.nansum(Nsumm_fertscale_sd[:,7:10]**2,axis=1)**0.5/np.nansum(Nsumm_fertscale[:,7:10],axis=1)*tmp
    if n > 0: # fert, dep, fix
        tmp = Nsumm_fertscale[:,6+n]/Nsumm_fertscale[:,n]*100
        tmp_sd = Nsumm_fertscale_sd[:,6+n]/Nsumm_fertscale[:,6+n]*tmp
        if n==1: tmp_sd = Nsumm_fertscale_sd[:,6+n]/Nsumm_fertscale[:,6+n]*tmp/2
    ax[2].plot(N_summary[:,0],tmp,c="C"+str(n))
    ax[2].plot(N_summary[:,0],tmp-tmp_sd,":",c="C"+str(n)); ax[2].plot(N_summary[:,0],tmp+tmp_sd,":",c="C"+str(n))
    ax[2].set_ylabel("NH3: % of N")
    ax[2].set_xlim(xlimit)
    # NO
    if n == 0: # all inputs
        tmp = np.nansum(Nsumm_fertscale[:,10:13],axis=1)/np.nansum(Nsumm_fertscale[:,1:4],axis=1)*100
        tmp_sd = np.nansum(Nsumm_fertscale_sd[:,10:13]**2,axis=1)**0.5/np.nansum(Nsumm_fertscale[:,10:13],axis=1)*tmp
    if n > 0: # fert, dep, fix
        tmp = Nsumm_fertscale[:,9+n]/Nsumm_fertscale[:,n]*100
        tmp_sd = Nsumm_fertscale_sd[:,9+n]/Nsumm_fertscale[:,9+n]*tmp
        if n==1: tmp_sd = Nsumm_fertscale_sd[:,9+n]/Nsumm_fertscale[:,9+n]*tmp/3
    ax[3].plot(N_summary[:,0],tmp,c="C"+str(n))
    ax[3].plot(N_summary[:,0],tmp-tmp_sd,":",c="C"+str(n)); ax[3].plot(N_summary[:,0],tmp+tmp_sd,":",c="C"+str(n))
    ax[3].set_ylabel("NO: % of N")
    ax[3].set_xlim(xlimit)
    # N2O
    if n == 0: # all inputs
        tmp = np.nansum(N_summary[:,2:5],axis=1)/np.nansum(Nsumm_fertscale[:,1:4],axis=1)*100
        tmp_sd = np.nansum(N_summary_sd[:,2:5]**2,axis=1)**0.5/np.nansum(N_summary[:,2:5],axis=1)*tmp
    if n > 0: # fert, dep, fix
        tmp = N_summary[:,1+n]/Nsumm_fertscale[:,n]*100
        tmp_sd = N_summary_sd[:,1+n]/N_summary[:,1+n]*tmp
        if n==1: tmp_sd = N_summary_sd[:,1+n]/N_summary[:,1+n]*tmp/3
    ax[4].plot(N_summary[:,0],tmp,c="C"+str(n))
    ax[4].plot(N_summary[:,0],tmp-tmp_sd,":",c="C"+str(n)); ax[4].plot(N_summary[:,0],tmp+tmp_sd,":",c="C"+str(n))
    ax[4].set_ylabel("N2O: % of N")
    ax[4].set_xlim(xlimit)
ax[0].legend((["all","fert","dep","fix"]))

# the same just for fert 
xlimit=((1960,2020))
fig, ax = plt.subplots(2,1)
sty = ("-","--",":")
for n in range(1,4):
    if n==1:
        leach = N_summary_full[:,15+n]/N_summary_full[:,n]*100
        NH3 = N_summary_full[:,12+n]/N_summary_full[:,n]*100
        N2 = N_summary_full[:,6+n]/N_summary_full[:,n]*100
        NO = N_summary_full[:,9+n]/N_summary_full[:,n]*100
        N2O = N_summary[:,1+n]/N_summary_full[:,n]*100
    else:
        leach = N_summary_full[:,15+n]/N_summary_full[:,n]*100
    ax[0].plot(N_summary[:,0],leach,sty[n-1],color="C0")
    ax[1].plot(N_summary[:,0],NH3,sty[n-1],color="C1")
    ax[1].plot(N_summary[:,0],N2,sty[n-1],color="C2")
    ax[1].plot(N_summary[:,0],NO,sty[n-1],color="C3")
    ax[1].plot(N_summary[:,0],N2O,sty[n-1],color="C4")
    ax[0].set_ylabel("Leaching: % of N")
    ax[0].set_xlim(xlimit)
    ax[1].set_ylabel("Other: % of N")
    ax[1].set_xlim(xlimit)
    ax[0].legend((["fert","dep","fix"]))
    ax[1].legend((["NH3","N2","NO","N2O"]))

#%% emission factors with time for leach, NH3, N2, NO, N2O for fert, dep, fix, all (anthrop only)

xlimit=((1960,2020))
x = (N_summary[:,0]>=xlimit[0]) & (N_summary[:,0]<=xlimit[1])
fig, ax = plt.subplots(5,1)
for n in range(0,4):
    # leaching
    if n == 0: # all inputs
        tmp = np.nansum(N_summary_full[:,16:19]-N_summary_full[0,16:19],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax[0].plot(N_summary[x,0],tmp[x]*100)
    if n > 0: # fert, dep, fix
        tmp = (N_summary_full[:,15+n]-N_summary_full[0,15+n])/(N_summary_full[:,n]-N_summary_full[0,n])
        ax[0].plot(N_summary[x,0],tmp[x]*100)
    ax[0].set_ylabel("Leaching: % of N")
    ax[0].set_xlim(xlimit)
    # NH3
    if n == 0: # all inputs
        tmp = np.nansum(N_summary_full[:,13:16]-N_summary_full[0,13:16],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax[1].plot(N_summary[x,0],tmp[x]*100)
    if n > 0: # fert, dep, fix
        tmp = (N_summary_full[:,12+n]-N_summary_full[0,12+n])/(N_summary_full[:,n]-N_summary_full[0,n])
        ax[1].plot(N_summary[x,0],tmp[x]*100)
    ax[1].set_ylabel("NH3: % of N")
    ax[1].set_xlim(xlimit)
    # N2
    if n == 0: # all inputs
        tmp = np.nansum(N_summary_full[:,7:10]-N_summary_full[0,7:10],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax[2].plot(N_summary[x,0],tmp[x]*100)
    if n > 0: # fert, dep, fix
        tmp = (N_summary_full[:,6+n]-N_summary_full[0,6+n])/(N_summary_full[:,n]-N_summary_full[0,n])
        ax[2].plot(N_summary[x,0],tmp[x]*100)
    ax[2].set_ylabel("N2: % of N")
    ax[2].set_xlim(xlimit)
    # NO
    if n == 0: # all inputs
        tmp = np.nansum(N_summary_full[:,10:13]-N_summary_full[0,10:13],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax[3].plot(N_summary[x,0],tmp[x]*100)
    if n > 0: # fert, dep, fix
        tmp = (N_summary_full[:,9+n]-N_summary_full[0,9+n])/(N_summary_full[:,n]-N_summary_full[0,n])
        ax[3].plot(N_summary[x,0],tmp[x]*100)
    ax[3].set_ylabel("NO: % of N")
    ax[3].set_xlim(xlimit)
    # N2O
    if n == 0: # all inputs
        tmp = np.nansum(N_summary[:,2:5]-N_summary[0,2:5],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax[4].plot(N_summary[x,0],tmp[x]*100)
    if n > 0: # fert, dep, fix
        tmp = (N_summary[:,1+n]-N_summary[0,1+n])/(N_summary_full[:,n]-N_summary_full[0,n])
        ax[4].plot(N_summary[x,0],tmp[x]*100)
    ax[4].set_ylabel("N2O: % of N")
    ax[4].set_xlim(xlimit)
ax[0].legend((["all","fert","dep","fix"]))

# N2O EFs vs. inputs (can spatial changes account for Rona's idea that EFs are higher at higher input levels?)
fig, ax = plt.subplots(1,1)
x = (N_summary[:,0]>=1900) & (N_summary[:,0]<=2020)
for n in range(0,4):
    if n == 0: # all inputs
        # anthrop only
        tmp = np.nansum(N_summary[:,2:5]-N_summary[0,2:5],axis=1)/np.nansum(N_summary_full[:,1:4]-N_summary_full[0,1:4],axis=1)
        ax.plot(np.nansum(N_summary_full[x,1:4]-N_summary_full[0,1:4],axis=1),tmp[x]*100,"o",markeredgecolor="k")
        # all
        ax.plot(np.nansum(N_summary_full[x,1:4],axis=1),np.nansum(N_summary[x,2:5],axis=1)/np.nansum(N_summary_full[x,1:4],axis=1)*100,"o",color="C"+str(n))
    if n > 0: # fert, dep, fix
        tmp = (N_summary[:,1+n]-N_summary[0,1+n])/(N_summary_full[:,n]-N_summary_full[0,n]) # anthrop only
        ax.plot((N_summary_full[x,n]-N_summary_full[0,n]),tmp[x]*100,"o",color="C"+str(n),markeredgecolor="k")
        ax.plot(N_summary_full[x,n],N_summary[x,1+n]/N_summary_full[x,n]*100,"o",color="C"+str(n)) # all
    ax.set_ylim((0,10))

#%% Inputs and N2O vs. various params

# conclusion: f_G and f_N2O increase for fert over the years, causing the EF to increase for fert. This is because of a movement into 
# 100 lon 25 lat bin eg. China, which has a greater proportion of emissions further south than Europe and US.
# combine this plot iwth a plot of k_G*f_N2O for the years shown...

param = k_G_post*f_N2O; param_name = "fraction of N as N2O gas"
#param = d15N_grid; param_name = "d15N"
param_rng = np.where((~np.isnan(param)))# & (WFPS_grid!=100)) # don't use the ==100 values if plotting cumsum figure
param_vec = param[param_rng]
factor_area = (area_grid/1000/1000/1000/1000)[param_rng] # convert to Tg-N y-1
# get input fields for diff years to compare
years_n = (1940,1980,2020)
results = np.zeros((len(param_rng[0]),4,len(years_n))) # index, (all, fert, dep, fix), year
results_cumu = results.copy()
resultsA = results.copy()
results_cumuA = results.copy()
index0 = int(np.where(datarange==1850)[0])
for n in range(0,len(years_n)):
    index = int(np.where(datarange==years_n[n])[0])
    fix = Ninputs.variables["fixation"][index,:,:].data[param_rng]*factor_area
    dep = Ninputs.variables["deposition"][index,:,:].data[param_rng]*factor_area
    fert = Ninputs.variables["fertilisation"][index,:,:].data[param_rng]*factor_area
    fixA = fix - Ninputs.variables["fixation"][index0,:,:].data[param_rng]*factor_area
    depA = dep - Ninputs.variables["deposition"][index0,:,:].data[param_rng]*factor_area
    fertA = fert - Ninputs.variables["fertilisation"][index0,:,:].data[param_rng]*factor_area
    allin = np.nansum((fix,dep,fert),axis=0)
    allinA = np.nansum((fixA,depA,fertA),axis=0)
    results[:,:,n] = np.vstack((allin,fert,dep,fix)).transpose()
    resultsA[:,:,n] = np.vstack((allinA,fertA,depA,fixA)).transpose()
    # order by the param and find cumulative sum
    tmp = param_vec.argsort()
    results_cumu[:,:,n] = np.cumsum(results[tmp,:,n],axis=0)
    results_cumuA[:,:,n] = np.cumsum(resultsA[tmp,:,n],axis=0)

# amount and % of N per bin
bins = 20
binmin = np.nanmin(param)
binmax = np.nanmax(param)
binhigh = np.nanmean(param) + 3*np.nanstd(param) # using binhigh is to avoid long tail at high values
binsize = (binhigh - binmin)/bins
binvalues = np.arange(binmin,binhigh,binsize)
N_perbin = np.zeros((bins,5,len(years_n))) # bins, (param bin mid, all, fert, dep, fix), years
relN_perbin = N_perbin.copy()
N_perbinA = N_perbin.copy()
relN_perbinA = N_perbin.copy()
r = 0
titles = ("all","fert","dep","fix")
for i in binvalues:
    rng = np.where((param_vec>=i) & (param_vec<(i+binsize))) # find data in range
    if i==binvalues[-1]: rng = np.where(param_vec>=i)
    N_perbin[r,0,:] = i+binsize/2 # bin mid
    relN_perbin[r,0,:] = i+binsize/2
    N_perbinA[r,0,:] = i+binsize/2 # bin mid
    relN_perbinA[r,0,:] = i+binsize/2
    for n in range(0,len(years_n)): # for each year
        N_perbin[r,1:5,n] = np.nansum(results[rng,:,n],axis=1)
        relN_perbin[r,1:5,n] = np.nansum(results[rng,:,n],axis=1)/np.nansum(results[:,:,n],axis=0)*100
        N_perbinA[r,1:5,n] = np.nansum(resultsA[rng,:,n],axis=1)
        relN_perbinA[r,1:5,n] = np.nansum(resultsA[rng,:,n],axis=1)/np.nansum(resultsA[:,:,n],axis=0)*100
    r = r+1
fig, ax = plt.subplots(4,2)    
styles = (":","--","-","-.")
for n in range(0,4): # loop through all fert dep fix
    for i in range(0,len(years_n)): # loop through years
        # LHS = absolute N per bin
        ax[n,0].plot(N_perbin[:,0,i],N_perbin[:,n+1,i],linestyle = styles[i],color = "C"+str(n),alpha=0.2)
        ax[n,0].plot(N_perbinA[:,0,i],N_perbinA[:,n+1,i],linestyle = styles[i],color = "C"+str(n))
        # RHS = relative N per bin
        ax[n,1].plot(N_perbin[:,0,i],relN_perbin[:,n+1,i],linestyle = styles[i],color = "C"+str(n),alpha=0.2)
        ax[n,1].plot(N_perbinA[:,0,i],relN_perbinA[:,n+1,i],linestyle = styles[i],color = "C"+str(n))
    ax[n,0].set_ylabel("N in bin ("+titles[n]+", Tg-N y-1)")
    ax[n,1].set_ylabel("%N in bin ("+titles[n]+")")
    ax[n,0].set_xlabel(param_name)
ax[0,0].legend((years_n))

#%% calc absolute difference in case emissions had stayed the same and also without temp sens

t1=1940; t2=2020

# get T sens
tindex = int(np.where(abs(Tanom.variables["time"][:]-t1) == np.min(abs(Tanom.variables["time"][:]-t1)))[0])
tanomaly_sens_t1 = (post[4]-1)*Tnorm[tindex,:,:] + 1
tindex = int(np.where(abs(Tanom.variables["time"][:]-t2) == np.min(abs(Tanom.variables["time"][:]-t2)))[0])
tanomaly_sens_t2 = (post[4]-1)*Tnorm[tindex,:,:] + 1

# change due to  changing distribution
i_t2 = int(np.where(datarange==t2)[0])
allin_t2 = Ninputs.variables["fertilisation"][i_t2,:,:]*post[3]# + Ninputs.variables["fixation"][i_t2,:,:] + Ninputs.variables["deposition"][i_t2,:,:]
i_t1 = int(np.where(datarange==t1)[0])
allin_t1 = Ninputs.variables["fertilisation"][i_t1,:,:]*post[3]# + Ninputs.variables["fixation"][i_t1,:,:] + Ninputs.variables["deposition"][i_t1,:,:]
allin_t2_t1d = allin_t1*np.nansum(allin_t2[param_rng]*factor_area)/np.nansum(allin_t1[param_rng]*factor_area)
t2_emissions = k_G_post*f_N2O_post*allin_t2*tanomaly_sens_t2
t2_emissions_t1d = k_G_post*f_N2O_post*allin_t2_t1d*tanomaly_sens_t2
relchange = t2_emissions - t2_emissions_t1d
print("t2:",str(np.nansum(t2_emissions[param_rng]*factor_area)),", t2_t1d:",str(np.nansum(t2_emissions_t1d[param_rng]*factor_area)),", change:",str(np.nansum(relchange[param_rng]*factor_area)))
plot_map(LON,LAT,relchange.clip(-0.15,0.15),"Change in fert N2O emissions since 1940 due to spatial distribution only",cmap="bwr")

# change only due to temp_sens
allin_t2 = Ninputs.variables["fertilisation"][i_t2,:,:]*post[3] + Ninputs.variables["fixation"][i_t2,:,:] + Ninputs.variables["deposition"][i_t2,:,:]
t2_emissions = k_G_post*f_N2O_post*allin_t2*tanomaly_sens_t2
t2_emissions_notempsens = k_G_post*f_N2O_post*allin_t2*tanomaly_sens_t1
relchange = t2_emissions - t2_emissions_notempsens
print("t2:",str(np.nansum(t2_emissions[param_rng]*factor_area)),", t2_t1d:",str(np.nansum(t2_emissions_notempsens[param_rng]*factor_area)),", change:",str(np.nansum(relchange[param_rng]*factor_area)))
plot_map(LON,LAT,relchange.clip(-0.15,0.15),"Change in fert N2O emissions since 1940 due to temp_sens only",cmap="bwr",vminmax=(-0.15,0.15)) 

#%% "Mask" and values for different geographical regions

# different regions, lon start end, lat start end
europe = np.array((-5,20,40,55))
northamerica = np.array((-125,-75,30,50))
southamerica = np.array((-75,-35,-30,5))
china = np.array((88,120,22,42))
india = np.array((70,88,10,30))
southasia = np.array((95,150,-10,22))
australia = np.array((115,155,-40,-15))
subsafrica = np.array((-15,50,-35,15))
regions = np.vstack((europe,northamerica,southamerica,china,india,southasia,australia,subsafrica))
regionnames = ("europe","northamerica","southamerica","china","india","southasia","australia","subsafrica")

# plot the regions
plot_map(LON,LAT,k_G_post*f_N2O_post*100,"fraction of Nr lost to N2O, post")
for n in range(0,len(regionnames)):
    plt.plot((regions[n,0],regions[n,0]),(regions[n,2],regions[n,3]),c="C"+str(n))
    plt.plot((regions[n,1],regions[n,1]),(regions[n,2],regions[n,3]),c="C"+str(n))
    plt.plot((regions[n,0],regions[n,1]),(regions[n,2],regions[n,2]),c="C"+str(n))
    plt.plot((regions[n,0],regions[n,1]),(regions[n,3],regions[n,3]),c="C"+str(n))
    plt.annotate(regionnames[n], (regions[n,1]+0.01, regions[n,3]+0.01))
    
# find modelled mean values of N2O total, fN2O, f_G, d15N grid, d15N actual, and stdev/n of each
years_n = np.array((1990,2020))
datarng = np.where((~np.isnan(preamble.d15N_grid)))# & (WFPS_grid!=100)) # don't use the ==100 values if plotting cumsum figure
factor_area = (area_grid/1000/1000/1000/1000)[datarng] # convert to Tg-N y-1
index = np.zeros(len(years_n))
for n in range(0,len(years_n)):  
    index[n] = int(np.where(datarange==years_n[n])[0]) # choose year
lon_vec = LON[datarng]
lat_vec = LAT[datarng]
for n in range(0,len(regionnames)):  
    bounding = np.where((lon_vec>=regions[n,0]) & (lon_vec<=regions[n,1]) & (lat_vec>=regions[n,2]) & (lat_vec<=regions[n,3])) # select the region  
    # total fert N input and N2O emitted for each region
    fertinput = np.zeros(len(years_n))
    N2Oout = np.zeros(len(years_n))
    for i in range(0,len(years_n)):  
        tmp = Ninputs.variables["fertilisation"][index[i],:,:].data[datarng]*factor_area # fert in Tg y-1
        fertinput[i] = np.nansum(tmp[bounding])
        tmp = tmp[bounding]*k_G_post[datarng][bounding]*f_N2O_post[datarng][bounding]*post[3]
        N2Oout[i] = np.nansum(tmp)
    # k_G for each region (area weighted)
    k_G_region = np.nansum((k_G_post[datarng][bounding]*factor_area[bounding])/np.nansum(factor_area[bounding]))
    # EF fert for each region
    fN2Ofert_region = np.nansum((k_G_post[datarng][bounding]*f_N2O_post[datarng][bounding]*post[3]*factor_area[bounding])/np.nansum(factor_area[bounding]))
    # d15N of all measurements in the region
    ### need to get the d15N data from the d15N script before running this!!
    y = 0
    if y==1: 
        bounding = np.where((soil15N[:,1]>=regions[n,0]) & (soil15N[:,1]<=regions[n,1]) & (soil15N[:,0]>=regions[n,2]) & (soil15N[:,0]<=regions[n,3])) # select the region  
        d15N = np.array((np.nanmean(soil15N[bounding,5]),np.nanstd(soil15N[bounding,5]),len(bounding[0])))
    else : d15N = np.array((1,1,1))
    # EF of all measurements in the region
    bounding = np.where((N2O_fluxes_short[:,7]>=regions[n,0]) & (N2O_fluxes_short[:,7]<=regions[n,1]) & (N2O_fluxes_short[:,6]>=regions[n,2]) & (N2O_fluxes_short[:,6]<=regions[n,3])) # select the region  
    EFregions = np.array((np.nanmean(N2O_fluxes_short[bounding,10]),np.nanstd(N2O_fluxes_short[bounding,10]),len(bounding[0])))
    # combine and output
    tmp = np.hstack((fertinput,N2Oout,k_G_region,fN2Ofert_region,d15N,EFregions))
    if n == 0: 
        output = tmp
    else: output = vstack((tmp,output))
# get the names
names = []
for i in range(0,len(years_n)): names.append("fert"+str(years_n[i]))
for i in range(0,len(years_n)): names.append("N2O"+str(years_n[i]))
names.append("k_G"); names.append("EF_fert"); names.append("d15N"); names.append("d15N_sd"); names.append("d15N_n");
names.append("EF"); names.append("EF_sd"); names.append("EF_n");
df = pd.DataFrame(output,columns=names,index=regionnames)
df.iloc[-1]

#%% Get the data exported for Stefan

boxes = np.array((-90,-30,0,30,90))
# find modelled annual emissions per box per year
years_n = np.arange(1850,2020)
datarng = np.where((~np.isnan(preamble.d15N_grid)))
factor_area = (area_grid/1000/1000/1000/1000)[datarng] # convert to Tg-N y-1
lon_vec = LON[datarng]
lat_vec = LAT[datarng]
results = np.zeros((len(years_n),7))
for n in range(0,len(years_n)):  
    tindex = int(np.where(abs(Tanom.variables["time"][:]-years_n[n]) == np.min(abs(Tanom.variables["time"][:]-years_n[n])))[0])
    tanomaly = Tnorm[tindex,:,:]
    tanomaly_sens = (post[4]-1)*tanomaly + 1
    index = np.where(Ninputs.variables["time"][:].data == years_n[n])[0] # select the year
    fix = ( Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens )[0,:,:]
    dep = ( Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens  )[0,:,:]
    fert = ( Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] )[0,:,:]
    TRO = TRO_grid[index,:,:][0,:,:]
    WWT = WWT_grid[index,:,:][0,:,:]
    CHE = CHE_grid[index,:,:][0,:,:]
    ENE = ENE_grid[index,:,:][0,:,:]
    N2O_per_grid = ((fix+dep+fert)*k_G*f_N2O_post + TRO + WWT + ENE + CHE)[datarng]*factor_area
    N2O_per_grid_nonsoil = (TRO + WWT + ENE + CHE)[datarng]*factor_area
    results[n,0] = years_n[n]
    results[n,1] = np.nansum(N2O_per_grid)
    results[n,2] = np.nansum(N2O_per_grid_nonsoil)
    for i in range(0,4):
        bounding = np.where((lat_vec>=boxes[i]) & (lat_vec<boxes[i+1])) # select the region 
        results[n,i+3] = np.nansum(N2O_per_grid[bounding])
df = pd.DataFrame(results,columns=("year","N2O_all_TgN2O-2_a-1","N2O_nonsoil_TgN2O-2_a-1","N2O_all_90S-30S_TgN2O-2_a-1","N2O_all_30S-0_TgN2O-2_a-1","N2O_all_0-30N_TgN2O-2_a-1","N2O_all_30N-90N_TgN2O-2_a-1"))
df.to_csv("StefanHenne_N2O_boxestimates.txt")

#%% Get the data exported forKenya

ssa = np.array((-15,50,-35,15))
# find modelled annual emissions per box per year
years_n = np.arange(1850,2049)
datarng = np.where((~np.isnan(preamble.d15N_grid)))# & (WFPS_grid!=100)) # don't use the ==100 values if plotting cumsum figure
factor_area = (area_grid/1000/1000/1000/1000)[datarng] # convert to Tg-N y-1
lon_vec = LON[datarng]
lat_vec = LAT[datarng]
results = np.zeros((len(years_n),7))
baseyear = 1995; finalyear = 2045 # Years between which to plot a gridded comparison
for n in np.arange(0,len(years_n)):  
    tindex = int(np.where(abs(Tanom.variables["time"][:]-years_n[n]) == np.min(abs(Tanom.variables["time"][:]-years_n[n])))[0])
    tanomaly = Tnorm[tindex,:,:]
    tanomaly_sens = (post[4]-1)*tanomaly + 1
    index = np.where(Ninputs.variables["time"][:].data == years_n[n])[0] # select the year
    fix = ( Ninputs.variables["fixation"][index,:,:].data *tanomaly_sens )[0,:,:]
    dep = ( Ninputs.variables["deposition"][index,:,:].data *tanomaly_sens  )[0,:,:]
    fert = ( Ninputs.variables["fertilisation"][index,:,:].data *tanomaly_sens *post[3] )[0,:,:]
    TRO = TRO_grid[index,:,:][0,:,:]
    WWT = WWT_grid[index,:,:][0,:,:]
    CHE = CHE_grid[index,:,:][0,:,:]
    ENE = ENE_grid[index,:,:][0,:,:]
    N2O_per_grid = ((fix+dep+fert)*k_G*f_N2O_post + TRO + WWT + ENE + CHE)[datarng]*factor_area
    results[n,0] = years_n[n]
    results[n,1] = np.nansum(N2O_per_grid)
    bounding = np.where((lon_vec>=ssa[0]) & (lon_vec<=ssa[1]) & (lat_vec>=ssa[2]) & (lat_vec<=ssa[3])) # select the region  
    results[n,2] = np.nansum(N2O_per_grid[bounding])
    # Look at fert only
    N2O_per_grid_fert = (fert*k_G*f_N2O_post)[datarng]*factor_area
    results[n,5] = np.nansum(N2O_per_grid_fert)
    results[n,6] = np.nansum(N2O_per_grid_fert[bounding])
    # Gridded comparison
    plotbounding = ((lon_vec>=-25) & (lon_vec<=60) & (lat_vec>=-40) & (lat_vec<=35)) # select the region 
    compgrid = d15N_grid*0 + np.nan; baseplot = d15N_grid*0 + np.nan; finalplot = d15N_grid*0 + np.nan; 
    if years_n[n] == baseyear: 
        basegrid = (fert*k_G*f_N2O_post)[datarng]
        basegrid[~plotbounding] = np.nan
        baseplot[datarng] = basegrid
        plot_map(LON,LAT,baseplot.clip(0,0.05),"N2O from fertilisation, 1995",cmap="viridis")
        plt.scatter(35+17/60,31/60,s=50, c="yellow",marker=(5, 1))
    if years_n[n] == finalyear: 
        finalgrid = (fert*k_G*f_N2O_post)[datarng]
        finalgrid[~plotbounding] = np.nan
        finalplot[datarng] = finalgrid
        plot_map(LON,LAT,finalplot.clip(0,0.05),"N2O from fertilisation, 2045",cmap="viridis")
        plt.scatter(35+17/60,31/60,s=50, c="yellow",marker=(5, 1))
        compdata = (finalgrid-basegrid)/basegrid*100
        compdata[basegrid<=0] = np.nan
        compgrid[datarng] = compdata
        plot_map(LON,LAT,compgrid.clip(0,5000),"Predicted % increase in N2O from fertilisation, 1995-2045",cmap="viridis")
        plt.scatter(35+17/60,31/60,s=50, c="yellow",marker=(5, 1))
results[:,3] = results[:,1] - results[0,1]
results[:,4] = results[:,2] - results[0,2]
df = pd.DataFrame(results,columns=("year","N2O_all","N2O_ssa","N2O_anthall","N2O_anthssa","N2O_fertall","N2O_fertssa"))

fig, ax = plt.subplots(2)   
ax[0].plot(df["year"],df["N2O_anthssa"])
ax[0].plot(df["year"],df["N2O_ssa"])
ax[0].set_xlim((2000,2040))
ax[0].set_ylim((4,10))
