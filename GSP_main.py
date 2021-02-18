#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@version: 1.4.0
@file: GSP_main.py
@time: 2021/1/26 10:50
@functions: graph signal processing main script
@update: support Yeo-ICN definition
@update: support ICN-level brain activity and connecitivty strength saving
"""

import numpy as np
import glob
import os
import time
import matplotlib.pyplot as plt
from pygsp import graphs, filters, plotting
from GSP_utilities import surrogate_BOLD_create, save_variable, load_variable
import pandas as pd
from dppd import dppd
dp, X = dppd()

# 1. path locations and parameters
start = time.time()
deriv_path = '/home/amax/data/cye/MScohort_BIDS_clean/derivatives'
connectome_path = os.path.join(deriv_path, 'mrtrix')
xcpengine_path = os.path.join(deriv_path, 'xcpengine')
network_assign_path = 'CAB-NP_v1.1_Labels-ReorderedbyNetworks_Yeo.csv' 
num_BOLD_timepoints = 180
num_rand = 100 # number of surrogates
functional_type = 'BOLD'
tract_type = 'meanlength' # one of the following: invlength, invnodevol, level-participant_connectome, meanlength
ICN_type = 'Yeo' # one of the following: 'Yeo', 'Cole'
normalize_type = 'both' # 'W': normalize W; 'L': normalize Laplacian (Preti method); 'both': normalize both W and Laplacian 


# 2. read network assignment for hcpmmp
network_assign_csv = pd.read_csv(network_assign_path)

network_assign_csv = dp(network_assign_csv).mutate(NETWORK=X.Yeo_NETWORK).pd
network_assign_csv = dp(network_assign_csv).mutate(NETWORKKEY=X.Yeo_NETWORKKEY).pd


num_network_df = dp(network_assign_csv).summarise((X.NETWORKKEY, np.max, 'hp_max')).pd
num_network = num_network_df.iloc[0,0]
network_rowindex_ls = []
for network_i in range(1,num_network+1):
    df_network = dp(network_assign_csv).filter_by(X.NETWORKKEY == network_i).pd
    network_rowindex_ls.append(df_network.index.values)
    network_unique_df = dp(network_assign_csv).distinct('NETWORKKEY').pd
    network_unique_df = network_unique_df.sort_values(by='NETWORKKEY',ascending = True)

network_unique_df = dp(network_unique_df).filter_by(-X.NETWORK.isin(['Undefine'])).pd # remove undefined ICN
network_unique_df = network_unique_df.reset_index()

# 3. define group of interests
cohort1 = 'ms'
cohort2 = 'nc'
cohort3 = 'nmo'
cohort4 = 'cis'

cohort1_connectome_ls = glob.glob(os.path.join(connectome_path, 'sub-' + cohort1 + '*'))
cohort2_connectome_ls = glob.glob(os.path.join(connectome_path, 'sub-' + cohort2 + '*'))
cohort3_connectome_ls = glob.glob(os.path.join(connectome_path, 'sub-' + cohort3 + '*'))
cohort4_connectome_ls = glob.glob(os.path.join(connectome_path, 'sub-' + cohort4 + '*'))
cohort_connectome_ls = cohort1_connectome_ls + cohort2_connectome_ls + cohort3_connectome_ls + cohort4_connectome_ls
cohort_connectome_ls.sort()

cohort1_fmri_ls = glob.glob(os.path.join(xcpengine_path, 'sub-' + cohort1 + '*'))
cohort2_fmri_ls = glob.glob(os.path.join(xcpengine_path, 'sub-' + cohort2 + '*'))
cohort3_fmri_ls = glob.glob(os.path.join(xcpengine_path, 'sub-' + cohort3 + '*'))
cohort4_fmri_ls = glob.glob(os.path.join(xcpengine_path, 'sub-' + cohort4 + '*'))
cohort_fmri_ls = cohort1_fmri_ls + cohort2_fmri_ls + cohort3_fmri_ls + cohort4_fmri_ls


cohort_name_ls = [os.path.basename(item) for item in cohort_connectome_ls]
remove_name_ls = ['sub-nc011','sub-nc039', 'sub-nmo002', 'sub-nmo019', 'sub-cis002','sub-cis015', 'sub-ms015'] # problematic cases
cohort_name_ls = list(set(cohort_name_ls) - set(remove_name_ls))  # remove problematic cases

for i in remove_name_ls: # remove problematic cases
    cohort_connectome_ls = [x for x in cohort_connectome_ls if i not in x]
    cohort_fmri_ls = [x for x in cohort_fmri_ls if i not in x]

cohort_name_ls.sort()
cohort_connectome_ls.sort()
cohort_fmri_ls.sort()

if len(cohort_connectome_ls) != len(cohort_fmri_ls):
    print('Number of connectome and xcpengine results not matched')

# 4. create a dataframe to store individual filepath 
path_dict = {'subname':cohort_name_ls, 'mrtrix_path': cohort_connectome_ls, 'xcp_path':cohort_fmri_ls}
path_df = pd.DataFrame(path_dict, columns=['subname','mrtrix_path','xcp_path'])
path_df = dp(path_df).mutate(connectome_path=X.mrtrix_path + '/connectome/' + X.subname +'_parc-hcpmmp1_'  + tract_type + '.csv').pd
path_df = dp(path_df).mutate(BOLD_series_path=X.xcp_path + '/fcon/hcpmmp/hcpmmp.1D').pd
path_df = dp(path_df).mutate(fmri_map_path=X.xcp_path + '/roiquant/hcpmmp/' + X.subname +'_hcpmmp_mean.csv').pd
print('finished step 4')

# 5. load individual connectome as ndarray
num_parcels = len(network_assign_csv)
num_sub = len(path_df)
path_df_nc = dp(path_df).filter_by(X.subname.str.contains('nc')).pd
num_nc = len(path_df_nc)
nc_idx = path_df_nc.index
connectome_array = np.zeros(shape=(num_parcels, num_parcels, num_sub))
for sub_idx in range(len(path_df)):
    indiviudal_connectome = np.genfromtxt(path_df.loc[sub_idx, 'connectome_path'], delimiter=',')
    connectome_array[:,:,sub_idx] = indiviudal_connectome

# 6. load individual BOLD series and fill missing part according to /fcon/hcpmmp/missing.txt
BOLD_series_3D = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_sub))
for sub_idx in range(len(path_df)):
    BOLD_series = np.genfromtxt(path_df.loc[sub_idx, 'BOLD_series_path'])
    BOLD_series = BOLD_series.T
    missing_path = os.path.join(path_df.loc[sub_idx, 'xcp_path'], 'fcon', 'hcpmmp', 'hcpmmp_missing.txt')
    if os.path.exists(missing_path):
        missing_parcel_id = np.genfromtxt(missing_path, dtype=int)
        if missing_parcel_id.size == 1: # only one parcel missing
            if BOLD_series[missing_parcel_id-1,:].sum() != 0:
                print("missing parcel not match for subject {}".format(sub_idx))
            network_key = network_assign_csv.loc[missing_parcel_id-1,'NETWORKKEY']
            network_parcel_idx = network_rowindex_ls[network_key-1]
            BOLD_series[missing_parcel_id-1,:] = np.mean(BOLD_series[network_parcel_idx,:])
        else: # multiple parcels missing
            for missing_idx in missing_parcel_id:
                network_key = network_assign_csv.loc[missing_idx-1,'NETWORKKEY']
                network_parcel_idx = network_rowindex_ls[network_key-1]
                BOLD_series[missing_idx-1,:] = np.mean(BOLD_series[network_parcel_idx,:])

    BOLD_series_3D[:,:,sub_idx] = BOLD_series
print('finished loading individual BOLD series and filling missing part')
            
# 7. load fmri parametric map and fill missing part according to /fcon/hcpmmp/missing.txt
fmri_paramap = np.zeros(shape=(num_parcels, num_sub))
paramap_str = 'mean_alffZ'
for sub_idx in range(len(path_df)):
    fmri_map = pd.read_csv(path_df.loc[sub_idx, 'fmri_map_path'],index_col=0)
    fmri_map = fmri_map.loc[:,paramap_str]
    missing_path = os.path.join(path_df.loc[sub_idx, 'xcp_path'], 'fcon', 'hcpmmp', 'hcpmmp_missing.txt')
    if os.path.exists(missing_path):
        missing_parcel_id = np.genfromtxt(missing_path, dtype=int)
        if missing_parcel_id.size == 1: # only one parcel missing
            if not np.isnan(fmri_map[missing_parcel_id]):
                print("missing parcel not match for subject {}".format(sub_idx))
            network_key = network_assign_csv.loc[missing_parcel_id-1,'NETWORKKEY']
            network_parcel_idx = network_rowindex_ls[network_key-1]
            fmri_map[int(missing_parcel_id)] = np.mean(fmri_map[network_parcel_idx])
            fmri_map = fmri_map.to_numpy() 
        else: # multiple parcels missing
            network_key = network_assign_csv.loc[missing_parcel_id-1,'NETWORKKEY']
            network_rowindex_ls = np.array(network_rowindex_ls, dtype=object)
            network_parcel_idx = network_rowindex_ls[network_key-1]
            for parcel_i in range(missing_parcel_id.size):
                fmri_map[int(missing_parcel_id[parcel_i])] = np.mean(fmri_map[network_parcel_idx[parcel_i]])
            fmri_map = fmri_map.to_numpy() 

    fmri_paramap[:,sub_idx] = fmri_map
print('finished loading fmri parametric map and fill missing part')

# 8. load connectome and functional signal and do GSP
if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    func_sig = BOLD_series_3D 
    s_head_cohort = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_sub))
    s_rand_cohort = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_sub, num_rand))
else:
    raise ValueError('undefined functional signal')

G_U_cohort = np.zeros(shape=(num_parcels, num_parcels, num_sub))
for sub_idx in range(len(path_df)):
    W = np.genfromtxt(path_df.loc[sub_idx, 'connectome_path'], delimiter=',')
    # Symmetric Normalization of adjacency matrix
    D = np.diag(np.sum(W,1)) #degree
    D_power = np.power(D, (-1/2))
    D_power[np.isinf(D_power)] = 0
    Wsymm = D_power @ W @ D_power
    
    #The eigenvector matrix G.U is used to define the Graph Fourier Transform of the graph signal S
    if normalize_type == 'W':
        G = graphs.Graph(Wsymm)
        G.compute_fourier_basis() 
        G_U_cohort[:,:,sub_idx] = G.U
        U = G.U
    elif normalize_type == 'L':
        G = graphs.Graph(W, lap_type = 'normalized')
        G.compute_fourier_basis() 
        G_U_cohort[:,:,sub_idx] = G.U
        U = G.U
    elif normalize_type == 'both':
        Wsymm = np.triu(Wsymm) + np.triu(Wsymm).T - np.diag(np.triu(Wsymm).diagonal()) # force symmetric
        G = graphs.Graph(Wsymm, lap_type = 'normalized')
        G.compute_fourier_basis() 
        G_U_cohort[:,:,sub_idx] = G.U
        U = G.U
        # L = np.eye(len(Wsymm)) - Wsymm
        # lamda, U = np.linalg.eig(L)
        # U = U[:, np.argsort(lamda)]

    if functional_type == 'BOLD': # func_sig is BOLD_series_3D
        s_head = U.T @ func_sig[:,:,sub_idx]
        s_head_cohort[:,:,sub_idx] = s_head
        # calcualte surrogate for individual
        s_rand_cohort[:,:,sub_idx,:] = surrogate_BOLD_create(U, func_sig[:,:,sub_idx], num_rand)
print('finished Graph Fourier Transform')    

# save_variable(G_U_cohort, 'G_U_cohort.pkl')
# save_variable(s_head_cohort, 's_head_cohort.pkl')
# save_variable(s_rand_cohort, 's_rand_cohort.pkl')

# G_U_cohort = load_variable('G_U_cohort.pkl')
# s_head_cohort = load_variable('s_head_cohort.pkl')
# s_rand_cohort = load_variable('s_rand_cohort.pkl')


# 8.5(optional). plot Sihag2020 plot 
# take nc001 as example
nc001_idx = path_df.subname[path_df.subname == 'sub-nc001'].index.tolist()[0]
s_low = G_U_cohort[:,0:4, nc001_idx] @ s_head_cohort[0:4,:,nc001_idx]
s_high = G_U_cohort[:,-55:-51, nc001_idx] @ s_head_cohort[-55:-51,:,nc001_idx]
np.savetxt("nc001_s_low_both.csv", s_low, delimiter=",")
np.savetxt("nc001_s_high_both.csv", s_high, delimiter=",")


# 9. calculate the median-split threshold 
NC_index = [cohort_name_ls.index(x) for x in cohort_name_ls if 'nc' in x]

if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    s_head_NC = s_head_cohort[:,:,NC_index]
    s_head_NC_square = np.power(s_head_NC, 2)
    #s_head_NC_square = np.power(s_head_NC_square, 1/2)
    s_head_NC_square_mean = np.mean(s_head_NC_square, (1,2)) # average for each timepoint and each subject


s_head_NC_AUCTOT = np.trapz(s_head_NC_square_mean)

i=0
AUC=0
while AUC < s_head_NC_AUCTOT/2:
    AUC = np.trapz(s_head_NC_square_mean[:i])
    i = i + 1

cutoff = i-1
print('finished calculating the median-split threshold')
print('cutoff = {}'.format(cutoff))

# 10. calculate decoupling index for empirical data
if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    s_aligned_cohort = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_sub))
    s_liberal_cohort = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_sub))
    for sub_idx in range(len(path_df)):
        s_aligned_cohort[:,:,sub_idx] = G_U_cohort[:,0:cutoff, sub_idx] @ s_head_cohort[0:cutoff,:,sub_idx]
        s_liberal_cohort[:,:,sub_idx] = G_U_cohort[:,cutoff-1:-1, sub_idx] @ s_head_cohort[cutoff-1:-1,:,sub_idx]

    s_aligned_individual = np.linalg.norm(s_aligned_cohort, ord=2, axis=1)
    s_liberal_individual = np.linalg.norm(s_liberal_cohort, ord=2, axis=1)
    s_deCoupIdx_individual = s_liberal_individual / s_aligned_individual

    s_aligned = np.mean(s_aligned_individual[:,nc_idx], axis=1)
    s_liberal = np.mean(s_liberal_individual[:,nc_idx], axis=1)
    s_deCoupIdx_node = s_liberal/s_aligned # only for NC
        

print('finished calculating decoupling index for empirical data')


# 11. calculate decoupling index for surrogate data only for NC
if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    s_aligned_cohort_rand = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_nc, num_rand))
    s_liberal_cohort_rand = np.zeros(shape=(num_parcels, num_BOLD_timepoints, num_nc, num_rand))
    for i, sub_idx in enumerate(nc_idx):
        for rand_idx in range(num_rand):
            s_aligned_cohort_rand[:,:,i,rand_idx] = G_U_cohort[:,0:cutoff, sub_idx] @ s_rand_cohort[0:cutoff,:,sub_idx,rand_idx]
            s_liberal_cohort_rand[:,:,i,rand_idx] = G_U_cohort[:,cutoff-1:-1, sub_idx] @ s_rand_cohort[cutoff-1:-1,:,sub_idx,rand_idx]
    # norm for BOLD timepoints
    s_aligned_norm_rand = np.linalg.norm(s_aligned_cohort_rand, ord=2, axis=1) 
    s_liberal_norm_rand = np.linalg.norm(s_liberal_cohort_rand, ord=2, axis=1)
    # average for cohorts
    s_aligned_rand = np.mean(s_aligned_norm_rand, axis=1)
    s_liberal_rand = np.mean(s_liberal_norm_rand, axis=1)
    # decoupling index
    s_deCoupIdx_node_rand = s_liberal_rand/s_aligned_rand

print('finished calculating decoupling index for surrogate data')


# 12. network-level harmonics for emperical and surrogate data
s_aligned_network = np.zeros(shape=(num_network))
s_liberal_network = np.zeros(shape=(num_network))
s_aligned_network_individual = np.zeros(shape=(num_network, num_sub))
s_liberal_network_individual = np.zeros(shape=(num_network, num_sub))
s_aligned_network_rand = np.zeros(shape=(num_network, num_rand))
s_liberal_network_rand = np.zeros(shape=(num_network, num_rand))
for i in range(num_network):
    s_aligned_network[i] = np.mean(s_aligned[network_rowindex_ls[i]])
    s_liberal_network[i] = np.mean(s_liberal[network_rowindex_ls[i]])

    s_aligned_network_individual[i,:] = np.mean(s_aligned_individual[network_rowindex_ls[i],:], axis=0)
    s_liberal_network_individual[i,:] = np.mean(s_liberal_individual[network_rowindex_ls[i],:], axis=0)

    s_aligned_network_rand[i,:] = np.mean(s_aligned_rand[network_rowindex_ls[i],:], axis=0)
    s_liberal_network_rand[i,:] = np.mean(s_liberal_rand[network_rowindex_ls[i],:], axis=0)

s_deCoupIdx_network = s_liberal_network/s_aligned_network
s_deCoupIdx_network_individual = s_liberal_network_individual/s_aligned_network_individual
s_deCoupIdx_network_rand = s_liberal_network_rand/s_aligned_network_rand


# 13. brain-level harmonics for emperical and surrogate data
s_aligned_brain = np.mean(s_aligned)
s_liberal_brain = np.mean(s_liberal)
s_deCoupIdx_brain = s_liberal_brain/s_aligned_brain

s_aligned_brain_individual = np.mean(s_aligned_individual, axis=0)
s_liberal_brain_individual = np.mean(s_liberal_individual, axis=0)
s_deCoupIdx_brain_individual = s_liberal_brain_individual/s_aligned_brain_individual

s_aligned_brain_rand = np.mean(s_aligned_rand, axis=0)
s_liberal_brain_rand = np.mean(s_liberal_rand, axis=0)
s_deCoupIdx_brain_rand = s_liberal_brain_rand/s_aligned_brain_rand

print('s_deCoupIdx_brain = {}'.format(s_deCoupIdx_brain))

# 14. significance of surrogate for plot
# node-level
s_deCoupIdx_node_significance = np.logical_or((np.percentile(s_deCoupIdx_node_rand, 5, axis=1) >= s_deCoupIdx_node), (np.percentile(s_deCoupIdx_node_rand, 95, axis=1) <= s_deCoupIdx_node))
s_deCoupIdx_node_significance = s_deCoupIdx_node_significance.astype(np.int)

# network-level
s_deCoupIdx_network_significance = np.logical_or((np.percentile(s_deCoupIdx_network_rand, 5, axis=1) >= s_deCoupIdx_network), (np.percentile(s_deCoupIdx_network_rand, 95, axis=1) <= s_deCoupIdx_network))
s_deCoupIdx_network_significance = s_deCoupIdx_network_significance.astype(np.int)

# brain-level
s_deCoupIdx_brain_significance = np.logical_or((np.percentile(s_deCoupIdx_brain_rand, 5, axis=0) >= s_deCoupIdx_brain), (np.percentile(s_deCoupIdx_brain_rand, 95, axis=0) <= s_deCoupIdx_brain))

# 15. save results to csv
if normalize_type == 'W':
    normalize_str = '_W'
elif normalize_type == 'L':
    normalize_str = '_L'
elif normalize_type == 'both':
    normalize_str = '_both'
if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    csv_folder = 'BOLD_4D_'  + tract_type + '_'  + normalize_str

if not os.path.exists(os.path.abspath(csv_folder)):
    os.mkdir(os.path.abspath(csv_folder))


# save surrogate (ndarray with num_rand × num_region)
s_deCoupIdx_node_rand_df = pd.DataFrame(data = s_deCoupIdx_node_rand.T, columns = network_assign_csv.loc[:,'LABEL'])
s_deCoupIdx_network_rand_df = pd.DataFrame(data = s_deCoupIdx_network_rand.T, columns = network_unique_df.loc[:,'NETWORK'])
s_deCoupIdx_brain_rand_df = pd.DataFrame(data = s_deCoupIdx_brain_rand)

s_deCoupIdx_node_rand_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_node_rand_df.csv'))
s_deCoupIdx_network_rand_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_' + '-network_rand_df.csv'))
s_deCoupIdx_brain_rand_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_brain_rand_df.csv'))

# save surrogate significance (ndarray with 1 × num_region)
s_deCoupIdx_node_significance_df = pd.DataFrame(data = np.expand_dims(s_deCoupIdx_node_significance, axis=0), columns = network_assign_csv.loc[:,'LABEL'])
s_deCoupIdx_network_significance_df = pd.DataFrame(data = np.expand_dims(s_deCoupIdx_network_significance, axis=0), columns = network_unique_df.loc[:,'NETWORK'])

s_deCoupIdx_node_significance_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_node_significance_df.csv'))
s_deCoupIdx_network_significance_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_' + '-network_significance_df.csv'))
with open(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_brain_significance.txt'), 'w') as output_file:
    output_file.write(str(s_deCoupIdx_brain_significance))

# save empirical harmonics for NC cohort (for plot usage, ndarray with 1 × num_region)
s_deCoupIdx_node_empirical_df = pd.DataFrame(data = np.expand_dims(s_deCoupIdx_node, axis=0), columns = network_assign_csv.loc[:,'LABEL'])
s_deCoupIdx_network_empirical_df = pd.DataFrame(data = np.expand_dims(s_deCoupIdx_network, axis=0), columns = network_unique_df.loc[:,'NETWORK'])

s_deCoupIdx_node_empirical_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_node_empirical_df.csv'))
s_deCoupIdx_network_empirical_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_' +'-network_empirical_df.csv'))
with open(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_brain_empirical.txt'), 'w') as output_file:
    output_file.write(str(s_deCoupIdx_brain))

# save subject-level harmonics (ndarray with num_sub × num_region)
s_deCoupIdx_node_individual_df = pd.DataFrame(data = s_deCoupIdx_individual.T, columns = network_assign_csv.loc[:,'LABEL'])
s_deCoupIdx_network_individual_df = pd.DataFrame(data = s_deCoupIdx_network_individual.T, columns = network_unique_df.loc[:,'NETWORK'])
s_deCoupIdx_brain_individual_df = pd.DataFrame(data = s_deCoupIdx_brain_individual)

s_deCoupIdx_node_individual_df = pd.concat([path_df.loc[:,'subname'],s_deCoupIdx_node_individual_df],axis=1)
s_deCoupIdx_network_individual_df = pd.concat([path_df.loc[:,'subname'],s_deCoupIdx_network_individual_df],axis=1)
s_deCoupIdx_brain_individual_df = pd.concat([path_df.loc[:,'subname'],s_deCoupIdx_brain_individual_df],axis=1)

s_deCoupIdx_node_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_node_individual_df.csv'))
s_deCoupIdx_network_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_' + '-network_individual_df.csv'))
s_deCoupIdx_brain_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 's_deCoupIdx_brain_individual_df.csv'))


# 16.(optional) save connectivity strength
# parcel-level
connectome_parcel_individual = np.zeros(shape=(num_sub, num_parcels))
# mean of nonzero
def non_zero_mean(np_arr):
    exist = (np_arr != 0)
    num = np_arr.sum(axis=1)
    den = exist.sum(axis=1)
    return num/den

for sub_idx in range(num_sub):
    connectome_parcel_individual[sub_idx,:] = non_zero_mean(connectome_array[:,:,sub_idx])

connectome_parcel_individual_df = pd.DataFrame(data = connectome_parcel_individual, columns = network_assign_csv.loc[:,'LABEL'])
connectome_parcel_individual_df = pd.concat([path_df.loc[:,'subname'], connectome_parcel_individual_df],axis=1)
connectome_parcel_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 'connectome_' +  '-parcel_individual_df.csv'))

# ICN-level
connectome_network_individual = np.zeros(shape=(num_network, num_sub))
for i in range(num_network):
    network_i = network_unique_df.loc[i,'NETWORK']
    parcel_network_df = dp(network_assign_csv).filter_by(X.NETWORK.isin([network_i])).pd
    parcel_network_id = parcel_network_df.loc[:,'INDEX'].to_numpy()
    connectome_network_individual[i,:] = np.mean(connectome_array[np.ix_(parcel_network_id-1,parcel_network_id-1)], axis=(0,1))


connectome_network_individual_df = pd.DataFrame(data = connectome_network_individual.T, columns = network_unique_df.loc[:,'NETWORK'])
connectome_network_individual_df = pd.concat([path_df.loc[:,'subname'], connectome_network_individual_df],axis=1)
connectome_network_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 'connectome_' +'-network_individual_df.csv'))

# 17.(optional) save ICN-level brain activity
if functional_type == 'BOLD': # func_sig is BOLD_series_3D
    BOLD_network_individual = np.zeros(shape=(num_network, num_sub))
    for i in range(num_network):
        network_i = network_unique_df.loc[i,'NETWORK']
        parcel_network_df = dp(network_assign_csv).filter_by(X.NETWORK.isin([network_i])).pd
        parcel_network_id = parcel_network_df.loc[:,'INDEX'].to_numpy()
        BOLD_series_norm = np.linalg.norm(BOLD_series_3D, ord=2, axis=1)
        BOLD_network_individual[i,:] =  np.mean(BOLD_series_norm[parcel_network_id-1,:], axis=0)
    BOLD_network_individual_df = pd.DataFrame(data = BOLD_network_individual.T, columns = network_unique_df.loc[:,'NETWORK'])
    BOLD_network_individual_df = pd.concat([path_df.loc[:,'subname'], BOLD_network_individual_df],axis=1)
    BOLD_network_individual_df.to_csv(os.path.join(os.path.abspath(csv_folder), 'BOLD_' +  '-network_individual_df.csv'))


end = time.time()
running_time = end-start
print('time cost : %.2f sec' %running_time)
print('done')
