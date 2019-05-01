"""

    Generate graphs with results for Arp2/3 project

    Input:  - 2nd order results for W and IT
            - 1st order results for W and IT

    Output: - Plots with the results

"""

################# Package import

import os
import pickle
import numpy as np
import scipy as sp
import sys
import time
from surf_dst import pexceptions, sub, disperse_io, surf
from surf_dst.globals import unpickle_obj, sort_dict
from surf_dst.surf.model import ModelCSRV
from surf_dst.surf.utils import list_tomoparticles_pvalues
from surf_dst.spatial.sparse import compute_hist
from matplotlib import pyplot as plt, rcParams

###### Global variables

__author__ = 'Antonio Martinez-Sanchez'

BAR_WIDTH = .35

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20

########################################################################################
# PARAMETERS
########################################################################################

ROOT_PATH = '/fs/pool/pool-lucic2/antonio/tomograms/marion/Arp23complex'

# Input STAR files
in_star_w = ROOT_PATH + '/ltomos_all/W_all_mask/W_all_mask_ltomos.star' 
in_star_it = ROOT_PATH + '/ltomos_all/IT_all_mask/IT_all_mask_ltomos.star' 

# Input matrices (optional - organization analysis is skipped)
in_mats_w2nd = ROOT_PATH + '/uni_sph_all_W/W_all_sim200_org_lists.pkl' 
in_mats_it2nd = ROOT_PATH + '/uni_sph_all_IT/IT_all_sim200_org_lists.pkl' 
in_sims_w2nd = ROOT_PATH + '/uni_sph_all_W/W_all_sim200_org_sims.pkl' 
in_sims_it2nd = ROOT_PATH + '/uni_sph_all_IT/IT_all_sim200_org_sims.pkl'
in_mats_w1st = ROOT_PATH + '/uni_sph_1st_all_W/W_all_sim20_wspace.pkl'
in_mats_it1st = ROOT_PATH + '/uni_sph_1st_all_IT/IT_all_sim20_wspace.pkl'  

# Output directory
out_dir = ROOT_PATH + '/plots'

# Analysis variables
ana_res = 1.684 # nm/voxel
ana_rg = np.arange(5, 1100, 10) # np.arange(5, 800, 10) # np.arange(4, 100, 2) # in nm
# P-value computation settings
# Simulation model (currently only CSRV)
p_per = 1 # 5 # %

# Firts order analysis
ana_nbins = 20
ana_rmax = 200 

# Figure saving options
fig_fmt = '.png' # if None they showed instead

# Plotting options
pt_xrange = [0, 1100] # [0, 800]
pt_dxrange = [0, 380]
pt_yrange = [-25, 25]
pt_cmap = plt.get_cmap('gist_rainbow')

###### Additional functionality

# Computes IC from a matrix of measurements (n_arrays, array_samples)
def compute_ic(per, sims):
    if len(sims.shape) == 1:
        return sims, sims, sims
    ic_low = np.percentile(sims, per, axis=0, interpolation='linear')
    ic_med = np.percentile(sims, 50, axis=0, interpolation='linear')
    ic_high = np.percentile(sims, 100 - per, axis=0, interpolation='linear')
    return ic_low, ic_med, ic_high

# Computes pvalue from a matrix of simulations (n_arrays, array_samples)
def compute_pvals(exp_med, sims):
    n_sims = float(sims.shape[0])
    p_vals_low, p_vals_high = np.zeros(shape=exp_med.shape, dtype=np.float32), \
                              np.zeros(shape=exp_med.shape, dtype=np.float32)
    for i, exp in enumerate(exp_med):
        sim_slice = sims[:, i]
        p_vals_high[i] = float((exp > sim_slice).sum()) / n_sims
        p_vals_low[i] = float((exp < sim_slice).sum()) / n_sims
    return p_vals_low, p_vals_high

########################################################################################
# MAIN ROUTINE
########################################################################################

print '\tLoading input STAR file...'
star_w, star_it = sub.Star(), sub.Star()
try:
    star_w.load(in_star_w)
    star_it.load(in_star_it)
except pexceptions.PySegInputError as e:
    print 'ERROR: input STAR file could not be loaded because of "' + e.get_message() + '"'
    print 'Terminated. (' + time.strftime("%c") + ')'
    sys.exit(-1)
set_lists_w, set_lists_it = surf.SetListTomoParticles(), surf.SetListTomoParticles()
for row in range(star_w.get_nrows()):
    ltomos_pkl = star_w.get_element('_psPickleFile', row)
    ltomos = unpickle_obj(ltomos_pkl)
    set_lists_w.add_list_tomos(ltomos, ltomos_pkl)
for row in range(star_it.get_nrows()):
    ltomos_pkl = star_it.get_element('_psPickleFile', row)
    ltomos = unpickle_obj(ltomos_pkl)
    set_lists_it.add_list_tomos(ltomos, ltomos_pkl)

print '\tComputing densities...'
dens_w = np.asarray(set_lists_w.density_by_tomos(surface=False).values(), dtype=float)
dens_it = np.asarray(set_lists_it.density_by_tomos(surface=False).values(), dtype=float)
means = np.asarray((dens_w.mean(), dens_it.mean()))
stds = np.asarray((dens_w.std()/np.sqrt(float(dens_w.shape[0])), dens_it.std()/np.sqrt(float(dens_it.shape[0]))))
plt.figure()
plt.ylabel('Density [particles/nm$^3$]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.bar(0.25, means[0], 0.5, color='orange', yerr=stds[0], ecolor='k', linewidth=4, error_kw={'elinewidth':4, 'capthick':4})
plt.bar(1.25, means[1], 0.5, color='springgreen', yerr=stds[1], ecolor='k', linewidth=4, error_kw={'elinewidth':4, 'capthick':4})
plt.xticks((0.5, 1.5), ('W', 'IT'))
plt.xlim(0, 2)
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/W_IT_density' + fig_fmt, dpi=600)
plt.close()

print '\tPickling input matrices...'
if in_mats_w2nd is not None:
    with open(in_mats_w2nd, 'r') as pkl:
        mats_w2nd = pickle.load(pkl)
if in_mats_it2nd is not None:
    with open(in_mats_it2nd, 'r') as pkl:
        mats_it2nd = pickle.load(pkl)
if in_sims_w2nd is not None:
    with open(in_sims_w2nd, 'r') as pkl:
        sims_w2nd = pickle.load(pkl)
if in_sims_it2nd is not None:
    with open(in_sims_it2nd, 'r') as pkl:
        sims_it2nd = pickle.load(pkl)

if mats_w2nd is not None:
    plt.figure()
    plt.ylabel('Ripley\'s L')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_w2nd.values()[0], sims_w2nd.values()[0]
    plt.plot(ana_rg, mat, linewidth=3, color='blue', label='W-experimental') 
    ic_low, ic_med, ic_high = np.percentile(sims, p_per, axis=0), np.median(sims, axis=0), np.percentile(sims, 100-p_per, axis=0)
    ic_low_5, ic_high_95 = np.percentile(sims, 5, axis=0), np.percentile(sims, 95, axis=0)
    plt.plot(ana_rg, ic_med, linewidth=3, color='gray', label='IC-simulated') 
    plt.plot(ana_rg, ic_low_5, linewidth=2, color='k', linestyle='--')
    plt.plot(ana_rg, ic_high_95, linewidth=2, color='k', linestyle='--') 
    plt.fill_between(ana_rg, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
    x_max, y_max = ana_rg[mat.argmax()], mat.max()
    plt.plot((x_max, x_max), (pt_yrange[0], y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_max), (y_max, y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, x_max, 400, 600, 800, 1000))
    plt.yticks((-20, -10, 0, 10, y_max))
    if pt_xrange is not None:
        plt.xlim(pt_xrange)
    if pt_yrange is not None:
        plt.ylim(pt_yrange)
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/W_IC_RipleyL' + fig_fmt, dpi=600)
    plt.close()

    flt_ncoefs = 9
    plt.figure()
    plt.ylabel('Ripley\'s L')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_w2nd.values()[0], sims_w2nd.values()[0]
    mat = sp.signal.savgol_filter(mat, flt_ncoefs, 2, mode='interp')
    plt.plot(ana_rg, mat, linewidth=3, color='blue', label='W-experimental') 
    ic_low, ic_med, ic_high = np.percentile(sims, p_per, axis=0), np.median(sims, axis=0), np.percentile(sims, 100-p_per, axis=0)
    ic_low_5, ic_high_95 = np.percentile(sims, 5, axis=0), np.percentile(sims, 95, axis=0)
    ic_low = sp.signal.savgol_filter(ic_low, flt_ncoefs, 2, mode='interp')
    ic_med = sp.signal.savgol_filter(ic_med, flt_ncoefs, 2, mode='interp')
    ic_high = sp.signal.savgol_filter(ic_high, flt_ncoefs, 2, mode='interp')
    ic_low_5 = sp.signal.savgol_filter(ic_low_5, flt_ncoefs, 2, mode='interp')
    ic_high_95 = sp.signal.savgol_filter(ic_high_95, flt_ncoefs, 2, mode='interp')
    plt.plot(ana_rg, ic_med, linewidth=3, color='gray', label='IC-simulated') 
    plt.plot(ana_rg, ic_low_5, linewidth=2, color='k', linestyle='--')
    plt.plot(ana_rg, ic_high_95, linewidth=2, color='k', linestyle='--') 
    plt.fill_between(ana_rg, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
    x_max, y_max = ana_rg[mat.argmax()], mat.max()
    plt.plot((x_max, x_max), (pt_yrange[0], y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_max), (y_max, y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, x_max, 400, 600, 800, 1000))
    plt.yticks((-20, -10, 0, 10, y_max))
    if pt_xrange is not None:
        plt.xlim(pt_xrange)
    if pt_yrange is not None:
        plt.ylim(pt_yrange)
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/W_IC_RipleyL_sgflt' + fig_fmt, dpi=600)
    plt.close()

    plt.figure()
    plt.ylabel('Ripley\'s L\'')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_w2nd.values()[0], sims_w2nd.values()[0]
    mat = np.gradient(sp.signal.savgol_filter(mat, flt_ncoefs, 2, mode='interp'), ana_rg[1]-ana_rg[0])
    ana_rg_c = ana_rg[(ana_rg >= pt_dxrange[0]) & (ana_rg <= pt_dxrange[1])]
    mat_c = mat[(ana_rg >= pt_dxrange[0]) & (ana_rg <= pt_dxrange[1])]
    plt.plot(ana_rg_c, mat_c, linewidth=3, color='blue', label='W-experimental') 
    x_min, y_min = ana_rg_c[mat_c.argmin()], mat_c.min()
    plt.plot((x_min, x_min), (-.1, y_min), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_min), (y_min, y_min), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, 100, x_min, 250))
    plt.yticks((-.1, y_min, 0, .2, .4))
    if pt_xrange is not None:
        plt.xlim(pt_dxrange)
    if pt_yrange is not None:
        plt.ylim((-.1, .5))
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/W_IC_RipleyL_sgdif' + fig_fmt, dpi=600)
    plt.close()
else:
    print 'ERROR: organization could not be computed'
    print 'Unsuccessfully terminated. (' + time.strftime("%c") + ')'
    sys.exit(-1)

if mats_it2nd is not None:
    plt.figure()
    plt.ylabel('Ripley\'s L')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_it2nd.values()[0], sims_it2nd.values()[0]
    plt.plot(ana_rg, mat, linewidth=3, color='blue', label='IT-experimental') 
    ic_low, ic_med, ic_high = np.percentile(sims, p_per, axis=0), np.median(sims, axis=0), np.percentile(sims, 100-p_per, axis=0)
    ic_low_5, ic_high_95 = np.percentile(sims, 5, axis=0), np.percentile(sims, 95, axis=0)    
    plt.plot(ana_rg, ic_med, linewidth=3, color='gray', label='IT-simulated') 
    plt.plot(ana_rg, ic_low_5, linewidth=2, color='k', linestyle='--')
    plt.plot(ana_rg, ic_high_95, linewidth=2, color='k', linestyle='--') 
    plt.fill_between(ana_rg, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
    x_max, y_max = ana_rg[mat.argmax()], mat.max()
    plt.plot((x_max, x_max), (pt_yrange[0], y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_max), (y_max, y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, x_max, 400, 600, 800, 1000))
    plt.yticks((-20, -10, 0, y_max, 20))
    if pt_xrange is not None:
        plt.xlim(pt_xrange)
    if pt_yrange is not None:
        plt.ylim(pt_yrange)
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/IT_IC_RipleyL' + fig_fmt, dpi=600)
    plt.close()

    flt_ncoefs = 9
    plt.figure()
    plt.ylabel('Ripley\'s L')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_it2nd.values()[0], sims_it2nd.values()[0]
    mat = sp.signal.savgol_filter(mat, flt_ncoefs, 2, mode='interp')
    plt.plot(ana_rg, mat, linewidth=3, color='blue', label='W-experimental') 
    ic_low, ic_med, ic_high = np.percentile(sims, p_per, axis=0), np.median(sims, axis=0), np.percentile(sims, 100-p_per, axis=0)
    ic_low_5, ic_high_95 = np.percentile(sims, 5, axis=0), np.percentile(sims, 95, axis=0)
    ic_low = sp.signal.savgol_filter(ic_low, flt_ncoefs, 2, mode='interp')
    ic_med = sp.signal.savgol_filter(ic_med, flt_ncoefs, 2, mode='interp')
    ic_high = sp.signal.savgol_filter(ic_high, flt_ncoefs, 2, mode='interp')
    ic_low_5 = sp.signal.savgol_filter(ic_low_5, flt_ncoefs, 2, mode='interp')
    ic_high_95 = sp.signal.savgol_filter(ic_high_95, flt_ncoefs, 2, mode='interp')
    plt.plot(ana_rg, ic_med, linewidth=3, color='gray', label='IC-simulated') 
    plt.plot(ana_rg, ic_low_5, linewidth=2, color='k', linestyle='--')
    plt.plot(ana_rg, ic_high_95, linewidth=2, color='k', linestyle='--') 
    plt.fill_between(ana_rg, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
    x_max, y_max = ana_rg[mat.argmax()], mat.max()
    plt.plot((x_max, x_max), (pt_yrange[0], y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_max), (y_max, y_max), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, x_max, 400, 600, 800, 1000))
    plt.yticks((-20, -10, 0, 10, y_max))
    if pt_xrange is not None:
        plt.xlim(pt_xrange)
    if pt_yrange is not None:
        plt.ylim(pt_yrange)
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/IT_IC_RipleyL_sgflt' + fig_fmt, dpi=600)
    plt.close()

    plt.figure()
    plt.ylabel('Ripley\'s L\'')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_it2nd.values()[0], sims_it2nd.values()[0]
    mat = np.gradient(sp.signal.savgol_filter(mat, flt_ncoefs, 2, mode='interp'), ana_rg[1]-ana_rg[0])
    ana_rg_c = ana_rg[(ana_rg >= pt_dxrange[0]) & (ana_rg <= pt_dxrange[1])]
    mat_c = mat[(ana_rg >= pt_dxrange[0]) & (ana_rg <= pt_dxrange[1])]
    plt.plot(ana_rg_c, mat_c, linewidth=3, color='blue', label='W-experimental') 
    x_min, y_min = ana_rg_c[mat_c.argmin()], mat_c.min()
    plt.plot((x_min, x_min), (-.1, y_min), linewidth=2, marker='o', linestyle='--', color='k')
    plt.plot((pt_xrange[0], x_min), (y_min, y_min), linewidth=2, marker='o', linestyle='--', color='k')
    plt.xticks((0, 100, x_min, 250))
    plt.yticks((-.1, y_min, 0, .2, .4))
    if pt_xrange is not None:
        plt.xlim(pt_dxrange)
    if pt_yrange is not None:
        plt.ylim((-.1, .5))
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/IT_IC_RipleyL_sgdif' + fig_fmt, dpi=600)
    plt.close()
else:
    print 'ERROR: organization could not be computed'
    print 'Unsuccessfully terminated. (' + time.strftime("%c") + ')'
    sys.exit(-1)

if (mats_w2nd is not None) and (mats_it2nd is not None):
    plt.figure()
    plt.ylabel('Ripley\'s L')
    plt.xlabel('Distance [nm]')
    mat, sims = mats_w2nd.values()[0], sims_w2nd.values()[0]
    plt.plot(ana_rg, mat, linewidth=4, color='orange', label='W-experimental') 
    mat, sims = mats_it2nd.values()[0], sims_it2nd.values()[0]
    plt.plot(ana_rg, mat, linewidth=4, color='springgreen', label='IT-experimental') 
    if pt_xrange is not None:
        plt.xlim(pt_xrange)
    if pt_yrange is not None:
        plt.ylim(pt_yrange)
    # plt.legend()
    plt.tight_layout()
    if fig_fmt is None:
        plt.show(block=True)
    else:
        plt.savefig(out_dir + '/W_IT_RipleyL' + fig_fmt, dpi=600)
    plt.close()
else:
    print 'ERROR: organization could not be computed'
    print 'Unsuccessfully terminated. (' + time.strftime("%c") + ')'
    sys.exit(-1)

print '\nLoading 1st order analysis...'
with open(in_mats_w1st, 'r') as pkl:
    wspace = pickle.load(pkl)
lists_count_w, tomos_count_w = wspace[0], wspace[1]
lists_hash_w, tomos_hash_w = wspace[2], wspace[3]
tomos_exp_dsts_w, tomos_sim_dsts_w, tomos_exp_fdsts_w, tomos_sim_fdsts_w = wspace[4], wspace[5], wspace[6], wspace[7]
lists_exp_dsts_w, lists_sim_dsts_w, lists_exp_fdsts_w, lists_sim_fdsts_w, lists_color_w = wspace[8], wspace[9], wspace[10], wspace[11], wspace[12]
with open(in_mats_it1st, 'r') as pkl:
    wspace = pickle.load(pkl)
lists_count_it, tomos_count_it = wspace[0], wspace[1]
lists_hash_it, tomos_hash_it = wspace[2], wspace[3]
tomos_exp_dsts_it, tomos_sim_dsts_it, tomos_exp_fdsts_it, tomos_sim_fdsts_it = wspace[4], wspace[5], wspace[6], wspace[7]
lists_exp_dsts_it, lists_sim_dsts_it, lists_exp_fdsts_it, lists_sim_fdsts_it, lists_color_it = wspace[8], wspace[9], wspace[10], wspace[11], wspace[12]

print '\t\t-Plotting Histogram...'
ltomo = lists_exp_dsts_w.values()[0]
hist_bins, hist_vals_w = compute_hist(np.concatenate(np.asarray(ltomo)), ana_nbins, ana_rmax)
list_sim_dsts = lists_sim_dsts_w.values()[0]
sims_hist_vals = list()
for sim_dsts in list_sim_dsts:
    sims_hist_vals.append(compute_hist(sim_dsts, ana_nbins, ana_rmax)[1])
if len(sims_hist_vals) > 0:
    ic_low, ic_med, ic_high = compute_ic(p_per, np.asarray(sims_hist_vals))
plt.figure()
# plt.ylabel('Nearest neighbor distribution')
plt.ylabel('Nearest Neighbor Probability')
plt.xlabel('Distance [nm]')
plt.plot(hist_bins, hist_vals_w, color='blue', linewidth=2, label='W-experimental', marker='o')
plt.plot(hist_bins, ic_med, 'gray', linewidth=2, label='W-simulated')
plt.fill_between(hist_bins, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
plt.xlim(0, ana_rmax)
plt.ylim(0, 0.035)
# plt.legend()
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/W_IC_histogram' + fig_fmt)
plt.close()

print '\t\t-Plotting Histogram...'
ltomo = lists_exp_dsts_it.values()[0]
hist_bins, hist_vals_it = compute_hist(np.concatenate(np.asarray(ltomo)), ana_nbins, ana_rmax)
list_sim_dsts = lists_sim_dsts_it.values()[0]
sims_hist_vals = list()
for sim_dsts in list_sim_dsts:
    sims_hist_vals.append(compute_hist(sim_dsts, ana_nbins, ana_rmax)[1])
if len(sims_hist_vals) > 0:
    ic_low, ic_med, ic_high = compute_ic(p_per, np.asarray(sims_hist_vals))
plt.figure()
plt.ylabel('Nearest Neighbor Probability')
plt.xlabel('Distance [nm]')
plt.plot(hist_bins, hist_vals_it, color='blue', linewidth=2, label='IT-experimental', marker='o')
plt.plot(hist_bins, ic_med, 'gray', linewidth=2, label='IT-simulated')
plt.fill_between(hist_bins, ic_low, ic_high, alpha=0.5, color='gray', edgecolor='w')
plt.xlim(0, ana_rmax)
plt.ylim(0, 0.035)
# plt.legend()
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/IT_IC_histogram' + fig_fmt)
plt.close()

plt.figure()
plt.ylabel('Nearest Neighbor Probability')
plt.xlabel('Distance [nm]')
plt.plot(hist_bins, hist_vals_w, color='blue', linestyle='-', marker='o', linewidth=2, label='W-experimental')
plt.plot(hist_bins, hist_vals_it, color='blue', linestyle='--', marker='s', linewidth=2, label='IT-experimental')
plt.xlim(0, ana_rmax)
plt.xticks(hist_bins)
# plt.legend()
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/W_IT_histogram' + fig_fmt)
plt.close()

plt.figure()
plt.ylabel('Probability')
plt.xlabel('Nearest Neighbor Distance [nm]')
H, B = np.histogram(np.concatenate(np.asarray(lists_exp_dsts_w.values()[0])), bins=np.arange(0,ana_rmax,10), range=None, normed=True)
plt.bar(B[:-1], 10.*H, width=10, color='b', linewidth=2)
# plt.hist(np.concatenate(np.asarray(lists_exp_dsts_w.values()[0])), bins=np.arange(0,ana_rmax,10), range=None, normed=10)
# plt.plot(hist_bins, hist_vals_w, color='blue', linestyle='-', marker='o', linewidth=2, label='W-experimental')
# plt.xlim(0, ana_rmax)
# plt.ylim(0, 0.02)
# plt.xticks(hist_bins)
# plt.legend()
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/W_histogram' + fig_fmt)
plt.close()

plt.figure()
plt.ylabel('Probability')
plt.xlabel('Nearest Neighbor Distance [nm]')
H, B = np.histogram(np.concatenate(np.asarray(lists_exp_dsts_it.values()[0])), bins=np.arange(0,ana_rmax,10), range=None, normed=True)
plt.bar(B[:-1], 10.*H, width=10, color='b', linewidth=2)
# plt.hist(np.concatenate(np.asarray(lists_exp_dsts_it.values()[0])), bins=np.arange(0,ana_rmax,10), range=None, normed=10)
# plt.plot(hist_bins, hist_vals_it, color='blue', linestyle='-', marker='s', linewidth=2, label='IT-experimental')
# plt.xlim(0, ana_rmax)
# plt.ylim(0, 0.02)
# plt.xticks(hist_bins)
# plt.legend()
plt.tight_layout()
if fig_fmt is None:
    plt.show(block=True)
else:
    plt.savefig(out_dir + '/IT_histogram' + fig_fmt)
plt.close()

print 'Successfully terminated. (' + time.strftime("%c") + ')'
