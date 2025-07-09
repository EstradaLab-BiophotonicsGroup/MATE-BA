# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:07:44 2024

@author: Ignacio Sallaberry
"""

import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)

from Functions import *

#%%
intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 5\Simulation Fig 5 - 1 mers.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
kimogram_1_mers = Kimogram(intensity_values,64)

intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 5\Simulation Fig 5 - 2 mers.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
kimogram_2_mers = Kimogram(intensity_values,64)

intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 5\Simulation Fig 5 - 4 mers.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
kimogram_4_mers = Kimogram(intensity_values,64)

intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 5\Simulation Fig 5 - 8 mers.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
kimogram_8_mers = Kimogram(intensity_values,64)

#%%
t_p=1e-6

#%%
## Molecular brightness calculated from N = 100 consecutive correlated lines                                                                                            
# ---------------------------------------------------------------------------------

window_points = 100 
window_shift = window_points//4 
delta_line = 1

Epsilon_kimograms_100_corr_lines = []
for i in (kimogram_1_mers, kimogram_2_mers, kimogram_4_mers, kimogram_8_mers):
    Epsilon_kimograms_100_corr_lines.append(B_apparent_brightness_consecutives_lines(i, t_p,window_points=window_points, window_shift=window_shift)[0].ravel())
Epsilon_kimograms_100_corr_lines = (np.array(Epsilon_kimograms_100_corr_lines)-1)/0.3536/t_p

#%%
## Molecular brightness calculated from N = 50,000 consecutive correlated lines                                                                                            
# ---------------------------------------------------------------------------------

window_points = 50000
window_shift = window_points//400 
delta_line = 1

Epsilon_kimograms_50k_corr_lines = []
for i in (kimogram_1_mers, kimogram_2_mers, kimogram_4_mers, kimogram_8_mers):
    Epsilon_kimograms_50k_corr_lines.append(B_apparent_brightness_consecutives_lines(i, t_p, window_points=window_points, window_shift=window_shift)[0].ravel())
Epsilon_kimograms_50k_corr_lines = (np.array(Epsilon_kimograms_50k_corr_lines)-1)/0.3536/t_p

#%%
## Molecular brightness calculated from N = 100 uncorrelated (interspersed) lines                                                                                            
# ---------------------------------------------------------------------------------

window_points = 100 
window_shift = 1 
delta_line = int(len(kimogram_1_mers)/window_points)    

Epsilon_kimograms_100_lines_MATE_BA = []
for i in (kimogram_1_mers, kimogram_2_mers, kimogram_4_mers, kimogram_8_mers):
    Epsilon_kimograms_100_lines_MATE_BA.append(MATE_BA(i, t_p, window_points=window_points, delta_line=delta_line, window_shift=window_shift)[0].ravel())
Epsilon_kimograms_100_lines_MATE_BA = (np.array(Epsilon_kimograms_100_lines_MATE_BA)-1)/0.3536/t_p

#%%
## Load save data

# Epsilon_kimograms_100_corr_lines = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 5 - Epsilon_kimograms_100_corr_lines.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Epsilon_kimograms_50k_corr_lines = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 5 - Epsilon_kimograms_50k_corr_lines.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Epsilon_kimograms_100_lines_MATE_BA = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 5 - Epsilon_kimograms_100_lines_MATE_BA.csv.gz', compression='gzip', dtype=np.float64).to_numpy()

# %%
plt.close('all')

## General data parameters for plotting
colors = ['#FA003F', '#03CEA4', '#8700F5', '#F9AF10']
edgecolor__ = [colors, colors, [None]*4]
letters = ['A', 'B', 'C']
simulated_epsilon = [0.5/t_p, 1/t_p, 2/t_p, 4/t_p] 
hatch_fill_ = [['//////']*4, ['//////']*4, [None]*4]
histtype_ =[['step']*4, ['step']*4, ['bar']*4]

bins=np.arange(-1.25/t_p,5.75/t_p,0.025/t_p)
density = False
log = False

## Fig 5
fig, axs = plt.subplots(3, 1, figsize=(3,3), sharex=True, sharey=True, dpi=300)

for i, ax, letter, hatch_fill__, histtype__, edgecolor_, epsilons in zip([0,1,2], axs.flatten(), letters, hatch_fill_, histtype_, edgecolor__,
                                                             [Epsilon_kimograms_100_corr_lines,
                                                              Epsilon_kimograms_50k_corr_lines, 
                                                              Epsilon_kimograms_100_lines_MATE_BA]):

    ax.text(6.15/t_p, 1.55e4, letter, fontsize=11)

    for e, c, sim_eps, h_f, h_type, e_c in zip(epsilons, colors, simulated_epsilon, hatch_fill__, histtype__, edgecolor_):

        ax.axvline(np.median(e), ls='-', color=c, lw=lw, alpha=1)
        ax.axvline(sim_eps, ls='--', color=c, lw=lw, alpha=0.75)
        ax.hist(e.ravel(), bins=bins, color= c, density=density, log=log, hatch=h_f, alpha=0.75, edgecolor=e_c, histtype=h_type)

    ax.set_ylim(0, 17500)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(4,4),  useMathText=True)  # Apply scientific notation to y-axis
    ax.yaxis.get_offset_text().set_fontsize(fontsize-6)


## Figure parameters     
plt.xlabel(r'$\epsilon$ (cpms)', fontsize = fontsize)
plt.ylabel('                            Frequency', loc='center', fontsize = fontsize, labelpad = labelpad*4.5)

for ax in axs.flatten():
    ax.tick_params(which='minor', length=1.25, width=1)
    ax.tick_params(which='major', length=2.25, width=1)
    ax.tick_params(axis='both', labelsize=fontsize-3)

## affect borders lines plot making theme thinners
for ax in axs.flatten():    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.6)

axs[2].ticklabel_format(style='scientific', axis='x', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
axs[2].xaxis.get_offset_text().set_fontsize(fontsize-6)


fig.subplots_adjust(top=0.8,
                    bottom=0.145,
                    left=0.15,
                    right=0.925,
                    hspace=0.3,
                    wspace=0.2)

## Figure legend
handles, labels = plt.gca().get_legend_handles_labels()
correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='//////', facecolor='white', edgecolor='k', alpha=0.75)
uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)

empty_line = Line2D([0], [0], ls='', color='k', lw=lw*0.75)
dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)
solid_line = Line2D([0], [0], ls='-', color='k', lw=lw)

handles = [correlated_legend_box, uncorrelated_legend_box, dash_line, solid_line]
labels = [r'Correlated lines', r'Uncorrelated (Interspersed) lines', r'$\epsilon_{sim}$', r'$\epsilon_{median}$']
fig.legend(handles, labels, loc='upper right', fontsize=fontsize-3.5, ncol=2, mode='expand')

plt.show()

