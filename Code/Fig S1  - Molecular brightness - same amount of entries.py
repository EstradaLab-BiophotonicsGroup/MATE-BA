# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:21:45 2024

@author: Ignacio Sallaberry
"""
import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)
    
from Functions import *

#%%
SEED = 2       ## Use SEED = 2 to reproduce Figure S1 or load saved data
# SEED = None  ## Use SEED = None to obtain new random values

points = np.logspace(2, 8, num=100)

N_1e2_index = min(range(len(points)), key=lambda i: abs(points[i]-1e2))
N_1e4_index = min(range(len(points)), key=lambda i: abs(points[i]-1e4))
N_1e6_index = min(range(len(points)), key=lambda i: abs(points[i]-1e6))

#%%
## Pick 10 random values from Epsilon values

Epsilon_distributions_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Epsilon distributions - correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
Epsilon_distributions_UN_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Epsilon distributions - UN correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()

k_random_Epsilon_correlated = []
k_random_Epsilon_UN_correlated = []
k=10
for i in zip([N_1e2_index, N_1e4_index, N_1e6_index]):

    eps_correlated = Epsilon_distributions_correlated_values[i]
    eps_UN_correlated = Epsilon_distributions_UN_correlated_values[i]

    random.seed(SEED)

    number_of_non_nan_elements = np.count_nonzero(~np.isnan(eps_correlated))

    k_random_Epsilon_correlated.append(random.sample(list(eps_correlated[0:number_of_non_nan_elements]),k=k))
    k_random_Epsilon_UN_correlated.append(random.sample(list(eps_UN_correlated[0:number_of_non_nan_elements]),k=k))



k_random_Epsilon_correlated = np.array(k_random_Epsilon_correlated)
k_random_Epsilon_UN_correlated = np.array(k_random_Epsilon_UN_correlated)

#%%
## Load save data

# k_random_Epsilon_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S1 - 10 random values - Epsilons distributions - correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))
# k_random_Epsilon_UN_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S1 - 10 random values - Epsilons distributions - UN correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))

#%% Comparison between CORR (hatched) and UNCORR (solid) distributions

Mann_Whitney_U_test_Fig_S1A = mannwhitneyu(k_random_Epsilon_correlated[0],
                                           k_random_Epsilon_UN_correlated[0], nan_policy='omit')[1]

Mann_Whitney_U_test_Fig_S1B = mannwhitneyu(k_random_Epsilon_correlated[1],
                                           k_random_Epsilon_UN_correlated[1], nan_policy='omit')[1]

Mann_Whitney_U_test_Fig_S1C = mannwhitneyu(k_random_Epsilon_correlated[2],
                                           k_random_Epsilon_UN_correlated[2], nan_policy='omit')[1]


print ('Fig S1 A: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s ' % Mann_Whitney_U_test_Fig_S1A)
print ('Fig S1 B: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s ' % Mann_Whitney_U_test_Fig_S1B)
print ('Fig S1 C: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s ' % Mann_Whitney_U_test_Fig_S1C)

#%%
## Plot Figure S1
## Molecular brightness distribution where each entry is calculated using N=100, 100k and 1M data points

tp = 1e-6 

fig, axs = plt.subplots(3, 1, figsize=(3, 3), sharex=True, dpi=300)

bins = np.arange(-0.5/tp,2/tp, 0.1/tp)
log=False
density = False

letters = ['A', 'B', 'C']

for i, ax, letter_, eps_corr, eps_UN_corr  in zip([0,1,2], axs.flatten(), letters,
                                                  k_random_Epsilon_correlated, 
                                                  k_random_Epsilon_UN_correlated):

    ax.text(2.55/tp, 8.85/tp, letter_, fontsize=11)
    ax.hist(k_random_Epsilon_UN_correlated[i], bins=bins, log=log, density = density, color= 'olivedrab', alpha=0.5, edgecolor='None')
    ax.hist(k_random_Epsilon_correlated[i], bins=bins, log=log, density = density, color= 'white', hatch="//////", alpha=0.5, edgecolor='lightgreen')
    ax.axvline(np.nanmedian(k_random_Epsilon_UN_correlated[i]), ls='-', color='olivedrab', lw=0.75*lw, alpha=0.75)
    ax.axvline(np.nanmedian(k_random_Epsilon_correlated[i]), ls='-', color='lightgreen', lw=lw)
    ax.axvline(0.5/tp, ls='--', color='k', lw=0.75*lw, alpha=0.75)
    

plt.xlabel(r'$\epsilon$ (cpms)', fontsize = fontsize)
plt.ylabel('                                   Frequency', loc='center', fontsize = fontsize, labelpad = labelpad*2)


for ax in axs.flatten():
    ax.tick_params(which='minor', length=1.25, width=1)
    ax.tick_params(which='major', length=2.25, width=1)
    
    ax.tick_params(axis='both', labelsize=fontsize-3)
    ax.tick_params(axis='both', labelsize=fontsize-3)

# ## affect borders lines plot making theme thinners
for ax in axs.flatten():    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.6)

axs[2].ticklabel_format(style='scientific', axis='x', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
axs[2].xaxis.get_offset_text().set_fontsize(fontsize-6)


axs[0].set_xlim(-1.5/tp, 2.5/tp)
axs[0].set_ylim(0, 11.05)
axs[1].set_ylim(0, 11.05)
axs[2].set_ylim(0, 11.05)

fig.subplots_adjust(top=0.825,
                    bottom=0.140,
                    left=0.155,
                    right=0.935,
                    hspace=0.220,
                    wspace=0.2)

## access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='//////', facecolor='white', edgecolor='k', alpha=0.75)
uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)
dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)

## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r'Correlated data', r'Uncorrelated data (shuffled)', r'$\epsilon_{sim}$']
fig.legend(handles, labels, loc='upper right', fontsize=fontsize-3, ncol=2)

plt.show()
