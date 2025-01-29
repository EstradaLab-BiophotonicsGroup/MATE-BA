# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:14:56 2024

@author: Ignacio Sallaberry
"""

import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)

from Functions import *

#%%
SEED = 0       ## Use SEED = 0 to reproduce Figure S2 or load saved data
# SEED = None  ## Use SEED = None to obtain new random values

tp = 1e-6
Epsilon_simulated = 500000

N = np.logspace(2, 8, num=100)

N_1e2_index = min(range(len(N)), key=lambda i: abs(N[i]-1e2))
N_1e4_index = min(range(len(N)), key=lambda i: abs(N[i]-1e4))
N_1e6_index = min(range(len(N)), key=lambda i: abs(N[i]-1e6))

#%%

# intensity_points = read_B64(path[:-5]+r'\Simulations\Fig 1 and 2 - Simulations\Simulation Fig 1 and 2.b64')
intensity_points = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 1 - 2 and 4\Simulation Fig 1 - 2 and 4.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
points_random = Random_data(intensity_points, seed=SEED)    

#%%
## Get distributions of Mean Intensity, Variance and Epsilon in cpms

Mean_I_distributions_correlated_values, Variance_distributions_correlated_values, Epsilon_distributions_correlated_values = [d/tp for d in Get_distributions(intensity_points, N)]

Mean_I_distributions_UN_correlated_values, Variance_distributions_UN_correlated_values, Epsilon_distributions_UN_correlated_values = [d/tp for d in Get_distributions(points_random, N)]

#%%
k_random_I_mean_correlated = []
k_random_VAR_correlated = []
k_random_Epsilon_correlated = []
k=100
for k_i, k_v, k_b in zip(Mean_I_distributions_correlated_values,
                         Variance_distributions_correlated_values,
                         Epsilon_distributions_correlated_values):
        
        random.seed(SEED)
    
        number_of_non_nan_elements = np.count_nonzero(~np.isnan(k_b))
    
        if number_of_non_nan_elements >=100:
            k_random_I_mean_correlated.append(random.sample(list(k_i[0:number_of_non_nan_elements]),k=k))
            k_random_VAR_correlated.append(random.sample(list(k_v[0:number_of_non_nan_elements]),k=k))
            k_random_Epsilon_correlated.append(random.sample(list(k_b[0:number_of_non_nan_elements]),k=k))
        
        else:
            break
k_random_I_mean_correlated = np.array(k_random_I_mean_correlated)
k_random_VAR_correlated = np.array(k_random_VAR_correlated)
k_random_Epsilon_correlated = np.array(k_random_Epsilon_correlated)


#%%
k=100
k_random_I_mean_UN_correlated = []
k_random_VAR_UN_correlated = []
k_random_Epsilon_UN_correlated = []
for k_i, k_v, k_b in zip(Mean_I_distributions_UN_correlated_values,
                         Variance_distributions_UN_correlated_values,
                         Epsilon_distributions_UN_correlated_values):
        
        random.seed(SEED)    
    
        number_of_non_nan_elements = np.count_nonzero(~np.isnan(k_b))
    
        if number_of_non_nan_elements >=100:
            k_random_I_mean_UN_correlated.append(random.sample(list(k_i[0:number_of_non_nan_elements]),k=k))
            k_random_VAR_UN_correlated.append(random.sample(list(k_v[0:number_of_non_nan_elements]),k=k))
            k_random_Epsilon_UN_correlated.append(random.sample(list(k_b[0:number_of_non_nan_elements]),k=k))
        else:
            break

k_random_I_mean_UN_correlated = np.array(k_random_I_mean_UN_correlated)
k_random_VAR_UN_correlated = np.array(k_random_VAR_UN_correlated)
k_random_Epsilon_UN_correlated = np.array(k_random_Epsilon_UN_correlated)

#%%
# Load save data

# k_random_I_mean_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - I mean distributions - correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))
# k_random_VAR_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - Variance distributions - correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))
# k_random_Epsilon_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - Epsilon distributions - correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))

# k_random_I_mean_UN_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - I mean distributions - UN correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))
# k_random_VAR_UN_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - Variance distributions - UN correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))
# k_random_Epsilon_UN_correlated = np.array(list(pd.read_csv(path[:-5]+r'\Data for Figures\Fig S2 - k random values - Epsilon distributions - UN correlated data.txt', delimiter='\t').apply(lambda x: np.array(x.to_list()), axis=1)))

#%%
plt.close('all')

MEDIAN_k_random_I_mean_correlacionados = np.nanmedian(k_random_I_mean_correlated, axis=1)
MAD_k_random_I_mean_correlacionados = median_abs_deviation(k_random_I_mean_correlated, axis=1, nan_policy='omit')  
MEDIAN_k_random_I_mean_NO_correlacionados = np.nanmedian(k_random_I_mean_UN_correlated, axis=1)
MAD_k_random_I_mean_NO_correlacionados = median_abs_deviation(k_random_I_mean_UN_correlated, axis=1, nan_policy='omit')

MEDIAN_k_random_VAR_correlacionados = np.nanmedian(k_random_VAR_correlated, axis=1)
MAD_k_random_VAR_correlacionados = median_abs_deviation(k_random_VAR_correlated, axis=1, nan_policy='omit')
MEDIAN_k_random_VAR_NO_correlacionados = np.nanmedian(k_random_VAR_UN_correlated, axis=1)
MAD_k_random_VAR_NO_correlacionados = median_abs_deviation(k_random_VAR_UN_correlated, axis=1, nan_policy='omit')  

MEDIAN_k_random_epsilon_correlacionados = np.nanmedian(k_random_Epsilon_correlated, axis=1)
MAD_k_random_epsilon_correlacionados = median_abs_deviation(k_random_Epsilon_correlated, axis=1, nan_policy='omit')  
MEDIAN_k_random_epsilon_NO_correlacionados = np.nanmedian(k_random_Epsilon_UN_correlated, axis=1)
MAD_k_random_epsilon_NO_correlacionados = median_abs_deviation(k_random_Epsilon_UN_correlated, axis=1, nan_policy='omit')

#%%
## Plot Figure S2
## Median of distributions as functions of number of data point use to calculate each entry

N_ = np.array(N[0:len(MEDIAN_k_random_I_mean_correlacionados)])

fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(3,3), sharex=True, dpi=300)

letters = ['A', 'B', 'C']
ax0.text(1.45e7, 2.4/tp, letters[0], fontsize=11)
ax1.text(1.45e7, 2.4/tp, letters[1], fontsize=11)
ax2.text(1.45e7, 0.8/tp, letters[2], fontsize=11)


ax0.semilogx(N_, MEDIAN_k_random_I_mean_correlacionados, '-', color='#0267C1', lw=lw)
ax0.plot(N_, MEDIAN_k_random_I_mean_NO_correlacionados, '-', color='#2A0C4E', lw=lw)
ax0.fill_between(N_,
                  MEDIAN_k_random_I_mean_correlacionados - MAD_k_random_I_mean_correlacionados, 
                  MEDIAN_k_random_I_mean_correlacionados + MAD_k_random_I_mean_correlacionados, 
                  color= 'none', hatch="//////", alpha=0.5, edgecolor='#0267C1', label= 'I correlated data')
ax0.fill_between(N_,
                  MEDIAN_k_random_I_mean_NO_correlacionados - MAD_k_random_I_mean_NO_correlacionados,
                  MEDIAN_k_random_I_mean_NO_correlacionados + MAD_k_random_I_mean_NO_correlacionados, 
                  color='#2A0C4E', alpha=0.3, label ='I uncorrelated data')


ax1.semilogx(N_, MEDIAN_k_random_VAR_correlacionados, '-', color='salmon', lw=lw)
ax1.plot(N_, MEDIAN_k_random_VAR_NO_correlacionados, '-', color='goldenrod', lw=lw)
ax1.fill_between(N_,
                  MEDIAN_k_random_VAR_correlacionados - MAD_k_random_VAR_correlacionados, 
                  MEDIAN_k_random_VAR_correlacionados + MAD_k_random_VAR_correlacionados, 
                  color= 'none', hatch="//////", alpha=0.5, edgecolor='salmon', label= 'VAR correlated data')
ax1.fill_between(N_,
                  MEDIAN_k_random_VAR_NO_correlacionados - MAD_k_random_VAR_NO_correlacionados,
                  MEDIAN_k_random_VAR_NO_correlacionados + MAD_k_random_VAR_NO_correlacionados, 
                  color='goldenrod', alpha=0.3, label ='VAR uncorrelated data')


ax2.semilogx(N_, MEDIAN_k_random_epsilon_correlacionados, '-', color='lightgreen', lw=lw)
ax2.plot(N_, MEDIAN_k_random_epsilon_NO_correlacionados, '-', color='olivedrab', lw=lw)
ax2.fill_between(N_,
                  MEDIAN_k_random_epsilon_correlacionados - MAD_k_random_epsilon_correlacionados, 
                  MEDIAN_k_random_epsilon_correlacionados + MAD_k_random_epsilon_correlacionados, 
                  color= 'none', hatch="//////", alpha=0.5, edgecolor='lightgreen', label= r'$\epsilon$ correlated data')
ax2.fill_between(N_,
                  MEDIAN_k_random_epsilon_NO_correlacionados - MAD_k_random_epsilon_NO_correlacionados,
                  MEDIAN_k_random_epsilon_NO_correlacionados + MAD_k_random_epsilon_NO_correlacionados, 
                  color='olivedrab', alpha=0.3, label =r'$\epsilon$ correlated data')

ax2.axhline(Epsilon_simulated, 0, N[-1]*1.25, ls='--', color='k', lw=lw*0.75, alpha=0.75)


ax0.set_ylabel(r'$<I>_{median}$', fontsize = fontsize-2, labelpad = labelpad)
ax1.set_ylabel(r'$\sigma^2_{median}$', fontsize = fontsize-2, labelpad = labelpad)
ax2.set_ylabel(r'$\epsilon_{median}$', fontsize = fontsize-2, labelpad = labelpad)
plt.xlabel('$N$', fontsize = fontsize, labelpad = labelpad)


ax0.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
ax0.yaxis.get_offset_text().set_fontsize(fontsize-6)

ax1.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
ax1.yaxis.get_offset_text().set_fontsize(fontsize-6)

ax2.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
ax2.yaxis.get_offset_text().set_fontsize(fontsize-6)


ax0.tick_params(which='minor', length=1.25, width=1)
ax0.tick_params(which='major', length=2.25, width=1)
ax1.tick_params(which='minor', length=1.25, width=1)
ax1.tick_params(which='major', length=2.25, width=1)
ax2.tick_params(which='minor', length=1.25, width=1)
ax2.tick_params(which='major', length=2.25, width=1)

ax0.tick_params(axis='both', labelsize=fontsize-3)
ax1.tick_params(axis='both', labelsize=fontsize-3)
ax2.tick_params(axis='both', labelsize=fontsize-3)

ax0.set_ylim(1.3/tp, 2.55/tp)
ax1.set_ylim(1.3/tp, 2.55/tp)
ax2.set_ylim(-0.45/tp, 0.95/tp)
ax0.set_xlim(57, 1.3e7)


## affect borders lines plot making it thinners
for axis in ['top','bottom','left','right']:
    ax0.spines[axis].set_linewidth(0.6)
    ax1.spines[axis].set_linewidth(0.6)
    ax2.spines[axis].set_linewidth(0.6)

fig.subplots_adjust(top=0.795,
                    bottom=0.135,
                    left=0.165,
                    right=0.935,
                    hspace=0.355,
                    wspace=0.2)

## access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='//////', facecolor='white', edgecolor='k', alpha=0.75)
uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)
dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)

## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r'Correlated data', r'Uncorrelated data (Mixed)', r'$\epsilon_{sim}$']
fig.legend(handles, labels, loc='upper right', fontsize=fontsize-3, ncol=2)

plt.show()

