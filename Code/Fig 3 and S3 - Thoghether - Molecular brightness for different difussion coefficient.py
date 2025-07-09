# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 20:53:52 2024

@author: Ignacio Sallaberry
"""
import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)

from Functions import *

#%%
SEED = 0       ## Use SEED = 0 to reproduce Figure 3
# SEED = None  ## Use SEED = None to obtain new random values

## Amount of data points used to calculate each entry
N = np.logspace(2, 8, num=100)

tp = 1e-6

#%%

point_D_0_5 = np.concatenate((pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_0_5_a.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy(),
                              pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_0_5_b.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy(),
                              pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_0_5_c.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()))

point_D_1 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 1 - 2 and 4\Simulation Fig 1 - 2 and 4.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_10 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_10.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_20 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_20.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_50 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_50.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_100 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_100.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_200 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_200.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
point_D_400 = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 3\Simulation Fig 3 - D_400.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()


#%%

points_random_D_0_5 = Random_data(point_D_0_5, seed=SEED)
points_random_D_1 = Random_data(point_D_1, seed=SEED)
points_random_D_10 = Random_data(point_D_10, seed=SEED)
points_random_D_20 = Random_data(point_D_20, seed=SEED)
points_random_D_50 = Random_data(point_D_50, seed=SEED)
points_random_D_100 = Random_data(point_D_100, seed=SEED)
points_random_D_200 = Random_data(point_D_200, seed=SEED)
points_random_D_400 = Random_data(point_D_400, seed=SEED)


#%%
Epsilon_distributions_correlated_values_for_each_D = [d/tp for d in Get_distributions_2 ([point_D_0_5, point_D_1, point_D_10, point_D_20, point_D_50, point_D_100, point_D_200, point_D_400], N)]

Epsilon_distributions_UN_correlated_values_for_each_D = [d/tp for d in Get_distributions_2 ([points_random_D_0_5, points_random_D_1, points_random_D_10, points_random_D_20, points_random_D_50, points_random_D_100, points_random_D_200, points_random_D_400], N)]

#%%

del point_D_0_5
del point_D_1
del point_D_10
del point_D_20
del point_D_50
del point_D_100
del point_D_200
del point_D_400

del points_random_D_0_5
del points_random_D_1
del points_random_D_10
del points_random_D_20
del points_random_D_50
del points_random_D_100
del points_random_D_200
del points_random_D_400


#%%
### Prepared vaiables only for ploting FIG 3
X_corr = np.array([Epsilon_distributions_correlated_values_for_each_D[0],
                   Epsilon_distributions_correlated_values_for_each_D[2],
                   Epsilon_distributions_correlated_values_for_each_D[5]])


X_UN_corr = np.array([Epsilon_distributions_UN_correlated_values_for_each_D[0],
                      Epsilon_distributions_UN_correlated_values_for_each_D[2],
                      Epsilon_distributions_UN_correlated_values_for_each_D[5]])

#%%
## Plot Fig 3
Epsilon_sim = 0.5/tp
colors = [   '#F51AA4', '#F3B61F', '#0A81D1']    
D = np.asarray([0.5,       10,        100]) ## um**2/seg

plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(3, 3), sharex=True, dpi=300)

letters = ['A', 'B', 'C']
for ax, letter in zip (axs.flatten(), letters):
    ax.text(2.55e8, 0.85/tp, letter, fontsize=11)

for i, ax, x_corr, x_UN_corr, c, d in zip([0,1,2], axs.flatten(), X_corr, X_UN_corr, colors, D):

    ## This is just for plotting reasons
    ## Because D=0.1um2/s simulations is longer than others difussion coefficient simulations, variables are calculated appart. 
    if i==0:
        MEDIAN_x_corr = np.nanmedian(x_corr, axis=1)
        MAD_x_corr = median_abs_deviation(x_corr, axis=1, nan_policy='omit')   
        MEDIAN_x_UN_corr = np.nanmedian(x_UN_corr, axis=1)
        MAD_x_UN_corr = median_abs_deviation(x_UN_corr, axis=1, nan_policy='omit')  
        points_ = N
        
    else:        
        MEDIAN_x_corr = np.nanmedian(x_corr[0:-15], axis=1)
        MAD_x_corr = median_abs_deviation(x_corr[0:-15], axis=1, nan_policy='omit')   
        MEDIAN_x_UN_corr = np.nanmedian(x_UN_corr[0:-15], axis=1)
        MAD_x_UN_corr = median_abs_deviation(x_UN_corr[0:-15], axis=1, nan_policy='omit')   
        points_ = N[0:-15]

    ax.semilogx(points_, MEDIAN_x_corr, '-', color='%s'%c, lw=lw)
    ax.fill_between(points_,
                      MEDIAN_x_corr - MAD_x_corr,
                      MEDIAN_x_corr + MAD_x_corr,
                      color= 'none', hatch="//////", alpha=0.5, edgecolor='%s'%c, label=r'correlated')

    ax.semilogx(points_, MEDIAN_x_UN_corr, '-', color='%s'%c, lw=lw)
    ax.fill_between(points_,
                      MEDIAN_x_UN_corr - MAD_x_UN_corr,
                      MEDIAN_x_UN_corr + MAD_x_UN_corr,
                      color='%s'%c, alpha=0.3, label=r'mixed - data')

    ax.axhline(Epsilon_sim, lw=lw, ls='--', color='%s'%(c))

    if d>1:
        d=int(d)
    ax.annotate(text=r'$D=$%s' % d +r' $\frac{\mu m^{2}}{sec}$', xy=(9.5e5, -0.25/tp),fontsize=fontsize-3, color='%s' % c)
    

    ## Locate value of N when difference between the median of molecular brightness distribution from correlated data and the simulated value of epsilon is equal to 5%
    d1=0
    delta_diferences_index = 0
    while d1 < len(x_corr):
        if abs(MEDIAN_x_corr[d1]-Epsilon_sim)<0.050*max(abs(MEDIAN_x_corr-Epsilon_sim)):  ## Me quedo con el indice cuando es (5 +- 0.5) %
            delta_diferences_index=d1
            d1=len(x_corr)
        d1+=1

    ## Locate value of N for molecular brightness distribution from UNcorrelated data present the same MAD as the correlated epsilon distribution whichs median presents 5% to simulated epsilon value    
    delta_diferences_UN_corr_index = min(range(len(MAD_x_UN_corr)), key=lambda i: abs(MAD_x_UN_corr[i]-MAD_x_corr [delta_diferences_index]))    
    
    print(delta_diferences_UN_corr_index, delta_diferences_index)
    print(MAD_x_UN_corr[delta_diferences_UN_corr_index], MAD_x_corr[delta_diferences_index])

    head_length=0.2
    head_width=0.05

    # Top right arrow
    arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_index], (MEDIAN_x_corr + MAD_x_corr)[delta_diferences_index]),
                  xytext=(points_[delta_diferences_index], (MEDIAN_x_corr + MAD_x_corr)[delta_diferences_index]+0.475/tp),
                  arrowprops=arrowprops)

    # Bottom right arrow    
    arrowprops = dict(arrowstyle='<-, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_index], (MEDIAN_x_corr - MAD_x_corr)[delta_diferences_index]-0.475/tp),
                  xytext=(points_[delta_diferences_index], (MEDIAN_x_corr - MAD_x_corr)[delta_diferences_index]),
                  arrowprops=arrowprops)
    
    # Top left arrow
    arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr + MAD_x_UN_corr)[delta_diferences_index]),
                  xytext=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr + MAD_x_UN_corr)[delta_diferences_index]+0.475/tp),
                  arrowprops=arrowprops)
    
    # Bottom left arrow
    arrowprops = dict(arrowstyle='<-, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr- MAD_x_UN_corr)[delta_diferences_index]-0.475/tp),
                  xytext=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr - MAD_x_UN_corr)[delta_diferences_index]),
                  arrowprops=arrowprops)


    ax.set_ylim(-0.4/tp, 0.99/tp)

    ax.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
    ax.yaxis.get_offset_text().set_fontsize(fontsize-5)


fig.supylabel(r'$\epsilon$ $_{median}$', fontsize = fontsize)
fig.supxlabel('$N$', fontsize = fontsize)


for ax in axs.flatten():
    ax.tick_params(which='minor', length=1.25, width=1)
    ax.tick_params(which='major', length=2.25, width=1)
    
    ax.tick_params(axis='both', labelsize=fontsize-3)
    ax.tick_params(axis='both', labelsize=fontsize-3)

## affect borders lines plot making theme thinners
for ax in axs.flatten():    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.6)

fig.subplots_adjust(top=0.79,
                    bottom=0.14,
                    left=0.175,
                    right=0.9,
                    hspace=0.33,
                    wspace=0.2)


correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='//////', facecolor='white', edgecolor='k', alpha=0.75)
uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)
dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)

## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r'Correlated data', r'Uncorrelated data (shuffled)', r'$\epsilon_{sim}$']
fig.legend(handles, labels, loc='upper center', fontsize=fontsize-3, ncol=2)

plt.show()

#%%
## Plot Fig S3

plt.close('all')

D = np.asarray([0.5,      1,         10,         20,        50,       100,       200,       400]) ## um**2/seg
colors = [  '#F51AA4', '#264000', '#F3B61F', '#FF4242', '#FB8B24', '#0A81D1', '#912F56', '#FF69EB']
letters = [    'A',       'B',       'C',       'D',        'E',      'F',       'G',       'H']

## https://coolors.co/efaac4-ffc4d1-ffe8e1-568259-912f56

fig, axs = plt.subplots(4, 2, figsize=(6, 3), sharex=True, dpi=300)

for i, ax in zip([0,4,1,5,2,6,3,7], axs.flatten()):
    ax.text(2.55e8, 0.8/tp, letters[i], fontsize=11)
    x_corr = Epsilon_distributions_correlated_values_for_each_D[i]
    x_UN_corr = Epsilon_distributions_UN_correlated_values_for_each_D[i]
    c = colors[i]
    d = D[i]
    
    MEDIAN_x_corr = np.nanmedian(x_corr, axis=1)
    MAD_x_corr  = median_abs_deviation(x_corr, axis=1, nan_policy='omit') 
    MEDIAN_x_UN_corr = np.nanmedian(x_UN_corr, axis=1)
    MAD_x_UN_corr = median_abs_deviation(x_UN_corr, axis=1, nan_policy='omit') 
    points_ = N
    
    
    ax.semilogx(points_, MEDIAN_x_corr, '-', color='%s'%c, lw=lw)
    ax.fill_between(points_,
                      MEDIAN_x_corr - MAD_x_corr ,
                      MEDIAN_x_corr + MAD_x_corr ,
                      color= 'none', hatch="//////", alpha=0.5, edgecolor='%s'%c, label=r'correlated')

    ax.semilogx(points_, MEDIAN_x_UN_corr, '-', color='%s'%c, lw=lw)
    ax.fill_between(points_,
                      MEDIAN_x_UN_corr - MAD_x_UN_corr,
                      MEDIAN_x_UN_corr + MAD_x_UN_corr,
                      color='%s'%c, alpha=0.3, label=r'mixed - data')

    ax.axhline(Epsilon_sim, lw=lw, ls='--', color='%s'%(c))

    if d>1:
        d=int(d)
    ax.annotate(text=r'$D=$%s' % d +r' $\frac{\mu m^{2}}{sec}$', xy=(1.7e6, -0.35/tp),fontsize=fontsize-3, color='%s' % c)

    
    ## Locate value of N when difference between the median of molecular brightness distribution from correlated data and the simulated value of epsilon is equal to 5%
    d1=0
    delta_diferences_index = 0
    while d1 < len(x_corr):
        if abs(MEDIAN_x_corr[d1]-Epsilon_sim)<0.050*max(abs(MEDIAN_x_corr-Epsilon_sim)):
            delta_diferences_index=d1
            d1=len(x_corr)
        d1+=1

    ## Locate value of N for molecular brightness distribution from UNcorrelated data present the same MAD as the correlated epsilon distribution whichs median presents 5% to simulated epsilon value    
    delta_diferences_UN_corr_index = min(range(len(MAD_x_UN_corr)), key=lambda i: abs(MAD_x_UN_corr[i]-MAD_x_corr [delta_diferences_index]))    
    
    print(delta_diferences_UN_corr_index, delta_diferences_index)
    print(MAD_x_UN_corr[delta_diferences_UN_corr_index], MAD_x_corr [delta_diferences_index])

    ## Head arrow parameters
    head_length=0.2
    head_width=0.075

    # Top right arrow
    arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_index], (MEDIAN_x_corr + MAD_x_corr )[delta_diferences_index]),
                  xytext=(points_[delta_diferences_index], (MEDIAN_x_corr + MAD_x_corr )[delta_diferences_index]+0.40/tp),
                  arrowprops=arrowprops)

    # Bottom right arrow
    arrowprops = dict(arrowstyle='<-, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_index], (MEDIAN_x_corr - MAD_x_corr )[delta_diferences_index]-0.40/tp),
                  xytext=(points_[delta_diferences_index], (MEDIAN_x_corr - MAD_x_corr )[delta_diferences_index]),
                  arrowprops=arrowprops)
    
    # Top left arrow
    arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr + MAD_x_UN_corr)[delta_diferences_index]),
                  xytext=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr + MAD_x_UN_corr)[delta_diferences_index]+0.40/tp),
                  arrowprops=arrowprops)
    
    # Bottom left arrow
    arrowprops = dict(arrowstyle='<-, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='%s' % c, edgecolor='%s' % c)
    ax.annotate('', xy=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr- MAD_x_UN_corr)[delta_diferences_index]-0.40/tp),
                  xytext=(points_[delta_diferences_UN_corr_index], (MEDIAN_x_UN_corr - MAD_x_UN_corr)[delta_diferences_index]),
                  arrowprops=arrowprops)

    
    ax.set_ylim(-0.48/tp, 0.9/tp)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
    ax.yaxis.get_offset_text().set_fontsize(fontsize-5)


fig.supylabel(r'$\epsilon$ $_{median}$', fontsize = fontsize, x=0.00001)
fig.supxlabel('$N$', fontsize = fontsize)


for ax in axs.flatten():
    ax.tick_params(which='minor', length=1.25, width=1)
    ax.tick_params(which='major', length=2.25, width=1)
    ax.tick_params(axis='both', labelsize=fontsize-3)


## affect borders lines plot making theme thinners
for ax in axs.flatten():    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.6)

fig.subplots_adjust(top=0.845,
                    bottom=0.125,
                    left=0.08,
                    right=0.955,
                    hspace=0.39,
                    wspace=0.255)

# fig.suptitle(r'$\epsilon: Molecular Brightness$', y=0.90, fontsize=fontsize)
correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='//////', facecolor='white', edgecolor='k', alpha=0.75)
uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)
dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)

## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r'Correlated data', r'Uncorrelated data (shuffled)', r'$\epsilon_{sim}$']
fig.legend(handles, labels, loc='upper center', fontsize=fontsize-3, ncol=3)

plt.show()




