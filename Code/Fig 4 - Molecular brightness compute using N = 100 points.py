# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:04:28 2024

@author: Ignacio Sallaberry
"""

import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)

from Functions import *

#%%

intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 1 - 2 and 4\Simulation Fig 1 - 2 and 4.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()


tp=1e-6
Epsilon_simulated=0.5/tp

## N = number_of_points_used_in_epsilon_calculation
N = 100

minumum_number_of_epsilons_per_distribution = 100

## Define the time interval distance "J" values
Num=50 ## Number of J values to evaluate
J = np.array(list(np.logspace(1, 4, num=Num, dtype=int))+
                             list(np.logspace(4.1, math.log10(len(intensity_values)/N), num=Num, dtype=int)))
J[0] = 1

#%%
MEDIAN_ = []
MAD_ = []
EPS_ = []

## Auxiliar list for analysis that will contain epsilon values for a given value of J
## From this auxiliar list medians, MAD and 
epsilons_ = []

for u in tqdm.trange(0,len(J)):    
    ## j is the temporal spacing value
    j = J[u]


    ## index_ is a set of list with the point's temporally spaced coordinates of data points.
    ## Each list in index will be use in the calculation of epsilon. if len(index_)=100, then 100 values of epsilon will be calculated
    index_ = []
    for i in range (0, int(len(intensity_values)/(j*N))):
        index_.append(list(np.arange(0, j*N, j) + j*N * i))

    ##------------------------------
    if len(index_)>=2:
    
        ## Complete number of index packages until complete the minumum_number_of_epsilons_per_distribution value
        if len(index_)<minumum_number_of_epsilons_per_distribution and len(index_)>=2:

            ## ind_ is the last pivot used
            ind_ = index_[-1][0] 
            
            ## Here will construct complementary index list until index's lenght is minumum_number_of_epsilons_per_distribution value
            while len(index_) < minumum_number_of_epsilons_per_distribution:
                    
                ## random_pivot is a possible first index coordinate value within 0 and the las pivot used (ind_)
                random_pivot = random.sample(list(np.arange(0,ind_)), k=1)[0]
                
                ## Make sure that the random pivot value has not been used as index value or to construct complementary index list
                while index_[0].count(random_pivot):
                    ## If random pivot is used, then a new random value is selected
                    random_pivot = random.sample(list(np.arange(1,ind_)), k=1)[0]
                
                ## Defined the pivot value a new list of point's temporally spaced coordinates of data points spaced in j is construct
                complementary_index_ = list(np.arange(random_pivot,
                                                      random_pivot+j*N,
                                                      j))
                
                index_.append(complementary_index_)
        
    ##------------------------------
    if len(index_)<2:
        
        ## Complete number of index packages until complete the minumum_number_of_epsilons_per_distribution value
        while len(index_)<minumum_number_of_epsilons_per_distribution:
            random_pivot = random.sample(list(np.arange(1,index_[0][1]-1)), k=1)[0]

            while index_[0].count(random_pivot):
                random_pivot = random.sample(list(np.arange(1,index_[0][1]-1)), k=1)[0]
            
            complementary_index_ = list(np.arange(random_pivot,
                                                  random_pivot+j*N,
                                                  j))

            index_.append(complementary_index_)


    ## Construct an array only with the data that will be used in the calculation  
    data_ = np.array([intensity_values[ind] for ind in list(np.array(index_).ravel())])

    ## Calculate molecular brightness value (epsilon). 
    ## The total number of epsilons will be len(data_)//N
    ## The minimum number of epsilons will be minumum_number_of_epsilons_per_distribution
    epsilons_ = (B_from_n_points(data_, N)[2] -1)/0.3536/tp

    ## Construct a vector of with the median of the epsilons distribution
    MEDIAN_.append(np.nanmedian(epsilons_))

    ## Construct a vector of with the std of the epsilons distribution
    MAD_.append(median_abs_deviation(epsilons_,nan_policy='omit'))

    ## Attach the apsilons values to a vector.
    EPS_.append(epsilons_.ravel())


#%%

MEDIAN_ = np.array(MEDIAN_)
MAD_ = np.array(MAD_)

#%%
## Load save data to reproduce Figure 4

#saved_data = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 4 - Epsilon values.txt', delimiter='\t')
#EPS_ = list(saved_data.apply(lambda x: np.array((x.dropna()).to_list()), axis=1))

#saved_data = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 4 - Median values.txt', delimiter='\t')
#MEDIAN_ = np.array(list(saved_data.apply(lambda x: np.array(x.dropna().to_list()), axis=1))).ravel()

#saved_data = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 4 - MAD values.txt', delimiter='\t')
#MAD_ = np.array(list(saved_data.apply(lambda x: np.array(x.dropna().to_list()), axis=1))).ravel()


#%%

## Look the position in MAD list where the difference between median of molecular brightness distribution and simulated epsilon value is 5%
d1=0
delta_diferences_index = 0
while d1 < len(MEDIAN_):
    if abs(MEDIAN_[d1]-Epsilon_simulated)<0.01*max(abs(np.array(MEDIAN_)-Epsilon_simulated)):
        delta_diferences_index=d1
        d1=len(MEDIAN_)
    d1+=1


#%%

fig, ax0 = plt.subplots(1, 1, figsize=(6,3), sharex=True, dpi=300)

ax0.axhline(Epsilon_simulated, 0, J[-1], ls='--', color='k', lw=lw*0.75, alpha=0.75)

ax0.semilogx(J[0:len(EPS_)], MEDIAN_, '-', color='r')

ax0.fill_between(J[0:len(EPS_)],
                 np.array(MEDIAN_) - np.array(MAD_),
                 np.array(MEDIAN_) + np.array(MAD_),
                 color= 'none', hatch="//////", alpha=0.5, edgecolor='r')

## Arrows parameters
head_length=0.4
head_width=0.5

## Top arrow
arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='r', edgecolor='r', lw=1.5)
ax0.annotate('', xy=(J[delta_diferences_index], (MEDIAN_ + MAD_)[delta_diferences_index]*1.05),
              xytext=(J[delta_diferences_index], (MEDIAN_ + MAD_)[delta_diferences_index]*1.5),
              arrowprops=arrowprops)

## Bottom arrow
arrowprops = dict(arrowstyle='->, head_length=%s, head_width=%s'%(head_length, head_width), facecolor='r', edgecolor='r', lw=1.5)
ax0.annotate('', xy=(J[delta_diferences_index], 0.009),
              xytext=(J[delta_diferences_index], (MEDIAN_ - MAD_)[delta_diferences_index]*-1.85),
              arrowprops=arrowprops)

ax0.set_ylabel(r'$\epsilon_{median}$', fontsize = fontsize, labelpad = labelpad)
ax0.set_xlabel('$J$', fontsize = fontsize, labelpad = labelpad)

ax0.tick_params(which='minor', length=1.25, width=1)
ax0.tick_params(which='major', length=2.25, width=1)
ax0.tick_params(axis='both', labelsize=fontsize-3)

ax0.set_ylim(-0.4/tp, 1.7/tp)

ax0.ticklabel_format(style='scientific', axis='y', scilimits=(6,6),  useMathText=True)  # Apply scientific notation to y-axis
ax0.yaxis.get_offset_text().set_fontsize(fontsize-6)

fig.subplots_adjust(top=0.940,
                    bottom=0.155,
                    left=0.115,
                    right=0.975,
                    hspace=0.2,
                    wspace=0.2)

## affect borders lines plot making theme thinners
for axis in ['top','bottom','left','right']:
    ax0.spines[axis].set_linewidth(0.6)
    
plt.show()


