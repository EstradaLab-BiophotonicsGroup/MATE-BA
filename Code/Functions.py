# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 09:03:14 2024

@author: Ignacio Sallaberry
"""
import numpy as np

import random

import math

import tqdm

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import axis
import matplotlib.patches as mpatches

from scipy.stats import mannwhitneyu 
from scipy.stats import norm
from scipy.stats import median_abs_deviation

import addcopyfighandler  ##allows to Ctrl+C matplotlib figures

import csv

# plt.rcParams["mathtext.fontset"] = "dejavusans" ## Default appearence of matplotlib typography
plt.rcParams["mathtext.fontset"] = "cm" ## Latex tipographic style
plt.rcParams['font.family'] = 'Arial'


fontsize=12
point_size = 4
lw=1
markeredgewidth=0.5  #ancho del borde del simbolo
labelpad = 2


#%%
def Open_figure_data(file_path):
    '''
    
    Parameters
    ----------
    file_path : str
        path to csv file.

    Returns
    -------
    nd-array
        array with intensity values from path file.

    '''

    data = []    
    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file)
        for row in reader:
            # Convert each value in the row to float or keep as np.nan if empty
            processed_row = [float(value) if value else np.nan for value in row]
            data.append(processed_row)

    return np.array(data)


def read_B64(File,Voltear=False):
    '''
    
    Parameters
    ----------
    File : str
        path to .b64 file.

    Returns
    -------
    TYPE 1D array
        numpy array with intensity values.

    '''

    Read=lfd(File)
    Matriz=Read.asarray()
    
    return Matriz


def Random_data (data, seed=None):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    seed : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    data_ : TYPE
        DESCRIPTION.

    '''
    ## If "seed=None" points will be randomlly mixed with a random seed. 
    ## Everytime "Random_data" is run, then a different distribution is obtained.
    rng = np.random.default_rng(seed)
    data_ = data.copy() ## To avoid modify original data, shuffle module must be apply to a copy of data values
    rng.shuffle(data_)    

    ## Obs: random.sample needs works with list. So if data is an array, random.sample needs to convert the array into list. 
    ## Also use memory to mantein the original data and create a new list with random sampled values.
    ## The function finally use memory to convert random sampled list into array
    
    ## Instead numpy.random accept an array as input and... 
    return data_


def Kimogram(Lines, pixel):
    '''
    
    Parameters
    ----------
    Lines : int16 1D-array.
        Lines is an 1D-array of consecutives lines.
    pixel : int
        number of pixel in a line.

    Returns
    -------
    kimogram : int16 2D-array
        Kimogram of lines.
        Each column is a pixel
        Each row is a line

    '''
    
    if type(len(Lines)/pixel) == float:   ## avoid possible uncompleted lines
        last_line = int(len(Lines)/pixel)*pixel
        Lines = Lines[0:last_line]                      
        
    kimogram = (Lines).reshape(int(len(Lines)/pixel),pixel)                                 
    
    return kimogram



#%%
def B_from_n_points(data, points):
    '''
    
    Parameters
    ----------
    data : 1D float array
        Array of Intensity values.
    points : int
        Number of data points use to each entry of a distribution.

    Returns
    -------
    List of floats
        Mean intensity, variance and apparent brightness distributions.

    '''
    I_mean = []
    VAR = []
    
    
    for i in tqdm.trange(0, len(data)-points, points):
        d = data[i:i+points]
        I_mean.append(np.mean(data[i:i+points]))
        VAR.append(abs(np.mean(data[i:i+points]**2)-(np.mean(data[i:i+points])**2)))

    B = np.zeros_like(I_mean)
    for k in range(0,len(I_mean)):
        if I_mean[k] == 0 or VAR[k]==I_mean[k]:
            B[k]=1  
        else:
            B[k] = VAR[k]/I_mean[k]

    return I_mean, VAR, np.asarray(B)


def Get_distributions(data_, N_, Vol_PSF=0.3536):
    '''
    
    Parameters
    ----------
    data_ : 1D float array
        Array of Intensity values.
    N_ : 1D int array
        Each value is the number of data points use to calculate each entries of a distribution.

    Returns
    -------
    I_distributions : 2D array
        Mean intensity distributions.
    Variance_distributions : 2D array
        Variance distributions.
    Epsilon_distributions : 2D array
        Molecular brightness distributions.

    '''

    I_distributions = []
    Variance_distributions = []
    B_distributions = []
    
    
    for p in N_:
        I_mean_ , VAR_, B_ = B_from_n_points(data_, int(p))
        
        I_distributions.append(list(I_mean_))
        Variance_distributions.append(list(VAR_))
        B_distributions.append(list(B_))
        
    ## To convert list into arrays fullfil with nan values until complete same lenth.
    ## Nan values are not taking into account when perform calculations
    
    b_ = len(B_distributions[0])
    for i, v, b in zip(I_distributions,
                       Variance_distributions,
                       B_distributions):
        if len(b) < b_:
            i+=[np.nan]*(abs(len(b) - b_))
            v+=[np.nan]*(abs(len(b) - b_))
            b+=[np.nan]*(abs(len(b) - b_))
    
    I_distributions = np.array(I_distributions)
    Variance_distributions = np.array(Variance_distributions)
    Epsilon_distributions = (np.array(B_distributions)-1)/Vol_PSF
    
    
    return I_distributions, Variance_distributions, Epsilon_distributions


def Get_distributions_2 (data_, N_):
    '''
    
    Parameters
    ----------
    data_ : ndarray
        Array with list of intensity values for different diffusion coefficients.
    N_ : 1D int array
        Each value is the number of data points use to calculate each entries of a distribution.

    Returns
    -------
    Epsilon_distributions : ndarray
        Molecular brightness distributions for different diffusion coefficients.


    '''

    B_values_for_each_D_ = []
    
    for d in data_:
            
        B_values_ = []

        for p in N_:
            B_ = B_from_n_points(d, int(p))[2]
        
            B_values_.append(list(B_))

        B_values_for_each_D_.append(B_values_)


    ## To convert list into arrays fullfil with nan values until complete same lenth.
    ## Nan values are not taking into account when perform calculations
    b_ = len(B_values_for_each_D_[0][0])

    for B_ in B_values_for_each_D_:

        for b in B_:
            if len(b) < b_:
                b+=[np.nan]*(abs(len(b) - b_))
    
    Epsilon_values_for_each_D_ = (np.array(B_values_for_each_D_)-1)/0.3536
    
    return np.array(Epsilon_values_for_each_D_)



#%%

def B_apparent_brightness_consecutives_lines(Lines, tp='', window_points = 100, delta_line=100, window_shift=1):
    '''
    Parameters
    ----------
    Lines : ndarray
        2D array. Each row is a Line
        
    tp : float
        pixel dwell time. Value must be in seconds

    window_points : int
        Number of lines for N&B analysis.
        By default 100 lines are consider for N&B aalysis.
    
    delta_line : TYPE, optional
        Distance between lines in a package. The default is 100.

    window_shift : TYPE, optional
        Displacement of the Package of line for epsilon calculation. The default is 1.
    
    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    B Apparent brightness kimogram.

    '''

    if window_shift=='':
        raise ValueError ('The parameter window shift can not be empty')

    if tp=='':
        raise ValueError ('The parameter pixel dwell time (tp) can not be empty')

        
    B=[]

    p=0


    l_f = len(Lines)-window_points ## l_f = final line parameter use for index. 

    for i in tqdm.trange(0, l_f, window_shift):

        pack_of_lines = Lines[0+i:window_points+i]

        k_medio = (np.mean(pack_of_lines, axis=0)).ravel()
        
        # Variance is calculated as intensity distribution 2nd moment --->  Var = <(I^2)> - (<I>)^2
        VAR = (np.mean((pack_of_lines)**2,axis=0)-((np.mean(pack_of_lines, axis=0))**2)).ravel() 

        B_line = np.array([1.0]*len(k_medio)).ravel()
        
        for k in range(0,len(k_medio)):
            if k_medio[k] == 0 or VAR[k]==k_medio[k]:
                B_line[k]= 1  
            else:
                B_line[k] = VAR[k]/k_medio[k]
        B.append(list(B_line))
   

        
        p+=1

    
    B = np.asarray(B)

    return B


def MATE_BA(Lines, tp='', window_points = 100, delta_line=100, window_shift=1, correlated_lines=True):

    '''
    Parameters
    ----------
    Lines : ndarray
        2D array. Each row is a Line
        
    tp : float
        pixel dwell time. Value must be in seconds.
        
    window_points : int
        Number of lines for epsilon calculation.
        By default 100 lines are consider.
    
    delta_line : int, optional
        Distance between lines in a package. The default is 100.

    window_shift : TYPE, optional
        Displacement of the Package of line for epsilon calculation. The default is 1.
    
    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    B Apparent brightness kimogram.
    
    '''
    
    if window_shift=='':
        raise ValueError ('The parameter window shift can not be empty. \n ¿How many lines must be jumped between intensity values?')

    if tp=='':
        raise ValueError ('The parameter pixel dwell time (tp) can not be empty')

    
    B=[]

    p=0

    l_f = window_points*delta_line  ## l_f = final line parameter use for index. 
    #==============================================================================
    #             Initilization of line's index values for B calculation
    #==============================================================================  
    index=(np.arange(0, l_f, delta_line))+window_shift*p

    # print(index)
    for p in tqdm.trange(0,delta_line,window_shift):

        pack_of_lines=[]    

        #==========================================================================
        #             Update of line's index values for B calculation
        #==========================================================================
        index=(np.arange(0, l_f, delta_line))+window_shift*p
            
        for i in index:
            pack_of_lines.append(list(Lines[i]))

        pack_of_lines=np.asarray(pack_of_lines)        
                
        k_medio = (np.mean(pack_of_lines, axis=0)).ravel()
        
        # Variance is calculated as intensity distribution 2nd moment --->  Var = <(I^2)> - (<I>)^2
        VAR = (np.mean((pack_of_lines)**2,axis=0)-((np.mean(pack_of_lines, axis=0))**2)).ravel() 

        B_line = np.array([1.0]*len(k_medio)).ravel()
        
        for k in range(0,len(k_medio)):
            if k_medio[k] == 0 or VAR[k]==k_medio[k]:
                B_line[k]= 1  

            else:
                B_line[k] = VAR[k]/k_medio[k]

        B.append(list(B_line))

        
    B = np.asarray(B)
    
    return B