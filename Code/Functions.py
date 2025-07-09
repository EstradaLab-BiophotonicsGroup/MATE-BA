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
from matplotlib import ticker, axis
import matplotlib.patches as mpatches

from scipy.stats import mannwhitneyu 
from scipy.stats import norm
from scipy.stats import median_abs_deviation

from scipy.ndimage import gaussian_filter

import addcopyfighandler  ##allows to Ctrl+C matplotlib figures

import csv

from lfdfiles import SimfcsB64 as lfd

# plt.rcParams["mathtext.fontset"] = "dejavusans" ## Default appearence of matplotlib typography
plt.rcParams["mathtext.fontset"] = "cm" ## Latex tipographic style
plt.rcParams['font.family'] = 'Arial'


fontsize=12
point_size = 4
lw=1
markeredgewidth=0.5  #width of symbol edge's
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
        VAR.append(np.var(data[i:i+points], ddof=1))

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
def Line_time (total_num_lines, total_num_columns, tp, tr=0):
    '''
    
    Parameters
    ----------
    total_num_lines : int

    total_num_columns : int

    tp : float
        pixel dwell time.
        
    tr : float, optional
        Return time in a real microscope. 
        If you are working with simulations tr is 0 by default. In a simulation the line time is tp*total_num_columns
        If you are working with microscope adquired data:
            -) tr would be 0 if tp already has incorporated the microscope pixel return time.
            -) tr would not be 0 if is just the pixel dwell time does not consider the microscope pixel return time.
    
    Returns
    -------
    1d-array: .

    '''
    sampling_frequency = 1/(tp*total_num_columns+tr)
    Time=[]
    
    for i in range (1,total_num_lines+1):
                     Time.append(i*1/sampling_frequency)
                     
    return Time



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

    Time = Line_time(Lines[:,0].size, Lines[0,:].size, tp)
    
    time=[]   
    B=[]

    p=0

    l_f = len(Lines)-window_points ## l_f = final line parameter use for index. 

    for i in tqdm.trange(0, l_f, window_shift):

        pack_of_lines = Lines[0+i:window_points+i]

        k_medio = (np.mean(pack_of_lines, axis=0)).ravel()
        
        # Variance is calculated as intensity distribution 2nd moment --->  Var = <(I^2)> - (<I>)^2 but with the denominator as N-1 (ie: ddof=1)
        VAR = np.var(pack_of_lines, axis=0, ddof=1).ravel()

        B_line = np.array([1.0]*len(k_medio)).ravel()
        
        for k in range(0,len(k_medio)):
            if k_medio[k] == 0 or VAR[k]==k_medio[k]:
                B_line[k]= 1  
                
            else:
                B_line[k] = VAR[k]/k_medio[k]
                
        B.append(list(B_line))
   
        time.append(Time[0]*(p+1))
        
        p+=1

    
    B = np.asarray(B)

    return B, time    

def MATE_BA(Lines, tp='', tr=0, window_points = 100, delta_line=100, window_shift=1, correlated_lines=True):

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

    Time = Line_time(Lines[:,0].size, Lines[0,:].size, tp, tr=tr)
    
    
    time = []

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
        
        # Variance is calculated as intensity distribution 2nd moment --->  Var = <(I^2)> - (<I>)^2 but with the denominator as N-1 (ie: ddof=1)
        VAR = np.var(pack_of_lines, axis=0, ddof=1).ravel() 

        B_line = np.array([1.0]*len(k_medio)).ravel()
        
        for k in range(0,len(k_medio)):
            if k_medio[k] == 0 or VAR[k]==k_medio[k]:
                B_line[k]= 1  

            else:
                B_line[k] = VAR[k]/k_medio[k]

        B.append(list(B_line))

        time.append(Time[0]*(p+1))


    B = np.asarray(B)
    
    return B, time

#%%

def pCOMB_analysis(carpet, C0, C1, Time, window_shift, reverse_PCOMB, return_time=0):

    '''
    Parameters
    ----------
    carpet : ndarray tipically (100k rows, 256 columns)
        carpet = Kimogram or B Carpet    
    
    C0 : int
        first column for correlate.
    C1 : TYPE
        Second column for correlate.

    tp: float
        pixel dwell time.
        
    return_time : float, optional
        line time return. The default is 0.
        
    Returns
    -------
    G : 1D-array
        Correlation between columns C0 and C1
        
    Tau : 1D-array
        Correlation time
    
    '''
    # First column
    C0=carpet [ : ,C0]
    # Second column
    C1=carpet [ : ,C1]    
    
    if reverse_PCOMB:
        C0, C1 = C1, C0

    ######################################
    
    """Return linear correlation of two arrays using DFT."""
    size = C0.size
    
    # subtract mean and pad with zeros to twice the size
    C0_mean = C0.mean()
    C1_mean = C1.mean()

    C0 = np.pad(C0-C0_mean, C0.size//2, mode='constant')
    C1 = np.pad(C1-C1_mean, C1.size//2, mode='constant')
    

    # forward DFT
    C0 = np.fft.rfft(C0)
    C1 = np.fft.rfft(C1)
    # multiply by complex conjugate
    G = C0.conj() * C1
    # reverse DFT
    G = np.fft.irfft(G)
    # positive delays only
    G = G[:size//2]
        
    # normalize with the averages of a and b
    
    G = G / (size * C0_mean * C1_mean)

    ######################################
    
    # #Time delay
    Tau=np.arange(1,len(G)+1)*Time[0]
    ######################################

    return np.array(G), np.array(Tau)

def pCOMB(B_carpet , Time, window_shift=1, dr=0, reverse_PCOMB=False, return_time=0):
    '''
    Parameters
    ----------
    B_carpet  : ndarray tipically (100k rows, 256 columns)
        carpet = Kimogram or B Carpet  
        
    dr : int
        pCF distance.
    
    tp : TYPE
        DESCRIPTION.

    return_time : float, optional
        line time return. The default is 0.
        

    Returns
    -------
    G : ndarray 
        Matrix of pCF analisys between at distance.
        -) Columns are the pixel position
        -) Rows are the correlation analysis
        
    Tau : ndarray
          Correlation time
          -) Columns are the pixel position
          -) Rows are the delay time
    '''

    Size = len(B_carpet [0]) #number of pixels in a line

    G=[]
    T=[]
    for i in tqdm.trange(Size-dr):
        pCOMB_results = pCOMB_analysis(B_carpet  ,i ,i+dr, Time, window_shift,
                                       reverse_PCOMB, return_time)

        G.append(pCOMB_results[0])
        T.append(pCOMB_results[1])

    return np.array(G).transpose(), np.array(T).transpose()

#%%
def plot_pCOMB_carpet_and_pCOMB_perfil(f, 
               T, G, dr ='',  xlabel='', ylabel='', label_colorbar='',
                                     vmin=None, vmax=None, sigma=[0,0],
                                     colorbar_x_position_in_figure=0,
                                     colorbar_height_in_figure = 1,
                                     letter ='',
                                     pCOMB_curve_color='C0', x_ticks=None, y_min_vertical_line=0, D_simulated=1, X_lim=None):
    """
    

    Parameters
    ----------
    f : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    G : TYPE
        DESCRIPTION.
    dr : TYPE, optional
        DESCRIPTION. The default is ''.
    xlabel : TYPE, optional
        DESCRIPTION. The default is ''.
    ylabel : TYPE, optional
        DESCRIPTION. The default is ''.
    label_colorbar : TYPE, optional
        DESCRIPTION. The default is ''.
    vmin : TYPE, optional
        DESCRIPTION. The default is None.
    vmax : TYPE, optional
        DESCRIPTION. The default is None.
    sigma : TYPE, optional
        DESCRIPTION. The default is [0,0].
    colorbar_x_position_in_figure : TYPE, optional
        DESCRIPTION. The default is 0.
    colorbar_height_in_figure : TYPE, optional
        DESCRIPTION. The default is 1.
    letter : TYPE, optional
        DESCRIPTION. The default is ''.
    pCOMB_curve_color : TYPE, optional
        DESCRIPTION. The default is 'C0'.
    x_ticks : TYPE, optional
        DESCRIPTION. The default is None.
    y_min_vertical_line : TYPE, optional
        DESCRIPTION. The default is 0.
    D_simulated : TYPE, optional
        DESCRIPTION. The default is 1.
    X_lim : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """

    
    
    ###### ----------------------------------------------------------
    ######  Prepare GB (Tau, dr) interpolation for plotting
    ###### ----------------------------------------------------------

    x1 = np.geomspace(1, len(G), 256, dtype=int, endpoint = False)  #
    t_lineal = T[:,0]
    t_log = np.geomspace(t_lineal[0], t_lineal[-1], 256)
    
    G_basura = []
    for i in x1:
        G_basura.append(list(G[i]))
    G = np.asarray(G_basura).transpose()

    
    t = []
    for i in x1:
        t.append(t_lineal[i])
    t_lineal = np.asarray(t)

    G_log = np.empty_like(G)
    for i, gi in enumerate(G):
        G_log[i] = np.interp(t_log, t_lineal, gi)
    G_log = gaussian_filter(G_log, sigma = sigma)


    ###### ----------------------------------------------------------
    ######  Plot interpoled GB(Tau, dr) in left axes
    ###### ----------------------------------------------------------

    y = np.arange(G_log.shape[0])

    ax = f.subplots(1, 2, width_ratios=[3, 1.5], sharey=True)

    ## Prepare colormap for plot
    cmap = plt.colormaps['jet'].copy()
    cmap.set_bad(color='k')

    
    if vmin is None:
        vmin=np.min(G_log)

    if vmax is None:
        vmax=np.max(G_log)
        
    if  vmin is not None:
        G_log_masked = np.asarray(G_log)
        G_log_masked = np.ma.masked_where(G_log_masked <= vmin, G_log_masked)
        ax[0].set_facecolor('k')
    
        ## pcolor = pcolor(X, Y, Z)        
        im = ax[0].pcolor(y.transpose(), t_log, (G_log_masked).transpose(), shading="nearest",
                       cmap=plt.cm.jet, vmin=vmin, vmax=vmax, rasterized=True)
      
    else:
        ## pcolor = pcolor(X, Y, Z)
        im = ax[0].pcolor(y.transpose(), t_log, (G_log).transpose(), shading="nearest",
                       cmap=plt.cm.jet, vmin=vmin, vmax=vmax, rasterized=True)

    
    ax[0].set_ylabel(r'$\tau$ (s)', fontsize=fontsize, labelpad=0)
    
    # Customize lines of left axes
    ax[0].set_xticks(ticks=np.arange(0, len(G_log[:,0]), len(G_log[:,0])/4))
    ax[0].tick_params(axis='both', labelsize=fontsize-4) ## this change the size of the number in the carpet
    ax[0].tick_params(which='minor', length=1.25, width=1)
    ax[0].tick_params(which='major', length=2.25, width=1)
    ax[0].xaxis.tick_top() ## This makes the numbers of x axis appear in the top of the plot
    ax[0].set_yscale("log")
    ax[0].invert_yaxis()


    ###### ----------------------------------------------------------
    ######  Colorbar
    ###### ----------------------------------------------------------
    ## Colorbar position and separation from axis
    colorbar_position_in_figure = ax[0].inset_axes([colorbar_x_position_in_figure, 0, 0.05, colorbar_height_in_figure])               
    colobar_set_offset_position = 'right'
    colorbar_pad = 0
    ## Colorbar numbers format
    cbarformat = ticker.ScalarFormatter()
    cbarformat.set_scientific('%.2e')
    cbarformat.set_powerlimits((0,0))
    cbarformat.set_useMathText(True)
    
    ## Plot colorbar
    cbar = f.colorbar(im, ax=ax[0], orientation='vertical', location='left', cax=colorbar_position_in_figure,
                      format=cbarformat, pad=0)
    ## Customize colorbar scitific notation and alignment
    cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize-5)   ## This change the size of the "x10-1" number outside the colorbar
    cbar.ax.yaxis.set_offset_position(colobar_set_offset_position) ## This change the position of the "x10-1" number outside the colorbar.
                                              ## It aligns the left size of the "x10-1" whit the left size of the colorbar. The alternative is '"right".
      
    ## Customize colorbar scitific notation and alignment
    cbar.ax.tick_params(labelsize=fontsize-4)   ## This change the size of the numbers in the colobar
    cbar.set_label(label=label_colorbar, size=2.5)   ## This set the label in the colorbar
    



    ###### ----------------------------------------------------------
    ######  Plot pCOMB columns perfil in right axes
    ###### ----------------------------------------------------------
    ax[1].semilogy(np.mean(G_log, axis=0), t_log, '-', color=pCOMB_curve_color, label=r'$G_{B}(\tau , \delta r = %s )$' % (dr), ms=1)
   
    if x_ticks == None:
        ax[1].set_xticks(ticks=np.arange(0,max(np.mean(G_log, axis=0))*1.2,0.004))
    
    else:
        ax[1].set_xticks(ticks=x_ticks)
        
        
    ax[1].ticklabel_format(style='scientific', axis='x', scilimits=(0,0),  useMathText=True)  # Apply scientific notation to y-axis
    ax[1].xaxis.get_offset_text().set_fontsize(fontsize-5)
    ax[1].yaxis.set_offset_position('right')

    ax[1].tick_params(axis='both', labelsize=fontsize-3) ## this change the size of the number in the perfil plot
    ax[1].tick_params(which='minor', length=1.25, width=1)
    ax[1].tick_params(which='major', length=2.25, width=1)

    ###### ----------------------------------------------------------
    ######  Estimate correlation characteristic time using interpolation between the upper top point and the bottom point from the np.e decay of GB
    ###### ----------------------------------------------------------
    # Calculate GB perfil
    G_B_log_perfil = np.mean(G_log, axis=0)

    # Calculate the nearest point of the 1/e decay value of GB
    nearest_G_value = G_B_log_perfil[min(range(len(G_B_log_perfil)), key = lambda i: abs(G_B_log_perfil.ravel()[i]-G_B_log_perfil[0]/np.e))]

    # Get index of the nearest point of the 1/e decay value of GB
    nearest_G_value = list(G_B_log_perfil).index(nearest_G_value)
    
    # Get the two points of GB nearest to the 1/e decay value of GB
    if G_B_log_perfil[nearest_G_value]>G_B_log_perfil[0]/np.e:
        G_B_log_perfil_first_value_to_interopole = nearest_G_value
        G_B_log_perfil_second_value_to_interopole = nearest_G_value+1
    else:        
        G_B_log_perfil_first_value_to_interopole = nearest_G_value-1
        G_B_log_perfil_second_value_to_interopole = nearest_G_value

    # Prepare x points for interpolation    
    x_vals_for_interpolation = np.arange(t[G_B_log_perfil_first_value_to_interopole], t[G_B_log_perfil_second_value_to_interopole],
                                                  (t[G_B_log_perfil_second_value_to_interopole]-t[G_B_log_perfil_first_value_to_interopole])/50)

    # Get GB's interpolated points between the first and second points of GB to interpolate
    interpolated_G_B_values = np.interp(x_vals_for_interpolation,
                                        [t[G_B_log_perfil_first_value_to_interopole], t[G_B_log_perfil_second_value_to_interopole]],
                                        [G_B_log_perfil[G_B_log_perfil_first_value_to_interopole], G_B_log_perfil[G_B_log_perfil_second_value_to_interopole]])

    # Search for the nearest GB interpolated points to the 1/e decay value of GB                                 
    nearest_G_value = interpolated_G_B_values[min(range(len(interpolated_G_B_values)), key = lambda i: abs(interpolated_G_B_values[i]-G_B_log_perfil[0]/np.e))]
    nearest_G_value = list(interpolated_G_B_values).index(nearest_G_value)

    print('pCOMB(0) characteristic time = %s' % x_vals_for_interpolation[nearest_G_value])
    

    if X_lim!=None:
        ax[1].set_xlim(X_lim[0], X_lim[1])    

    ###### ----------------------------------------------------------
    ######  Plot lines of the theorical 1/e decay of GB value and the estimated characteristic time by interpolation of GB points
    ###### ----------------------------------------------------------
    X_lenght_horizonatal_line = interpolated_G_B_values[nearest_G_value]/ax[1].get_xlim()[1]
    ax[1].axhline(y=x_vals_for_interpolation[nearest_G_value], xmin=0, xmax=X_lenght_horizonatal_line, ls='--', color=pCOMB_curve_color, lw=lw)
    ax[1].axvline(interpolated_G_B_values[nearest_G_value],ymin=y_min_vertical_line, ymax=1, ls='--', color=pCOMB_curve_color, lw=lw)

    ax[1].axhline(y=((0.3)**2)/6*D_simulated, xmax=1, ls='--', color='k', lw=lw)

    # Customize lines of right axes
    ax[1].yaxis.set_inverted(True)  # inverted axis with autoscaling
    ax[1].xaxis.set_ticks_position('top')
    ax[1].xaxis.set_label_position('top') 
    ax[1].text(1.45, 1.15, letter, transform=ax[1].transAxes,
        fontsize=fontsize-1, va='top', ha='right')

    # Customize the figure
    f.subplots_adjust(top=0.85,
                      bottom=0.05,
                      left=0.40,
                      right=0.9,
                      hspace=0.9,
                      wspace=0.075)
        


def plot_hist_fig(f, B_1, B_2, bins, color_B_1 = '#ff000d', color_B_2 = '#ff000d', letter='', B_simulated = 1.53,
                      X_label = r'B(cp $t_{p}$)',
                      legend_sim_label=r'B$_{sim}$', legend_median_label=r'B$_{median}$'):
    """
    

    Parameters
    ----------
    f : TYPE
        DESCRIPTION.
    B_1 : TYPE
        DESCRIPTION.
    B_2 : TYPE
        DESCRIPTION.
    bins : TYPE
        DESCRIPTION.
    color_B_1 : TYPE, optional
        DESCRIPTION. The default is '#ff000d'.
    color_B_2 : TYPE, optional
        DESCRIPTION. The default is '#ff000d'.
    letter : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    
    
    ax = f.subplots()
    
    bins=bins
    
    ###### ----------------------------------------------------------
    ######  Plot brightness histogram in axes
    ###### ----------------------------------------------------------
    ax.hist(B_1.ravel(), bins=bins, color=color_B_1, log=True, hatch='//', alpha=1, edgecolor=color_B_1, histtype='step', label='Consecutive lines')
    ax.hist(B_2.ravel(), color=color_B_2 ,bins=bins, log=True, hatch=None, alpha=0.65, edgecolor=None, histtype='bar', label='Interpersed lines')
    
    ax.axvline(np.median(B_1.ravel()), ls='-', color=color_B_1, lw=lw, alpha=1)
    ax.axvline(np.median(B_2.ravel()), ls='-', color=color_B_2, lw=lw, alpha=1)
    ax.axvline(B_simulated, ls='--', color='k', lw=lw, alpha=0.75)
    
    ax.set_ylim(0.9,5e8)

    ax.text(1.045, 1, letter, transform=ax.transAxes,
        fontsize=fontsize-1, va='top', ha='right')

    # Customize lines in axes
    ax.tick_params(which='minor', length=1.25, width=1)
    ax.tick_params(which='major', length=2.25, width=1)
    ax.tick_params(axis='both', labelsize=fontsize-3)

    ## affect borders lines plot making theme thinners
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.6)

    ## Figure legend
    handles, labels = plt.gca().get_legend_handles_labels()
    correlated_legend_box = mpatches.Patch(lw=lw*0.75, hatch='///', facecolor='white', edgecolor='k', alpha=0.75)
    uncorrelated_legend_box = mpatches.Patch(lw=lw*0.75, color='k', alpha=0.75)
    
    empty_line = Line2D([0], [0], ls='', color='k', lw=lw*0.75)
    dash_line = Line2D([0], [0], ls='--', color='k', lw=lw*0.75)
    solid_line = Line2D([0], [0], ls='-', color='k', lw=lw)
    
    handles = [correlated_legend_box, uncorrelated_legend_box, dash_line, solid_line]
    labels = [r'Correlated lines', r'Uncorrelated (Interspersed) lines', legend_sim_label, legend_median_label]
    f.legend(handles, labels, loc='upper center', fontsize=fontsize, ncol=2)
    
    plt.xlabel(X_label, fontsize = fontsize)
    plt.ylabel(' Frequency', loc='center', fontsize = fontsize, labelpad = labelpad*4.5)
    

    