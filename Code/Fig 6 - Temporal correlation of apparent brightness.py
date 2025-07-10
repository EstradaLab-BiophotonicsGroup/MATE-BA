# -*- coding: utf-8 -*-
"""
Created on Wed May  7 12:21:37 2025

@author: Ignacio Sallaberry
"""

import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path)
    
from Functions import *

#%%
intensity_values = pd.read_csv(path[:-5]+r'\Simulation Data\Fig 6\Simulation 6.csv.gz', compression='gzip', dtype=np.int16)['0'].to_numpy()
kimogram_1_mers = Kimogram(intensity_values,64)

t_p=5e-6

#%%
## Molecular brightness calculated from N = 100 consecutive correlated lines                                                                                            
# ---------------------------------------------------------------------------------

window_points = 100 
window_shift = 1 
delta_line = 1    

B_carpet_1_mers_correlated_lines, B_t_1_mers_correlated_lines = B_apparent_brightness_consecutives_lines(kimogram_1_mers, t_p,
                                                                                                         window_points=window_points,
                                                                                                         delta_line=delta_line,
                                                                                                         window_shift=window_shift)

#%%
## Molecular brightness calculated from N = 100 uncorrelated (interspersed) lines                                                                                            
# ---------------------------------------------------------------------------------

window_points = 100
window_shift = 1 
delta_line = int(len(kimogram_1_mers)/window_points)    

B_carpet_1_mers_interspersed_lines, B_t_1_mers_interspersed_lines = MATE_BA(kimogram_1_mers, t_p,
                                                                            window_points=window_points,
                                                                            window_shift=window_shift,
                                                                            delta_line=delta_line)

#%%
# # ---------------------------------------------------------------------------------

# # ██████       ██████      ██████      ███    ███     ██████  
# # ██   ██     ██          ██    ██     ████  ████     ██   ██ 
# # ██████      ██          ██    ██     ██ ████ ██     ██████  
# # ██          ██          ██    ██     ██  ██  ██     ██   ██ 
# # ██           ██████      ██████      ██      ██     ██████  
                                                                                            
# # ---------------------------------------------------------------------------------
epsilons_1_mers_correlated_lines = (B_carpet_1_mers_correlated_lines[75000:100000]-1)/0.3536/t_p
epsilons_1_mers_interspersed_lines = (B_carpet_1_mers_interspersed_lines-1)/0.3536/t_p
dr=0 #Correlation distance

print(len(epsilons_1_mers_correlated_lines[0:epsilons_1_mers_interspersed_lines.shape[0]]))
sigma=[5, 15]
#=============================================================
#                   G(B,T)
#=============================================================
G_eps_1_mers_correlated_lines, T_1_mers_correlated_lines = pCOMB(epsilons_1_mers_correlated_lines[0:epsilons_1_mers_interspersed_lines.shape[0]],
                                                               Time=B_t_1_mers_correlated_lines[0:epsilons_1_mers_interspersed_lines.shape[0]], 
                                                               window_shift=window_shift, dr=dr)

G_eps_1_mers_interspersed_lines, T_1_mers_interspersed_lines = pCOMB(epsilons_1_mers_interspersed_lines, 
                                                                   Time=B_t_1_mers_interspersed_lines, 
                                                                   window_shift=window_shift, dr=dr)


#%%
###### ----------------------------------------------------------
######  Plot Figure 6 - EPSILONS
###### ----------------------------------------------------------

fig = plt.figure(figsize=(6, 6))
   
colors = ['#009FFD', '#F40076'] 
# Figure consists on 3 subfigures (first row: 2 columns, second row: 1 column full width)

# First we create two subfigures (sfigs) in Figure (fig). Two rows and on column per row
sfigs = fig.subfigures(2, 1, height_ratios=[1, 1])  

# In top subfigure (sfigs[0]) we create two subfigures. One row and two columns
top_row = sfigs[0].subfigures(1, 2)

# Plot in left top sub-subfigure
X_TICKS = list(np.arange(0, 2.2, 1))
plot_pCOMB_carpet_and_pCOMB_perfil(top_row[0], T_1_mers_correlated_lines, G_eps_1_mers_correlated_lines,
                                     sigma=[5, 15], pCOMB_curve_color=colors[0], dr=0,
                                     colorbar_x_position_in_figure=-0.7,
                                      letter='A',
                                     y_min_vertical_line=0.45, x_ticks=X_TICKS, X_lim=[-0.2, 2.2])

# Plot in right top sub-subfigure
X_TICKS = list(np.arange(0, 0.022, 0.01))
plot_pCOMB_carpet_and_pCOMB_perfil(top_row[1], T_1_mers_interspersed_lines, G_eps_1_mers_interspersed_lines,
                                      sigma=[5, 15], pCOMB_curve_color=colors[1], dr=0,
                                      colorbar_x_position_in_figure=-0.7,
                                      letter='B',
                                      y_min_vertical_line=0.46, x_ticks=X_TICKS, X_lim=[-0.0025, 0.022])


# Plot full-width bottom subfigure (sfigs[1])
plot_hist_fig(sfigs[1],
              epsilons_1_mers_correlated_lines,
              epsilons_1_mers_interspersed_lines,
              bins=np.arange(-1.5,6,0.05)/t_p,
              color_B_1=colors[0], color_B_2=colors[1],
              letter='C', B_simulated=3e5,
              X_label=r'$\epsilon$ (cpms)',
              legend_sim_label=r'$\epsilon_{sim}$',
              legend_median_label=r'$\epsilon_{median}$')

# Customize bottom subfigure
sfigs[1].subplots_adjust(top=1,
                  bottom=0.15,
                  left=0.095,
                  right=0.95,
                  hspace=0,
                  wspace=0)

plt.show()
