# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:23:06 2024

@author: Ignacio Sallaberry
"""

# ruff: noqa

import os

path = os.path.dirname(__file__)
import sys

sys.path.append(path)

from Functions import *

# %%
SEED = 0  ## Use SEED = 0 to reproduce Figure 1 and 2 or load saved data
# SEED = None  ## Use SEED = None to obtain new random values


## Amount of data points used to calculate each entry
N = np.logspace(2, 8, num=100)

N_1e2_index = min(range(len(N)), key=lambda i: abs(N[i] - 1e2))
N_1e4_index = min(range(len(N)), key=lambda i: abs(N[i] - 1e4))
N_1e6_index = min(range(len(N)), key=lambda i: abs(N[i] - 1e6))

tp = 1e-6

Epsilon_simulated = 500000
# %%
intensity_points = pd.read_csv(
    path[:-5] + r"\Simulation Data\Fig 1 - 2 and 4\Simulation Fig 1 - 2 and 4.csv.gz",
    compression="gzip",
    dtype=np.int16,
)["0"].to_numpy()

## Generate randomly mixed intensity values
points_random = Random_data(intensity_points, seed=SEED)


# %%
## Get distributions of Mean Intensity, Variance and Epsilon in cpms

(
    Mean_I_distributions_correlated_values,
    Variance_distributions_correlated_values,
    Epsilon_distributions_correlated_values,
) = [d / tp for d in Get_distributions(intensity_points, N)]

(
    Mean_I_distributions_UN_correlated_values,
    Variance_distributions_UN_correlated_values,
    Epsilon_distributions_UN_correlated_values,
) = [d / tp for d in Get_distributions(points_random, N)]

# %%
## Load save data to reproduce Figure 1

# Mean_I_distributions_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - I mean distributions - correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Variance_distributions_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Variance distributions - correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Epsilon_distributions_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Epsilon distributions - correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()

# Mean_I_distributions_UN_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - I mean distributions - UN correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Variance_distributions_UN_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Variance distributions - UN correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()
# Epsilon_distributions_UN_correlated_values = pd.read_csv(path[:-5]+r'\Data for Figures\Fig 1 - Epsilon distributions - UN correlated data.csv.gz', compression='gzip', dtype=np.float64).to_numpy()

# %%
## Comparison between CORR (hatched) and UNCORR (solid) distributions performin a Mann Whitney test

Mann_Whitney_U_test_Fig_1A = mannwhitneyu(
    Epsilon_distributions_correlated_values[N_1e2_index],
    Epsilon_distributions_UN_correlated_values[N_1e2_index],
    nan_policy="omit",
)[1]

Mann_Whitney_U_test_Fig_1B = mannwhitneyu(
    Epsilon_distributions_correlated_values[N_1e4_index],
    Epsilon_distributions_UN_correlated_values[N_1e4_index],
    nan_policy="omit",
)[1]

Mann_Whitney_U_test_Fig_1C = mannwhitneyu(
    Epsilon_distributions_correlated_values[N_1e6_index],
    Epsilon_distributions_UN_correlated_values[N_1e6_index],
    nan_policy="omit",
)[1]

print(
    "Fig 1A: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s "
    % Mann_Whitney_U_test_Fig_1A
)
print(
    "Fig 1B: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s "
    % Mann_Whitney_U_test_Fig_1B
)
print(
    "Fig 1C: p value of Mann-Whitney test between correlated and uncorrelated molecular brightness distributions is: %s "
    % Mann_Whitney_U_test_Fig_1C
)


# %%

MEDIAN_I_mean_correlated = np.nanmedian(Mean_I_distributions_correlated_values, axis=1)
MAD_I_mean_correlated = median_abs_deviation(
    Mean_I_distributions_correlated_values, axis=1, nan_policy="omit"
)
MEDIAN_I_mean_UN_correlated = np.nanmedian(
    Mean_I_distributions_UN_correlated_values, axis=1
)
MAD_I_mean_UN_correlated = median_abs_deviation(
    Mean_I_distributions_UN_correlated_values, axis=1, nan_policy="omit"
)


MEDIAN_VAR_correlated = np.nanmedian(Variance_distributions_correlated_values, axis=1)
MAD_VAR_correlated = median_abs_deviation(
    Variance_distributions_correlated_values, axis=1, nan_policy="omit"
)
MEDIAN_VAR_UN_correlated = np.nanmedian(
    Variance_distributions_UN_correlated_values, axis=1
)
MAD_VAR_UN_correlated = median_abs_deviation(
    Variance_distributions_UN_correlated_values, axis=1, nan_policy="omit"
)


MEDIAN_epsilon_correlated = np.nanmedian(
    Epsilon_distributions_correlated_values, axis=1
)
MAD_epsilon_correlated = median_abs_deviation(
    Epsilon_distributions_correlated_values, axis=1, nan_policy="omit"
)
MEDIAN_epsilon_UN_correlated = np.nanmedian(
    Epsilon_distributions_UN_correlated_values, axis=1
)
MAD_epsilon_UN_correlated = median_abs_deviation(
    Epsilon_distributions_UN_correlated_values, axis=1, nan_policy="omit"
)


# %%
## Plot Figure 1
## Molecular brightness distribution where each entry is calculated using N=100, 100k and 1M data points
fig, axs = plt.subplots(3, 1, figsize=(3, 3), sharex=True, sharey=False, dpi=300)

letters = ["A", "B", "C"]
axs[0].text(2.55 / tp, 4400, letters[0], fontsize=11)
axs[1].text(2.55 / tp, 440, letters[1], fontsize=11)
axs[2].text(2.55 / tp, 13.5, letters[2], fontsize=11)

axs[0].text(
    0.85 / tp,
    3500,
    "$\#_{entries}= %s$"
    % str(
        len(Epsilon_distributions_UN_correlated_values[N_1e2_index])
        - np.count_nonzero(
            np.isnan(Epsilon_distributions_UN_correlated_values[N_1e2_index])
        )
    ),
    fontsize=9,
)
axs[1].text(
    0.85 / tp,
    350,
    "$\#_{entries}= %s$"
    % str(
        len(Epsilon_distributions_UN_correlated_values[N_1e4_index])
        - np.count_nonzero(
            np.isnan(Epsilon_distributions_UN_correlated_values[N_1e4_index])
        )
    ),
    fontsize=9,
)
axs[2].text(
    0.85 / tp,
    11.5,
    "$\#_{entries}= %s$"
    % str(
        len(Epsilon_distributions_UN_correlated_values[N_1e6_index])
        - np.count_nonzero(
            np.isnan(Epsilon_distributions_UN_correlated_values[N_1e6_index])
        )
    ),
    fontsize=9,
)


bins = np.arange(-2.5 / tp, 4 / tp, 0.03 / tp)
log = False
density = False

for i, ax in zip([N_1e2_index, N_1e4_index, N_1e6_index], axs.flatten()):
    ax.hist(
        Epsilon_distributions_UN_correlated_values[i],
        bins=bins,
        log=log,
        density=density,
        color="olivedrab",
        alpha=0.5,
        edgecolor="None",
    )
    ax.hist(
        Epsilon_distributions_correlated_values[i],
        bins=bins,
        log=log,
        density=density,
        color="white",
        hatch="//////",
        alpha=0.5,
        edgecolor="lightgreen",
    )
    ax.axvline(
        np.nanmedian(Epsilon_distributions_UN_correlated_values[i]),
        ls="-",
        color="olivedrab",
        lw=0.75 * lw,
        alpha=0.75,
    )
    ax.axvline(
        np.nanmedian(Epsilon_distributions_correlated_values[i]),
        ls="-",
        color="lightgreen",
        lw=lw,
    )
    ax.axvline(Epsilon_simulated, ls="--", color="k", lw=0.75 * lw, alpha=0.75)

for ax in axs.flatten():
    ax.tick_params(which="minor", length=1.25, width=1)
    ax.tick_params(which="major", length=2.25, width=1)

    ax.tick_params(axis="both", labelsize=fontsize - 3)
    ax.tick_params(axis="both", labelsize=fontsize - 3)

    ## affect borders lines plot making theme thinners
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(0.6)

plt.xlabel("$\epsilon$ (cpms)", fontsize=fontsize)
plt.ylabel(
    "                                Frequency",
    loc="center",
    fontsize=fontsize,
    labelpad=labelpad * 5,
)

axs[2].ticklabel_format(
    style="scientific", axis="x", scilimits=(6, 6), useMathText=True
)  # Apply scientific notation to y-axis
axs[2].xaxis.get_offset_text().set_fontsize(fontsize - 5)

axs[0].set_xlim(-1.51 / tp, 2.5 / tp)
axs[0].set_ylim(0, 5000)
axs[1].set_ylim(0, 500)
axs[2].set_ylim(0, 15)


fig.subplots_adjust(
    top=0.825, bottom=0.15, left=0.175, right=0.935, hspace=0.22, wspace=0.2
)

## access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

correlated_legend_box = mpatches.Patch(
    lw=lw * 0.75, hatch="//////", facecolor="white", edgecolor="k", alpha=0.75
)
uncorrelated_legend_box = mpatches.Patch(lw=lw * 0.75, color="k", alpha=0.75)
dash_line = Line2D([0], [0], ls="--", color="k", lw=lw * 0.75)
solid_line = Line2D([0], [0], ls=":", color="k", lw=lw * 0.75)


## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r"Correlated data", r"Uncorrelated data (shuffled)", r"$\epsilon_{sim}$"]
fig.legend(handles, labels, loc="upper right", fontsize=fontsize - 3, ncol=2)


plt.show()

# %%
## Locate value of N when difference between the median of molecular brightness distribution from correlated data and the simulated value of epsilon is equal to 5%

d1 = 0
delta_diferences_index = 0
while d1 < len(Epsilon_distributions_correlated_values):
    if abs(MEDIAN_epsilon_correlated[d1] - Epsilon_simulated) < 0.050 * max(
        abs(MEDIAN_epsilon_correlated - Epsilon_simulated)
    ):
        delta_diferences_index = d1
        d1 = len(Epsilon_distributions_correlated_values)
    d1 += 1

## Locate value of N for molecular brightness distribution from UNcorrelated data present the same MAD as the correlated epsilon distribution whichs median presents 5% to simulated epsilon value
delta_diferences_UN_corr_index = min(
    range(len(MAD_epsilon_UN_correlated)),
    key=lambda i: abs(
        MAD_epsilon_UN_correlated[i] - MAD_epsilon_correlated[delta_diferences_index]
    ),
)

print(
    MEDIAN_epsilon_correlated[delta_diferences_index],
    "\n",
    MAD_epsilon_correlated[delta_diferences_index],
    "\n",
    MEDIAN_epsilon_UN_correlated[delta_diferences_UN_corr_index],
    "\n",
    MAD_epsilon_UN_correlated[delta_diferences_UN_corr_index],
)

print(N[delta_diferences_index])
print(N[delta_diferences_UN_corr_index])


# %%
## Plot Figure 2
## Median of distributions as functions of number of data point use to calculate each entry

fig, (ax0, ax1, ax2) = plt.subplots(
    3, 1, figsize=(3, 3), sharex=True, sharey=False, dpi=300
)

letters = ["A", "B", "C"]
ax0.text(2.5e7, 2.4 / tp, letters[0], fontsize=11)
ax1.text(2.5e7, 2.4 / tp, letters[1], fontsize=11)
ax2.text(2.5e7, 0.8 / tp, letters[2], fontsize=11)

ax0.semilogx(N, MEDIAN_I_mean_correlated, "-", color="#0267C1", lw=lw)
ax0.plot(N, MEDIAN_I_mean_UN_correlated, "-", color="#2A0C4E", lw=lw)

ax0.fill_between(
    N,
    MEDIAN_I_mean_correlated - MAD_I_mean_correlated,
    MEDIAN_I_mean_correlated + MAD_I_mean_correlated,
    color="none",
    hatch="//////",
    alpha=0.5,
    edgecolor="#0267C1",
    label="I correlated data",
)

ax0.fill_between(
    N,
    MEDIAN_I_mean_UN_correlated - MAD_I_mean_UN_correlated,
    MEDIAN_I_mean_UN_correlated + MAD_I_mean_UN_correlated,
    color="#2A0C4E",
    alpha=0.3,
    label="I uncorrelated data",
)


ax1.semilogx(N, MEDIAN_VAR_correlated, "-", color="salmon", lw=lw)
ax1.plot(N, MEDIAN_VAR_UN_correlated, "-", color="goldenrod", lw=lw)

ax1.fill_between(
    N,
    MEDIAN_VAR_correlated - MAD_VAR_correlated,
    MEDIAN_VAR_correlated + MAD_VAR_correlated,
    color="none",
    hatch="//////",
    alpha=0.3,
    edgecolor="salmon",
    label="VAR correlated data",
)
ax1.fill_between(
    N,
    MEDIAN_VAR_UN_correlated - MAD_VAR_UN_correlated,
    MEDIAN_VAR_UN_correlated + MAD_VAR_UN_correlated,
    color="goldenrod",
    alpha=0.3,
    label="VAR uncorrelated data",
)


ax2.semilogx(N, MEDIAN_epsilon_correlated, "-", color="lightgreen", lw=lw)
ax2.plot(N, MEDIAN_epsilon_UN_correlated, "-", color="olivedrab", lw=lw)

ax2.fill_between(
    N,
    MEDIAN_epsilon_correlated - MAD_epsilon_correlated,
    MEDIAN_epsilon_correlated + MAD_epsilon_correlated,
    color="none",
    hatch="//////",
    alpha=0.5,
    edgecolor="lightgreen",
    label="B correlated data",
)
ax2.fill_between(
    N,
    MEDIAN_epsilon_UN_correlated - MAD_epsilon_UN_correlated,
    MEDIAN_epsilon_UN_correlated + MAD_epsilon_UN_correlated,
    color="olivedrab",
    alpha=0.3,
    label="B uncorrelated data",
)

ax2.axhline(
    Epsilon_simulated, 0, N[-1] * 1.25, ls="--", color="k", lw=lw * 0.75, alpha=0.75
)

head_length = 0.2
head_width = 0.05

# Top right arrow
arrowprops = dict(
    arrowstyle="->, head_length=%s, head_width=%s" % (head_length, head_width),
    facecolor="lightgreen",
    edgecolor="lightgreen",
)
ax2.annotate(
    "",
    xy=(
        N[delta_diferences_index],
        (MEDIAN_epsilon_correlated + MAD_epsilon_correlated * 0.05)[
            delta_diferences_index
        ],
    ),
    xytext=(
        N[delta_diferences_index],
        (MEDIAN_epsilon_correlated + MAD_epsilon_correlated)[delta_diferences_index]
        * 1.65,
    ),
    arrowprops=arrowprops,
)

# Bottom right arrow
arrowprops = dict(
    arrowstyle="<-, head_length=%s, head_width=%s" % (head_length, head_width),
    facecolor="lightgreen",
    edgecolor="lightgreen",
)
ax2.annotate(
    "",
    xy=(N[delta_diferences_index], 0.1),
    xytext=(
        N[delta_diferences_index],
        (MEDIAN_epsilon_correlated - MAD_epsilon_correlated)[delta_diferences_index]
        * 1.15,
    ),
    arrowprops=arrowprops,
)

# Top left arrow
arrowprops = dict(
    arrowstyle="->, head_length=%s, head_width=%s" % (head_length, head_width),
    facecolor="olivedrab",
    edgecolor="olivedrab",
)
ax2.annotate(
    "",
    xy=(
        N[delta_diferences_UN_corr_index],
        (MEDIAN_epsilon_correlated + MAD_epsilon_correlated * 0.05)[
            delta_diferences_index
        ],
    ),
    xytext=(
        N[delta_diferences_UN_corr_index],
        (MEDIAN_epsilon_correlated + MAD_epsilon_correlated)[delta_diferences_index]
        * 1.65,
    ),
    arrowprops=arrowprops,
)

# Bottom left arrow
arrowprops = dict(
    arrowstyle="<-, head_length=%s, head_width=%s" % (head_length, head_width),
    facecolor="olivedrab",
    edgecolor="olivedrab",
)
ax2.annotate(
    "",
    xy=(N[delta_diferences_UN_corr_index], 0.1),
    xytext=(
        N[delta_diferences_UN_corr_index],
        (MEDIAN_epsilon_correlated - MAD_epsilon_correlated)[delta_diferences_index]
        * 1.15,
    ),
    arrowprops=arrowprops,
)


ax0.set_ylabel(r"$<I>_{median}$", fontsize=fontsize - 2, labelpad=labelpad)
ax1.set_ylabel("$\sigma^2_{median}$", fontsize=fontsize - 2, labelpad=labelpad)
ax2.set_ylabel(r"$\epsilon_{median}$ ", fontsize=fontsize - 2, labelpad=labelpad)
plt.xlabel("$N$", fontsize=fontsize - 2, labelpad=labelpad)

ax0.tick_params(which="minor", length=1.25, width=1)
ax0.tick_params(which="major", length=2.25, width=1)
ax1.tick_params(which="minor", length=1.25, width=1)
ax1.tick_params(which="major", length=2.25, width=1)
ax2.tick_params(which="minor", length=1.25, width=1)
ax2.tick_params(which="major", length=2.25, width=1)

ax0.tick_params(axis="both", labelsize=fontsize - 3)
ax1.tick_params(axis="both", labelsize=fontsize - 3)
ax2.tick_params(axis="both", labelsize=fontsize - 3)

ax2.ticklabel_format(
    style="sci", axis="y", scilimits=(6, 6)
)  # Apply scientific notation to y-axis


ax0.ticklabel_format(
    style="scientific", axis="y", scilimits=(6, 6), useMathText=True
)  # Apply scientific notation to y-axis
ax0.yaxis.get_offset_text().set_fontsize(fontsize - 5)

ax1.ticklabel_format(
    style="scientific", axis="y", scilimits=(6, 6), useMathText=True
)  # Apply scientific notation to y-axis
ax1.yaxis.get_offset_text().set_fontsize(fontsize - 5)

ax2.ticklabel_format(
    style="scientific", axis="y", scilimits=(6, 6), useMathText=True
)  # Apply scientific notation to y-axis
ax2.yaxis.get_offset_text().set_fontsize(fontsize - 5)

ax0.set_ylim(1.3 / tp, 2.55 / tp)
ax1.set_ylim(1.3 / tp, 2.55 / tp)
ax2.set_ylim(-0.45 / tp, 0.95 / tp)

## affect borders lines plot making it thinners
for axis in ["top", "bottom", "left", "right"]:
    ax0.spines[axis].set_linewidth(0.6)
    ax1.spines[axis].set_linewidth(0.6)
    ax2.spines[axis].set_linewidth(0.6)

fig.subplots_adjust(
    top=0.790, bottom=0.13, left=0.165, right=0.935, hspace=0.325, wspace=0.2
)

## access legend objects automatically created from data
handles, labels = plt.gca().get_legend_handles_labels()

correlated_legend_box = mpatches.Patch(
    lw=lw * 0.75, hatch="//////", facecolor="white", edgecolor="k", alpha=0.75
)
uncorrelated_legend_box = mpatches.Patch(lw=lw * 0.75, color="k", alpha=0.75)
dash_line = Line2D([0], [0], ls="--", color="k", lw=lw * 0.75)

## add manual symbols to auto legend
handles = [correlated_legend_box, uncorrelated_legend_box, dash_line]
labels = [r"Correlated data", r"Uncorrelated data (shuffled)", r"$\epsilon_{sim}$"]
fig.legend(handles, labels, loc="upper right", fontsize=fontsize - 3, ncol=2)

plt.show()
