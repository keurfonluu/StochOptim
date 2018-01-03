# -*- coding: utf-8 -*-

"""
Author: Keurfon Luu <keurfon.luu@mines-paristech.fr>
License: MIT
"""

import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt


def fitness(Z, method = "mean", threshold = 0.):
    if method == "mean":
        return np.mean(Z, axis = -1) 
    elif method == "median":
        return np.median(Z, axis = -1)
    elif method == "min":
        return np.min(Z, axis = -1)
    elif method == "max":
        return np.max(Z, axis = -1)
    elif method == "sr":
        return np.sum(Z > threshold, axis = -1)
    else:
        raise ValueError("Error: unknown method '%s'" % method)


if __name__ == "__main__":
    # Parameters
    dtype = "float64"
    solver = "cpso"
    func = "rastrigin"
    nx, ny = 50, 50
    max_run = 50
    n_dim = 5
    indir = "../save"
    Z_file = "%s/%s_%s_%s_Z.bin" % (indir, solver, func, n_dim)
    xaxis_file = "%s/%s_%s_%s_xaxis.bin" % (indir, solver, func, n_dim)
    yaxis_file = "%s/%s_%s_%s_yaxis.bin" % (indir, solver, func, n_dim)
    
    # Load files
    Z = np.fromfile(Z_file, dtype = dtype).reshape((nx, ny, max_run), order = "F")
    xaxis = np.fromfile(xaxis_file, dtype = dtype)
    yaxis = np.fromfile(yaxis_file, dtype = dtype)
    
    # Clip Z
    fit = fitness(Z, method = "sr", threshold = 1.)
    fit = fit / max_run * 100

    # Plot
#    sns.set_style("ticks")
    fig = plt.figure(figsize = (8, 8), facecolor = "white")
    ax1 = fig.add_subplot(1, 1, 1)
    fig.patch.set_alpha(0.)
    c = ax1.contourf(yaxis, xaxis, fit, 100, cmap = "viridis_r", vmin = 0, vmax = 100)
#    cb = fig.colorbar(c)
    
    # Save figures
#    ax1.xaxis.set_visible(False)
#    ax1.yaxis.set_visible(False)
#    fig.tight_layout()
#    fig.savefig("%s/%s_%s_%s_Z.png" % (indir, solver, func, n_dim), dpi = 150)