# -*- coding: utf-8 -*-

"""
Author: Keurfon Luu <keurfon.luu@mines-paristech.fr>
License: MIT
"""

import numpy as np
import seaborn as sns


if __name__ == "__main__":
    # Parameters
    dtype = "float64"
    max_iter, n_dim = 10000, 2
    models_file = "../output/models_%d_%d.bin" % (max_iter, n_dim)
    energy_file = "../output/energy_%d.bin" % max_iter
    
    # Load files
    models = np.fromfile(models_file, dtype = dtype).reshape((max_iter, n_dim), order = "F")
    energy = np.fromfile(energy_file, dtype = dtype)
    
    # Plot
    sns.set_style("ticks")
    sns.jointplot(models[:,0], models[:,1], kind = "kde")