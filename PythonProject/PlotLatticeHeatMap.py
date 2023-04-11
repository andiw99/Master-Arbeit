import csv
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib
from FunctionsAndClasses import *

# import matplotlib; matplotlib.use("TkAgg")

def main():
    filepath = "../LearningProject/Brownian Motion.csv"
    root = "../../Generated content/Domain Size Test/"
    proj = False
    nr_parameters = 8
    nr_plots_per_dimension = 3
    plot_all = False
    plot_root = "../../Generated content/Domain Size Test/plots/"

    # TODO very unefficient to read in all trajectories first and then plot
    filepaths = new_files_in_dir(root, root, plot_all=plot_all)
    print(filepaths)
    trajectories = read_multiple_csv(filepaths, nrows=-nr_parameters)
    parameters = read_multiple_parameters(filepaths, nr_parameters)
    # plot_colormesh(trajectories[0])
    for traj_param_tuple in zip(trajectories, parameters):
        plot_multiple_times(traj_param_tuple[0], traj_param_tuple[1],
                            nr_plots_per_dimension ** 2, proj,
                            storage_root=plot_root)
    # at the end (no error) we write the plotted files into the tracking file
    plotted_files = open(root + "/plotted_files.txt", "a")
    for f in filepaths:
        plotted_files.write(f + "\n")
    # newline for visibility
    plotted_files.write("\n")

if __name__ == "__main__":
    main()
