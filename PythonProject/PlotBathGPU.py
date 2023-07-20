import csv
import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib
from FunctionsAndClasses import *

# import matplotlib; matplotlib.use("TkAgg")

def main():
    # some minor changes because of the file that is now generated by gpu and which doesn't save
    # the impulses

    root = "../../Generated content/Fit testing/"
    if len(sys.argv) >= 2:
        # Access and print the command-line arguments
        root = sys.argv[1]
    else:
        print("Please provide file path")
    proj = False
    nr_parameters = 8
    nr_plots_per_dimension = 4
    plot_all = False
    txt = True                  # whether the parameters are saved in a txt file
    chess_trafo = False
    v1 = -1
    plot_root = os.path.join(root, "plots/")

    # TODO very unefficient to read in all trajectories first and then plot
    filepaths = new_files_in_dir(root, root, plot_all=plot_all)
    print(filepaths)
    #trajectories = read_multiple_csv(filepaths)
    #parameters = read_multiple_parameters(filepaths, txt=txt)
    # plot_colormesh(trajectories[0])
    for filepath in filepaths[::-1]:
        # plot only if it is the first csv?
        # TODO kind of fishy just asking for 0.csv since if there is no 0.csv it wont plot anything
        if os.path.basename(filepath) == "0.csv":
            para_filepath = os.path.splitext(filepath)[0] + ".txt"
            plot_multiple_times(read_csv(filepath), read_parameters_txt(para_filepath),
                                nr_plots_per_dimension ** 2, proj,
                                storage_root=plot_root, p=False, show=False, chess_board=chess_trafo, v1=v1)
        # We just plot twice, one time for the large plot folder which lets us do easy comparisons and one
        # for all the repetitions
        if os.path.splitext(filepath)[1] == ".csv":
            print(filepath + "   done")
            para_filepath = os.path.splitext(filepath)[0] + ".txt"
            # storage path
            storage_path = os.path.join(os.path.dirname(filepath), "plots/")
            # name
            name = os.path.splitext(os.path.basename(filepath))[0] + ".png"
            plot_multiple_times(read_csv(filepath), read_parameters_txt(para_filepath),
                                nr_plots_per_dimension ** 2, proj,
                                storage_root=storage_path, p=False, name=name, v1=v1)
    # at the end (no error) we write the plotted files into the tracking file
    plotted_files = open(root + "/plotted_files.txt", "a")
    for f in filepaths:
        plotted_files.write(f + "\n")
    # newline for visibility
    plotted_files.write("\n")

if __name__ == "__main__":
    main()
