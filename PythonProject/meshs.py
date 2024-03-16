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
    root = "../../Generated content/Final/Amplitude/J_J=30/Lx_Ly=0.25/Amplitude/128/"
    root = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=2-3/"
    root = "../../Generated content/Final/Amplitude/J_J=60/final/Amplitude/4096"
    fig_format = "png"
    dpi = 200
    #root = "../../Generated content/Silicon Test/Anisotrop Antisymmetric Rectangular Subsystems Small2"

    plot_root = os.path.join(root, "plots/")

    config = {"nr_of_meshs": 2,
              "cell_L": 1024,
              "cell_nr": 0,
              "chess_trafo": -1,
              "nr_colorbar_ticks": 7,
              "angle": 3,       # 2 is -np.pi / 2 to np.pi / 2 , 3 is with mean of the system
              "subsystem": 1,
              "colormap": 'viridis'} # 'PiYG', 'viridis'

    folders = list_folders_and_subfolders(root)
    print(folders)
    # check which usecase, every csv in the given root or one for every folder in root:
    if folders:
        for folder in folders:
            filepath = find_first_csv_file(folder)
            if filepath:
                print("reading: ", filepath)
                para_filepath = os.path.splitext(filepath)[0] + ".txt"
                try:
                    fig, axes, name = plot_multiple_times(filepath, config)
                    create_directory_if_not_exists(root + "/plots/")
                    fig.savefig(root + "/plots/" + name + f".{fig_format}", format=fig_format, dpi=dpi)
                except ValueError as e:
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    print(e)
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        for file in os.listdir(root):
            if file.endswith(".csv"):
                filepath = os.path.join(root, file)
                print("reading: ", filepath)
                fig, axes, name = plot_multiple_times(filepath, config)
                create_directory_if_not_exists(root + "/plots/")
                fig.savefig(root + "/plots/" + file[0] + "-" + name + f".{fig_format}", format=fig_format, dpi=dpi)

    plt.show()

if __name__ == "__main__":
    main()
