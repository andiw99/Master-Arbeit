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
    root = "../../Generated content/Silicon/Amplitude/50000"
    #root = "../../Generated content/Silicon Test/Anisotrop Antisymmetric Rectangular Subsystems Small2"

    plot_root = os.path.join(root, "plots/")

    config = {"nr_of_meshs": 4,
              "cell_L": 500,
              "cell_nr": 0,
              "chess_trafo": -1,
              "nr_colorbar_ticks": 7,
              "angle": 2,
              "subsystem": 0}

    folders = list_folders_and_subfolders(root)
    print(folders)
    for folder in folders:
        filepath = find_first_csv_file(folder)
        if filepath:
            print("reading: ", filepath)
            para_filepath = os.path.splitext(filepath)[0] + ".txt"
            fig, axes, name = plot_multiple_times(filepath, config)
            create_directory_if_not_exists(root + "/plots/")
            fig.savefig(root + "/plots/" + name, format="png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
