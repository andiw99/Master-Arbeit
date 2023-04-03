import csv
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib

# import matplotlib; matplotlib.use("TkAgg")

def read_csv(filepath, nrows=None):
    df = pd.read_csv(filepath, header=None, index_col=0)
    return df.iloc[:nrows]


def read_multiple_csv(filepaths, nrows=None):
    trajectories = []
    for filepath in filepaths:
        df = read_csv(filepath, nrows=nrows)
        trajectories.append(df)

    return trajectories


def plot_colormesh(df, fig=None, ax=None, title=None, proj=False):
    """
    plots the dataframe as a colormesh
    :param df:
    :return:
    """
    # actually we woulndt even have to read in all values since we only
    # plot the last one, but for now thats okay
    # We have to construct a n x n matrix out of the values that we extract
    try:
        n = (int(np.sqrt(df.shape[1] / 2 - 1)))
        z_values = df.iloc[-1, 1:n*n+1]
    except IndexError:
        # catch case that it is 1D-array
        n = (int(np.sqrt(df.shape[0] / 2 - 1)))
        z_values = df.iloc[1:n*n+1]
    if proj:
        z_values = np.sign(z_values)
    # for now i just want to know in which minimum i fall so i make this
    # projection
    z_values = np.array(z_values, dtype=float).reshape((n, n))

    x = np.arange(0.5, n+1, 1)
    y = np.arange(0.5, n + 1, 1)
    if fig is None:
        fig, ax = plt.subplots()
        ax.set_title(title)
        cf = ax.pcolormesh(x, y, z_values, vmin=-np.max(np.abs(z_values)),
                           vmax=np.max(np.abs(z_values)))
        fig.colorbar(cf, ax=ax)
        plt.show()
    else:
        cf = ax.pcolormesh(x, y, z_values,  vmin=-np.max(np.abs(z_values)),
                           vmax=np.max(np.abs(z_values)))
        ax.set_title(title)
        fig.colorbar(cf, ax=ax)



def plot_multiple_times(df, paras, n, proj=False):
    # find out number of rows
    nr_rows = df.shape[0]

    # equidistant row numbers to use
    rows = np.linspace(0, nr_rows-1, n, endpoint=True)
    df_rows = df.iloc[rows]
    # create fig and axes
    fig, axes = plt.subplots(int(np.sqrt(n)), int(np.sqrt(n)), figsize =[12, 10])
    i = 0
    for axs in axes:
        for ax in axs:
            plot_colormesh(df_rows.iloc[i], fig, ax,
                           title=f"t = {df_rows.iloc[i, 0]}", proj=proj)
            i +=1
    # insert parameters
    textstr = ''
    for para in paras:
        textstr += para + "=" + str(paras[para]) + "\n"
    textstr = textstr[:-1]
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.02, 0.8, textstr, fontsize=14, bbox=props)
    plt.tight_layout()
    plt.show()


def get_plotted_files(root):
    try:
        with open(root + "plotted_files.txt") as f:
            # readlines returns list where each line is an item in the list
            plotted_files = f.read().splitlines()
        return plotted_files
    # if its a new root and we haven't plotted anything yet:
    except FileNotFoundError:
        # write the plotted_files.txt in itself lol
        with open(root + "plotted_files.txt", "a") as f:
            f.write(root + "plotted_files.txt\n")
        return [root + "plotted_files.txt"]


def new_files_in_dir(cur_dir, root, old_files=None, plot_all=False):
    filepaths = []
    if not old_files:
        old_files = get_plotted_files(root)
    if plot_all:
        old_files=[root + "plotted_files.txt"]
    for filename in os.listdir(cur_dir):
        f = os.path.join(cur_dir, filename)
        # check if f is file and if it wasnt already plotted
        if os.path.isfile(f):
            if not f in old_files:
                filepaths.append(f)
                plotted_files = open(root + "/plotted_files.txt", "a")
                plotted_files.write(f + "\n")
        elif os.path.isdir(f):
            filepaths += new_files_in_dir(f, root, old_files, plot_all)

    return filepaths


def read_parameters(filepath, nr_parameters):
    # TODO inefficient since i need to read the file two times but
    # therefore better for flexibility
    df = pd.read_csv(filepath, header=None, index_col=0)
    df = df.iloc[-nr_parameters:, 0:1]
    para_set = {}
    for label in df.index:
        para_set[label] = df.loc[label, 1]
    return para_set

def read_multiple_parameters(filepaths, nr_parameters):
    parameters = []
    for filepath in filepaths:
        para_set = read_parameters(filepath, nr_parameters=nr_parameters)
        parameters.append(para_set)

    return parameters


def main():
    filepath = "../LearningProject/Brownian Motion.csv"
    root = "../LearningProject/Brownian Motion Full Interaction/"
    proj = False
    nr_parameters=8
    nr_plots_per_dimension = 3

    # TODO very unefficient to read in all trajectories first and then plot
    filepaths = new_files_in_dir(root, root, plot_all=True)
    trajectories = read_multiple_csv(filepaths, nrows=-nr_parameters)
    parameters = read_multiple_parameters(filepaths, nr_parameters)
    # plot_colormesh(trajectories[0])
    for traj_param_tuple in zip(trajectories, parameters):
        plot_multiple_times(traj_param_tuple[0], traj_param_tuple[1],
                            nr_plots_per_dimension ** 2, proj)

if __name__ == "__main__":
    main()
