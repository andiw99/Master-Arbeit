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


def plot_colormesh(df, fig=None, ax=None, title=None, proj=False, p=True):
    """
    plots the dataframe as a colormesh
    :param df:
    :return:
    """
    # actually we woulndt even have to read in all values since we only
    # plot the last one, but for now thats okay
    # We have to construct a n x n matrix out of the values that we extract
    try:
        # this is if I am saving the impuls values
        if p:
            n = (int(np.sqrt(df.shape[1] / 2 - 1)))
        else:
            n = int(np.sqrt(df.shape[1] - 1))
            print(n)
        z_values = df.iloc[-1, 1:n*n+1]
    except IndexError:
        # catch case that it is 1D-array
        if p:
            n = (int(np.sqrt(df.shape[0] / 2 - 1)))
        else:
            n = int(np.sqrt(df.shape[0] -1))
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
                           vmax=np.max(np.abs(z_values)), cmap="greys")
        fig.colorbar(cf, ax=ax)
        plt.show()
    else:
        cf = ax.pcolormesh(x, y, z_values,  vmin=-np.max(np.abs(z_values)),
                           vmax=np.max(np.abs(z_values)))
        ax.set_title(title)
        fig.colorbar(cf, ax=ax)


def plot_multiple_times(df, paras, n, proj=False, storage_root="plots/", p=True):
    # find out number of rows
    nr_rows = df.shape[0]
    # equidistant row numbers to use
    rows = np.linspace(0, nr_rows-1, n, endpoint=True)
    # Select the rows with the row equidistant row numbers
    df_rows = df.iloc[rows]
    # create fig and axes
    fig, axes = plt.subplots(int(np.sqrt(n)), int(np.sqrt(n)), figsize =[12, 10])
    i = 0
    for axs in axes:
        for ax in axs:
            # plot with one row out of the row numbers
            plot_colormesh(
                df_rows.iloc[i], fig, ax,
                title=f"t = {df_rows.iloc[i, 0]:.2f}, T = {df_rows.iloc[i, -1]}",
                proj=proj, p=p)
            i +=1
    # insert parameters
    textstr = ''
    for para in paras:
        textstr += para + "=" + str(paras[para]) + "\n"
    textstr = textstr[:-1]
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.02, 0.8, textstr, fontsize=14, bbox=props)
    plt.tight_layout()
    # if you want to save the pics somewhere else
    try:
        plt.savefig(storage_root + plot_name_paras(paras), format="png")
    except FileNotFoundError:
        os.makedirs(storage_root)
        plt.savefig(storage_root + plot_name_paras(paras), format="png")
    plt.show()


def plot_name_paras(paras):
    fname = ""
    for key in paras.keys():
        fname += key + "=" + str(paras[key])
    return fname



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
        # and check if it is a csv
        if os.path.isfile(f):
            if os.path.splitext(f)[1] == ".csv":
                if not f in old_files:
                    filepaths.append(f)
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


def read_parameters_txt(filepath):
    df = pd.read_csv(filepath, delimiter=",", header=None, index_col=0)
    para_set = {}
    for label in df.index:
        para_set[label] = df.loc[label, 1]
    return para_set


def read_multiple_parameters(filepaths, nr_parameters=8, txt=True):
    parameters = []
    for filepath in filepaths:
        # filepaths are the paths with .csv extension since .txt files are
        # ignored
        if txt:
            filepath = os.path.splitext(filepath)[0] + ".txt"
            para_set = read_parameters_txt(filepath)
        else:
            para_set = read_parameters(filepath, nr_parameters=nr_parameters)
        parameters.append(para_set)

    return parameters


def main():
    print("This file is made to import, not to execute")

if __name__ == "__main__":
    main()
