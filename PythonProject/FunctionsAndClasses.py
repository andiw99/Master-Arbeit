import csv
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib

# import matplotlib; matplotlib.use("TkAgg")


def read_csv(filepath, nrows=None):
    df = pd.read_csv(filepath, header=None, index_col=None)
    return df.iloc[:nrows]


def read_multiple_csv(filepaths, nrows=None):
    trajectories = []
    for filepath in filepaths:
        df = read_csv(filepath, nrows=nrows)
        trajectories.append(df)

    return trajectories


def read_struct_func(filepath):
    df = pd.read_csv(filepath, delimiter=",", index_col=False)
    return df

def plot_colormesh(df, fig=None, ax=None, title=None, proj=False, p=True, beta=2, chess_board=False):
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
    if(chess_board):
        z_values = chess_board_trafo(z_values)


    x = np.arange(0.5, n+1, 1)
    y = np.arange(0.5, n + 1, 1)
    if fig is None:
        fig, ax = plt.subplots()
        ax.set_title(title)
        vmin = -np.max(np.abs(z_values))
        vmax = np.max(np.abs(z_values))
        # find out how many multiples of the minimum we need
        min_pos = np.sqrt(beta/2)
        print(min_pos)
        nr_ticks = vmax // min_pos + 1
        tick_labels = np.arange(-nr_ticks, nr_ticks + 1)
        # calculate tick positions
        tick_positions = tick_labels * min_pos
        # make labels to strings
        tick_labels = [str(label) for label in tick_labels]
        cf = ax.pcolormesh(x, y, z_values, vmin=vmin,
                           vmax=vmax)
        print(tick_labels)
        cbar = fig.colorbar(cf, ax=ax, ticks=tick_positions)
        cbar.ax.set_yticklabels(tick_labels)
        plt.show()
    else:
        ax.set_title(title)
        vmin = -np.max(np.abs(z_values))
        vmax = np.max(np.abs(z_values))
        # find out how many multiples of the minimum we need
        min_pos = np.sqrt(beta/2)
        max_tick_nr = vmax // min_pos + 1
        # always 5 ticks
        tick_labels = np.int32(np.linspace(-max_tick_nr, max_tick_nr, 7))
        # could be to long

        # calculate tick positions
        tick_positions = tick_labels * min_pos
        # make labels to strings

        tick_labels = [str(int(label)) for label in tick_labels]
        cf = ax.pcolormesh(x, y, z_values / min_pos, vmin=-1, vmax=1, cmap="copper")
        cbar = fig.colorbar(cf, ax=ax, ticks=tick_positions)
        cbar.ax.set_yticklabels(tick_labels)
        """
        cf = ax.pcolormesh(x, y, z_values,  vmin=-np.max(np.abs(z_values)),
                           vmax=np.max(np.abs(z_values)))
        ax.set_title(title)
        fig.colorbar(cf, ax=ax)
        """



def plot_multiple_times(df, paras, n, proj=False, storage_root="plots/", p=True, chess_board=False, show=False, name="", v1=1):
    # find out number of rows
    nr_rows = df.shape[0]
    # equidistant row numbers to use
    rows = np.linspace(0, nr_rows-1, n, endpoint=True)
    # Select the rows with the row equidistant row numbers
    df_rows = df.iloc[rows]
    # read beta
    beta = paras["beta"]
    # create fig and axes
    fig, axes = plt.subplots(int(np.sqrt(n)), int(np.sqrt(n)), figsize =[12, 10])
    i = 0
    for axs in axes:
        for ax in axs:
            # plot with one row out of the row numbers
            try:
                plot_colormesh(
                df_rows.iloc[i], fig, ax,
                title=f"t = {df_rows.iloc[i, 0]:.2f}, T = {df_rows.iloc[i, v1 * -1]}",
                proj=proj, p=p, beta=beta, chess_board=chess_board)
            except ValueError:
                plot_colormesh(
                df_rows.iloc[i][1: -1], fig, ax,
                title=f"t = {df_rows.iloc[i, 1]:.2f}, T = {df_rows.iloc[i, -1]}",
                proj=proj, p=p, beta=beta, chess_board=chess_board)

            i +=1
            # scale colormap in multiples of the position of the minimum

    # insert parameters
    textstr = ''
    wanted_paras= ["dt", "J", "alpha", "beta", "eta", "tau"]
    for para in wanted_paras:
        try:
            textstr += para + "=" + str(paras[para]) + "\n"
        except:
            pass
    textstr = textstr[:-1]
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.02, 0.8, textstr, fontsize=14, bbox=props)
    plt.tight_layout()
    # if you want to save the pics somewhere else

    if name == "":
        name = plot_name_paras(paras)

    try:
        plt.savefig(storage_root + name, format="png")
    except FileNotFoundError:
        os.makedirs(storage_root)
        plt.savefig(storage_root + name, format="png")
    if show:
        plt.show()


def make_dir(path):
    try:
        os.makedirs(path)
    except FileExistsError:
        pass


def save_plot(path, name, format="png"):
    make_dir(path)
    plt.savefig(path + name, format=format)

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


def chess_board_trafo(x):
    """
    Just flips the sign of every second value of x
    :param x: array
    :return: void
    """
    # in place trafo vermutlich langsamer als neues array aber egal

    for i in range(x.shape[0] // 2):
        for j in range(x.shape[1] // 2):
            # we iterate over all even values with even i, j
            x[2 * i][2 * j] *= (-1)
            # we iterate over lattice sites with odd indices
            x[2 * i + 1][2 * j + 1] *= (-1)

    # Lol we have to return this since it is not in place
    return x


def pd_chess_board_trafo(x):
    for i in range(x.shape[0]//2):
        for j in range(x.shape[1]//2):
            # we iterate over all even values with even i, j
            x.iloc[i, j] *= (-1)
            # we iterate over lattice sites with odd indices
            x.iloc[2 * i + 1, 2 * j + 1] += (-1)

def linear_fit(x, m, a):
    return a + m*x
def poly(x, exp, ampl):
    return ampl * x ** exp

def linear_corr(x, m, a, b):
    return a + m * x + b * np.exp(-x)
def main():
    print("This file is made to import, not to execute")

if __name__ == "__main__":
    main()
