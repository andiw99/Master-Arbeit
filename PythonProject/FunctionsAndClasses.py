import csv
import math
import os
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
from itertools import product, combinations
from IPython import embed
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.stats import linregress
import statsmodels.tsa.stattools as st
from pathlib import Path



# import matplotlib; matplotlib.use("TkAgg")

colors = ["#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356", "#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356"]
colors += colors + colors + colors + colors
colors += colors + colors + colors + colors
z_colors = ["#00305d", "#009de0", "#65b32e", "C1"]
markers = ["o", "s", "^", "v", "D", "p", "1", "2","*", "x", "+", "v", "^"]
blue_point_kwargs = {"linestyle": "None", "markerfacecolor": "none", "markeredgecolor": colors[0]}
blue_square_kwargs = {"linestyle": "None", "markerfacecolor": "none", "markeredgecolor": colors[0], "marker": "s"}
errorbar_kwargs = {"color": "black", "elinewidth":1,  "capsize": 2}
standard_figsize = (10, 4.8 / 6.4 * 10)

def get_point_kwargs_color(color, markeredgewidth=1):
    return {"linestyle": "None", "markerfacecolor": "none", "markeredgecolor": color, "markeredgewidth": markeredgewidth}

def read_csv(filepath, nrows=None):
    df = pd.read_csv(filepath, header=None, index_col=None)
    return df.iloc[:nrows]

def create_directory_if_not_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
def read_last_line(path):
    with open(path, 'rb') as f:
        try:  # catch OSError in case of a one line file
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
    print(last_line[-200:-2])

    return last_line[:-2]   # strip last comma

def write_lists_to_file(x, y, filename, header):
    with open(filename, 'w') as file:
        file.write(header)
        for x_val, y_val in zip(x, y):
            file.write(f'{x_val},{y_val}\n')
def string_to_array(input_string):
    # Split the input string by commas and convert each element to a float
    # input_string = input_string.strip()
    values = string_to_list(input_string)

    # Create a NumPy array from the list of float values
    array = np.array(values)

    return array

def string_to_list(input_string):
    return [float(x) for x in input_string.split(',')]

def read_multiple_csv(filepaths, nrows=None):
    trajectories = []
    for filepath in filepaths:
        df = read_csv(filepath, nrows=nrows)
        trajectories.append(df)

    return trajectories

def read_struct_func(filepath):
    df = pd.read_csv(filepath, delimiter=",", index_col=False)
    return df

def corr_scaling_right(T, Tc, nu, xi0):
    eps = (Tc - T) / Tc         # so negative temps are above Tc
    return xi0 / ((-eps) ** nu)

def corr_scaling_left(T, Tc, nu, xi0):
    eps = (Tc - T) / Tc         # so negative temps are above Tc
    return xi0 / (eps ** nu)

def read_large_df_array(filepath, rows=None, skiprows=0, sep=None):
    df = []
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i >= skiprows:
                if rows:
                    if i in rows:
                        if sep:
                            line = line.split(";")[1]
                        df.append(string_to_array(line[:-2]))
                else:
                    if sep:
                        line = line.split(";")[1]
                    df.append(string_to_array(line[:-2]))
    return df

def read_large_df(filepath, rows=None, skiprows=0, sep=None, cut_endline=2):
    df = []
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i >= skiprows:
                if rows:
                    if i in rows:
                        if sep:
                            line = line.split(";")[1]
                        df.append(string_to_list(line[:-2]))
                else:
                    if sep:
                        line = line.split(";")[1]
                    df.append(string_to_list(line[:-cut_endline]))
    return df

def read_large_df_with_times(filepath, rows=None, skiprows=0, sep=";", cut_endline=2, max_row=-1):
    df = []
    times = []
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i == max_row:
                break
            if i >= skiprows:
                if rows:
                    if i in rows:
                        t, line = line.split(sep)
                        df.append(string_to_list(line[:-2]))
                        times.append(float(t))
                else:
                    t, line = line.split(sep)
                    df.append(string_to_list(line[:-cut_endline]))
                    times.append(float(t))
    return df, times

def read_line_large_df_with_times(filepath, line_nr, rows=None, skiprows=0, sep=";", cut_endline=2):
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i == skiprows + line_nr:
                t, line = line.split(sep)
                line = string_to_list(line[:-2][:-cut_endline])
                return line, float(t)
    return [], None

def plot_multiple_times(filepath, config={"nr_of_meshs": 16, "cell_L": 128, "cell_nr": 1, "chess_trafo": 1, "nr_colorbar_ticks": 5}):

    para_filepath = os.path.splitext(filepath)[0] + ".txt"
    parameters = read_parameters_txt(para_filepath, skipfooter=1)

    # extract parameters out of config
    nr_of_meshs = config["nr_of_meshs"]
    cell_L = config["cell_L"]
    parameters["cell_L"] = cell_L
    with open(filepath, "rbU") as f:
        nr_rows = sum(1 for _ in f)
    print("nr of rows: ", nr_rows)
    Lx = int(parameters["dim_size_x"])
    Ly = int(parameters["dim_size_y"])
    # by now we always need subsystem because we only use the subsystem simulation
    subsystem_Lx = parameters["subsystem_Lx"]
    subsystem_Ly = parameters["subsystem_Ly"]
    x_y_factor = float(parameters["x_y_factor"])
    if nr_rows < nr_of_meshs:
        print(f"{nr_of_meshs} snapshots not available, using {nr_rows} instead.")
        nr_of_meshs = nr_rows
    # equidistant row numbers to use
    # Select the rows with the row equidistant row numbers
    rows = np.linspace(0, nr_rows - 1, nr_of_meshs, endpoint=True)
    rows = [int(row) for row in rows]
    df = read_large_df(filepath, rows)
    stretch = subsystem_Lx / subsystem_Ly
    print(stretch)
    if config["subsystem"]:
        stretch = np.minimum(parameters["subsystem_Lx"], cell_L) / parameters["subsystem_Ly"]
        print("stretch = " , stretch)
    if nr_of_meshs == 2:
        fig, axes = plt.subplots(2, 1, figsize=[2 + 10 * stretch/2, 10])
    else:
        fig, axes = plt.subplots(int(np.sqrt(nr_of_meshs)), int(np.sqrt(nr_of_meshs)), figsize=[2 + 10 * stretch, 10])
    i = 0
    if nr_of_meshs == 1:
        axes = [axes]
    for axs in (axes):
        if nr_of_meshs <= 2:
            axs = [axs]
        for ax in (axs):
            ind = i
            print("ind = ", ind)
            # read time and temperature and extract the values
            parameters["t"] = df[ind][0]
            parameters["T"] = df[ind][1]
            row = np.array(df[ind][2:])

            actual_cell_L = config["cell_L"]
            if config["cell_L"]:

                if config["subsystem"]:
                    print(np.sqrt(parameters["nr_subsystems"] * parameters["x_y_factor"]))
                    max_nr_subsystems_per_row = int(np.ceil(np.sqrt(parameters["nr_subsystems"] * parameters["x_y_factor"]) + 0.5))
                    max_nr_subsystems_per_col = int(np.ceil(parameters["nr_subsystems"] / max_nr_subsystems_per_row + 0.5))
                    nr_subsystems_per_row = np.minimum(int(config["cell_L"] / subsystem_Lx), max_nr_subsystems_per_row)

                    if nr_subsystems_per_row == 0:
                        nr_subsystems_per_row = 1
                        nr_subsystems_per_col = 1
                    else:
                        nr_subsystems_per_col = np.minimum(int(np.ceil(config["cell_L"] / subsystem_Ly + 0.5)), max_nr_subsystems_per_col)
                    print(nr_subsystems_per_col)
                    print(nr_subsystems_per_row)
                    nr_subsystems = nr_subsystems_per_row * nr_subsystems_per_col
                    # we now need to extract a rectangular cell that is subsystem_Ly high and nr_subsystems * subsystem_Lx wide
                    row = extract_rectangle_from_rectangle(row, config["cell_nr"], nr_subsystems, subsystem_Lx, subsystem_Ly, Lx, Ly)
                    if config["cell_L"] < subsystem_Lx:
                        print(row.shape)
                        row = extract_rectangle_from_rectangle(row.flatten(), 0, 1, config["cell_L"], subsystem_Ly, subsystem_Lx, subsystem_Ly)
                        print(row.shape)
                    else:
                        config["cell_L"] = nr_subsystems_per_row * subsystem_Lx
                else:
                    nr_cells = 1    # If not subsystems the nr of cells is always 1
                    row = extract_rectangle_from_rectangle(row, config["cell_nr"], nr_cells, subsystem_Lx, subsystem_Ly, Lx, Ly)
            if config["horizontal_center"]:
                print("centering horizontally...")
                row = center_horizontal(row)
            im = plot_rectangular_colormesh(ax, row, parameters, config)
            config["cell_L"] = actual_cell_L
            i += 1
    plt.tight_layout()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.04, 0.02, 0.7])
    nr_ticks = config["nr_colorbar_ticks"]
    if config["angle"] == 1:
        ticks = np.linspace(0, 2 * np.pi, nr_ticks, endpoint=True)
        tick_labels = np.linspace(0, 2 * np.pi, nr_ticks, endpoint=True)
        tick_labels = [f"{tick_label:.2f}" for tick_label in tick_labels]
    elif config["angle"] == 2:
        ticks = np.linspace(- np.pi / 2, np.pi / 2, nr_ticks, endpoint=True)
        tick_labels = np.linspace(-np.pi / 2, np.pi / 2, nr_ticks, endpoint=True)
        tick_labels = [f"{tick_label:.2f}" for tick_label in tick_labels]
    elif config["angle"] == 3:
        # this one calculates the position of the minimum and then sets those
        # plus 5% to be the limits
        # if i use the row here it should be the last row so the equilibrated one?
        min_pos = np.mean(np.abs(row))
        print("min_pos = ", min_pos)
        min_legend = - 1.15 * min_pos
        max_legend = 1.15 * min_pos
        ticks = np.linspace(min_legend, max_legend, nr_ticks, endpoint=True)
        tick_labels = np.linspace(min_legend, max_legend, nr_ticks,
                                  endpoint=True)
        tick_labels = [f"{tick_label:.2f}" for tick_label in tick_labels]
    else:
        well_pos = np.sqrt(parameters["beta"] / 2)
        ticks = np.linspace(- 1.5 * well_pos, 1.5 * well_pos, nr_ticks, endpoint=True)
        tick_labels = np.linspace(-1.5, 1.5, nr_ticks, endpoint=True)
        tick_labels = [str(tick_label) for tick_label in tick_labels]
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=ticks)
    cbar.ax.set_yticklabels(tick_labels)

    textstr = ''
    wanted_paras=["dt", "J", "alpha", "beta", "eta", "tau", "Jy", "lat_dim"]
    if config["subsystem"]:
        wanted_paras.append("subsystem_Lx")
    for para in wanted_paras:
        try:
            textstr += para + "=" + str(parameters[para]) + "\n"
        except:
            pass
    textstr = textstr[:-1]
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.83, 0.8, textstr, fontsize=16, bbox=props)
    plt.show()
    name = plot_name_paras(parameters)
    return fig, axes, name


def center_horizontal(lattice):
    # I think we just assume -pi/2 pi/2 interval here
    lattice += np.pi / 2
    lattice[lattice > np.pi/2] -= np.pi
    # is it really that easy?
    return lattice

def plot_colormesh(ax, row, parameters, config):
    title=f"t = {parameters['t']:.2f}, T = {parameters['T']}"
    ax.set_title(title)
    L = int(np.sqrt(row.size))
    print(row)
    row = row.reshape((L, L))
    chess_trafo = 0
    if config["chess_trafo"] == 1:
        chess_trafo = True
    elif config["chess_trafo"] == -1:
        # 'auto', so infering fom J
        print("J = ", parameters["J"])
        chess_trafo = parameters["J"] < 0
        print("auto chess trafo yields: chess_trafo=", chess_trafo)
    if chess_trafo:
        print("doing chess trafo")
        row = chess_board_trafo(row)
    print(row)
    well_pos = np.sqrt(parameters["beta"] / 2)
    cf = ax.pcolormesh(row, cmap="viridis_r", vmax=2 * well_pos, vmin=-2 * well_pos)

    return cf

def plot_rectangular_colormesh(ax, row, parameters, config):
    title=f"t = {parameters['t']:.2f}, T = {parameters['T']}"
    ax.set_title(title)
    Lx = int(parameters["dim_size_x"])
    Ly = int(parameters["dim_size_y"])
    J_para = np.abs(parameters["J"])
    J_perp = np.abs(parameters["Jy"])
    h = parameters["alpha"]
    try:
        p = parameters["p"]
        if p == 2:
            print("Guessing p = 2 is the wrong p, using p = 2.57")
            p = 2.57
    except KeyError:
        print("No p available in parameter file, using p = 2.57")
        p = 2.57
    subsystem_Lx = int(parameters["subsystem_Lx"])
    subsystem_Ly = int(parameters["subsystem_Ly"])
    cell_L = int(config["cell_L"])
    chess_trafo = 0
    if config["chess_trafo"] == 1:
        chess_trafo = True
    elif config["chess_trafo"] == -1:
        # 'auto', so infering fom J
        chess_trafo = parameters["J"] < 0
        print("auto chess trafo yields: chess_trafo=", chess_trafo)
    if chess_trafo:
        print("doing chess trafo")
        if config["subsystem"]:
            print("cell_L = ", cell_L, "  Lx = ", Lx)
            if cell_L < Lx:
                # sometimes it is cell_L, sometimes it is subsystem_Lx
                # I think we always want the smaller one?
                if cell_L > subsystem_Lx:
                    row = chess_board_trafo_rectangular_subsystems(row, subsystem_Lx, Ly)
                else:
                    row = chess_board_trafo_rectangular_subsystems(row, cell_L, Ly)
            else:
                print("cell_L > Lx")
                row = chess_board_trafo_rectangular_subsystems(row, Lx, Ly)
        else:
            row = chess_board_trafo(row)
    print(title)
    print(row)
    colormap = config["colormap"]
    if config["angle"] == 1:
        cf = ax.pcolormesh(row, cmap=colormap, vmax=2 * np.pi, vmin=0)
    elif config["angle"] == 4:
        cf = ax.pcolormesh(row, cmap=colormap, vmax=np.pi/2 / 1.2, vmin=-np.pi/2/ 1.2, rasterized=True)
    elif config["angle"] == 2:
        cf = ax.pcolormesh(row, cmap=colormap, vmax=np.pi / 2, vmin=-np.pi / 2, rasterized=True)
    elif config["angle"] == 3:
        # this should actually use the minimum of the system, but well I think
        # we can calculate it?
        # print("p = ", p)
        # print("J_para = ", J_para)
        # print("J_perp = ", J_perp)
        # print("h = ", h)
        theta_equil = get_equilibrium_position(J_para, J_perp, h, p)
        # print("equilibrium angle equation: ", equilibrium_angle_equation(theta_equil, J_para, J_perp, h, p))
        # print(theta_equil)
        cf = ax.pcolormesh(row, cmap=colormap, vmax= 1.15 * theta_equil,
                           vmin=-1.15 * theta_equil, rasterized=True, antialiased=True)
    else:
        well_pos = np.sqrt(parameters["beta"] / 2)
        cf = ax.pcolormesh(row, cmap=colormap, vmax=2 * well_pos, vmin=-2 * well_pos, linewidth=0)
    cf.set_edgecolor('face')

    return cf

def find_first_csv_file(folder_path):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".csv"):
            return os.path.join(folder_path, file_name)

    return None

def find_first_txt_file(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".txt"):
                return os.path.join(root, file)
    return None  # Return None if no .txt file is found

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
    wanted_paras=["cell_L","subsystem_Lx", "T", "J", "alpha", "beta", "eta", "tau", "Jy", "lat_dim"]
    for key in wanted_paras:
        try:
            if key == "T":
                fname += key + "=" + str(round(paras[key], 3))
            else:
                fname += key + "=" + str(round(paras[key], 2))
        except:
            pass
    return fname

def list_folders_and_subfolders(directory_path):
    folder_list = []

    for root, dirs, files in os.walk(directory_path):
        for dir_name in dirs:
            if (dir_name[0] != ".") & (dir_name != "plots"):
                folder_list.append(os.path.join(root, dir_name))

    return folder_list

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

def read_parameters_txt(filepath, skipfooter=1):
    print(filepath)
    if not filepath.endswith(".txt"):
        filepath = os.path.splitext(filepath)[0] + ".txt"
    df = pd.read_csv(filepath, delimiter=",", header=None, index_col=0)
    para_set = {}
    for label in df.index:
        try:
            para_set[label] = float(df.loc[label, 1])
        except ValueError:
            pass
        except TypeError:
            para_set[label] = float(df.loc[label, 1][-1])       # take the last one as it should be the most recent
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
    nr_rows = int(np.ceil(x.shape[0] / 2))
    nr_cols = int(np.ceil(x.shape[1] / 2))
    for i in range(nr_rows):
        for j in range(nr_cols):
            # we iterate over all even values with even i, j
            try:
                x[2 * i][2 * j] *= (-1)
                # we iterate over lattice sites with odd indices
                x[2 * i + 1][2 * j + 1] *= (-1)
            except IndexError:
                pass

    # Lol we have to return this since it is not in place
    return x

def chess_board_trafo_rectangular_subsystems(x, subsystem_Lx, subsystem_Ly):
    """
    We expect the ready 2D array?
    :param x:
    :return:
    """
    subsystems_per_row = x.shape[1] // subsystem_Lx
    subsystems_per_col = x.shape[0] // subsystem_Ly
    for col in range(subsystems_per_row):
        for row in range(subsystems_per_col):
            # extract i-th subsystem
            #embed()
            subsystem = x[row * subsystem_Ly : (row+1) * subsystem_Ly, col * subsystem_Lx : (col+1) * subsystem_Lx]
            subsystem = chess_board_trafo(subsystem)
            #if ((subsystem_Lx % 2 != 0) & (row % 2 != 0)) | ((subsystem_Ly % 2 != 0) & (col % 2 != 0)):
            #    subsystem = (-1) * subsystem
            #if (subsystem_Ly % 2 != 0) & (subsystem_Lx % 2 == 0) & (row % 2 != 0):
            #    # case uneven y, flip uneven rows
            #    subsystem = (-1) * subsystem
            #elif (subsystem_Lx % 2 != 0) & (subsystem_Ly % 2 == 0) & (col % 2 != 0):
            #    # case uneven x, flip uneven cols
            #    subsystem = (-1) * subsystem
            #elif (subsystem_Lx % 2 != 0) & (subsystem_Ly % 2 != 0) & ((col + row) % 2 != 0):
            #    # case both uneven, flip uneven cols
            #    print("flipping systems with col + row uneven ", row, " ", col)
            #    subsystem = (-1) * subsystem
            if (subsystem_Lx % 2 != 0):
                if ((subsystems_per_row % 2 != 0) & ((row + col) % 2 != 0)):
                    subsystem = (-1) * subsystem
                elif ((subsystems_per_row % 2 == 0) & ((col) % 2 != 0)):
                    print("flipping...")
                    subsystem = (-1) * subsystem

            x[row * subsystem_Ly : (row+1) * subsystem_Ly, col * subsystem_Lx : (col+1) * subsystem_Lx] = subsystem
    return x

from decimal import Decimal, ROUND_HALF_UP

def round_to_first_non_zero(number):
    dec_number = Decimal(str(number))

    # Find the position of the first non-zero digit
    nonzero_position = None
    for i, digit in enumerate(str(dec_number)):
        if digit != '0' and digit != '.':
            nonzero_position = i
            break
    decimal_position = None
    for i, digit in enumerate(str(number)):
        if digit == ".":
            decimal_position = i
    if nonzero_position is not None:
        # Round to the next digit after the first non-zero digit
        # decimal position is 1 if we have a float 0.something or 1.something
        # decimal position is greater than one if we have 10.something
        if(nonzero_position - decimal_position < 0):
            # means we have something in front of the comma
            rounded = round(number / (10 ** decimal_position), 1) * (10 ** decimal_position)
        else:
            rounded = round(number, nonzero_position - decimal_position)

        return rounded

    # If all digits are zeros, return the original number
    return number

def create_config(given_config, default_config):
    config = {}
    if given_config == None:
        return default_config
    else:
        for key in default_config.keys():
            try:
                config[key] = given_config[key]
            except KeyError:
                config[key] = default_config[key]
        # fontsizes are dealt with a bit different
        for font in default_config["fontsizes"]:
            config[font] = int(config[font] * (1 + config["increasefontsize"]))
        return config

PLOT_DEFAULT_CONFIG = {
    "ylabelsize": 11,
    "xlabelsize": 11,
    "titlesize": 14,
    "xtickfontsize": 9,
    "ytickfontsize": 9,
    "legendfontsize": 11,
    "increasefontsize": 0.0,
    "legendlocation": "best",
    "legendtitle": "",
    "ticklength": 6,
    "tickwidth": 2,
    "x_ticklength": 6,
    "x_tickwidth": 2,
    "y_ticklength": 6,
    "y_tickwidth": 2,
    "nr_y_minor_ticks": 5,
    "nr_y_major_ticks": 5,
    "nr_x_minor_ticks": 5,
    "nr_x_major_ticks": 5,
    "gridalpha": 0.5,
    "labelrotation": 0,
    "labelhorizontalalignment": "center",
    "labelverticalalignment": "center",
    "grid": True,
    "tight_layout": True,
    "legend": True,
    "spansfromdata" : False,
    "fontsizes": ["ylabelsize", "xlabelsize", "titlesize", "xtickfontsize", "ytickfontsize", "legendfontsize"]
}

def delete_last_line(filename):
    # Read all lines from the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Remove the last line
    if len(lines) > 0:
        lines.pop()

    # Write the remaining lines  back to the file
    with open(filename, 'w') as file:
        file.writelines(lines)

def get_10_power(x):
    if x == 0:
        return 0
    else:
        return int(np.floor(np.log10(abs(x))))

def configure_ax(fig, ax, config=None):
    """
    Takes a fig and an axes and configures the axes and stuff. If no config map is provided standard config is used
    :param fig:
    :param ax:
    :param config:
    :return: void
    """

    config = create_config(config, PLOT_DEFAULT_CONFIG)

    if config["spansfromdata"]:
        x_span, y_span, (xmin, xmax, ymin, ymax) = get_spans(ax) # the spans are at the moment dependent on the actuayl plotted values, but why?
    else:
        x_span = round_to_first_non_zero(ax.get_xlim()[1] - ax.get_xlim()[0])
        y_span = round_to_first_non_zero(ax.get_ylim()[1] - ax.get_ylim()[0])
    nr_y_minor_ticks = config["nr_y_minor_ticks"]
    nr_y_major_ticks = config["nr_y_major_ticks"]
    nr_x_minor_ticks = config["nr_x_minor_ticks"]
    nr_x_major_ticks = config["nr_x_major_ticks"]
    if ax.get_yscale() != "log":
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=y_span/nr_y_major_ticks))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=y_span/nr_y_major_ticks / nr_y_minor_ticks))
    if ax.get_xscale() != "log":
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_span/nr_x_major_ticks))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_span/nr_x_major_ticks / nr_x_minor_ticks))

    # We want to have inline ticks
    ax.tick_params(axis="x", direction='in', which='both', length=config["x_ticklength"], width=config["x_tickwidth"], labelsize=config["xtickfontsize"])
    ax.tick_params(axis="y", direction='in', which='both',
                   length=config["y_ticklength"], width=config["y_tickwidth"],
                   labelsize=config["ytickfontsize"])
    ax.tick_params(axis="x", direction='in', which='minor', length=int(config["x_ticklength"] * 0.75),
                   width=int(config["x_tickwidth"] * 0.75))
    ax.tick_params(axis="y", direction='in', which='minor', length=int(config["y_ticklength"] * 0.75),
                   width=int(config["y_tickwidth"] * 0.75))

    print("Y span = ", y_span)
    if y_span < 1e-3:
        print("Setting scientific mode?")
        ax.ticklabel_format(axis="y", style="sci")
    if x_span < 1e-3:
        ax.ticklabel_format(axis="x", style="sci")

    if ax.get_xscale() != "log":
        try:
            remove_origin_ticks(ax)
        except IndexError:
            pass

    # Füge Gitterlinien hinzu
    if config["grid"]:
        ax.grid(which='major', linestyle='--', alpha=config["gridalpha"])

    # title, achsenbeschriftungen, legend
    get_functions = [ax.get_xlabel, ax.get_ylabel]        # haha functions in a list, just python things
    set_functions = [ax.set_xlabel, ax.set_ylabel]
    default_texts = ["x", "y"]
    for i, get in enumerate(get_functions):
        if get() == "":
            # If we have empty string as title or stuff
            set_functions[i](default_texts[i])

    # rotate the y label
    ax.set_ylabel(ax.get_ylabel(), rotation=config["labelrotation"],
                  fontsize=config["ylabelsize"], ha=config["labelhorizontalalignment"],
                  va=config["labelverticalalignment"])
    ax.set_xlabel(ax.get_xlabel(), fontsize=config["xlabelsize"])
    ax.set_title(ax.get_title(), fontsize=config["titlesize"])
    #legend
    if config["legend"]:
        ax.legend(title=config["legendtitle"], fontsize=config["legendfontsize"], loc=config["legendlocation"],
                  title_fontsize=config["legendfontsize"])
    if config["tight_layout"]:
        plt.tight_layout()

def det_intersection(x, y_dic):
    min_total_dist = np.max(list(y_dic.values())[0])        # starting value for the min
    x_inter = 0
    y_inter = 0
    i_inter = 0
    for (i, x_val) in enumerate(x):
        total_dist = 0
        for key1 in y_dic.keys():
            for key2 in y_dic.keys():
                total_dist += np.abs(y_dic[key1][i] - y_dic[key2][i])
        total_dist /= len(y_dic.keys()) ** 2
        if (total_dist < min_total_dist) & (i > 1):
            min_total_dist = total_dist
            x_inter = x_val
            y_inter = list(y_dic.values())[0][i]
            i_inter = i
    return x_inter, y_inter, i_inter

def calc_diff_at(T_inter, T, value_dic):
    size_arr = []

    diff_beta_arr = []
    num_diff_arr = []
    for size in value_dic.keys():
        # okay we now think of a better method for numerical differentiation
        # first we look for the data point that is the closest to the intersection
        nearest_T_for_size, index_nearest_T = \
            find_nearest_value_and_index(T, T_inter)
        # now we calculate a central difference with the nearest value
        # being in the center
        # check whether the index is at least greater than zero...
        if index_nearest_T == 0:
            # vorwärtsdifferenz?
            num_diff = (value_dic[size][index_nearest_T + 1] - value_dic[size][index_nearest_T]) / (T[index_nearest_T + 1] - T[index_nearest_T])

        elif index_nearest_T == (len(T) - 1):
            # rückwärtsdiff
            num_diff = (value_dic[size][index_nearest_T] - value_dic[size][index_nearest_T - 1]) / (T[index_nearest_T] - T[index_nearest_T - 1])
        else:
            num_diff = (value_dic[size][index_nearest_T + 1] -
                        value_dic[size][index_nearest_T - 1]) \
                       / (2 * (T[index_nearest_T + 1] - T[index_nearest_T]))
        num_diff_arr.append(num_diff)
        size_arr.append(int(size))
    num_diff_arr = np.array(num_diff_arr)[np.argsort(size_arr)]
    size_arr = np.sort(size_arr)
    return num_diff_arr, size_arr

def simple_diff(x_eval, x_range, y_range):
    nearest_x, index = \
        find_nearest_value_and_index( x_range, x_eval)

    if nearest_x > x_eval:
        # Rückwärtsdiff
        num_diff = (y_range[index] - y_range[index - 1]) / (x_range[index] - x_range[index-1])
    elif nearest_x < x_eval:
        # Vorwärtsdiff
        num_diff = (y_range[index + 1] - y_range[index]) / (x_range[index + 1] - x_range[index])
    else:
        # Central diff
        num_diff = (y_range[index + 1] - y_range[index - 1]) / (x_range[index + 1] - x_range[index - 1])

    return num_diff

def interpolate_and_minimize(data_sets, res=100):
    """
    Determine the intersection of different sets of discrete x-y values.
    """
    # Extract x and y values from each data set
    x_values = [data_sets[size][list(data_sets[size].keys())[0]] for size in sorted(list(data_sets.keys()))]
    y_values = [data_sets[size][list(data_sets[size].keys())[1]] for size in sorted(list(data_sets.keys()))]
    # Find common x range TODO thats totally wrong right?
    common_x = np.unique(np.concatenate(x_values))
    x_range = np.linspace(common_x[0], common_x[-1], res)
    # Interpolate y values for each data set on the common x range

    interpolated_y_values = [interp1d(x, y, kind='linear', fill_value='extrapolate')(x_range) for x,y in zip(x_values, y_values)]
    dif = []
    y_mean = []

    for i in range(len(x_range)):
        # we just f-in calculate the error for every x
        y_vals = [y_set[i] for y_set in interpolated_y_values]
        diff = 0
        for j in range(len(y_vals)):
            for i in range(len(y_vals)):
                diff += abs(y_vals[i] - y_vals[j]) ** 2
        dif.append(diff)
        y_mean.append(np.mean(y_vals))
    min_ind = np.argmin(dif)
    y_intersect = y_mean[min_ind]
    x_intersect = x_range[min_ind]


    # Return the common x values and intersection y values
    return x_range, y_intersect, x_intersect, interpolated_y_values

def mark_point(ax, x, y, c="C0", label=""):
    ax.scatter(x, y, c=c)
    y_low = ax.get_ylim()[0]
    x_low = ax.get_xlim()[0]
    print("mark point called")
    ax.plot([x, x], [y_low, y], c=c, ls='--', lw=2, alpha=0.75, label=label)
    ax.plot([x_low, x], [y, y], c=c, ls='--', lw=2, alpha=0.75)

def remove_origin_ticks(ax):
    # Remove ticks in the origin of the plot
    get_functions = [ax.get_xticks, ax.get_yticks]
    set_functions = [ax.set_xticks, ax.set_yticks]
    origins = [ax.get_xlim()[0], ax.get_ylim()[0]]
    ends = [ax.get_xlim()[1], ax.get_ylim()[1]]
    for getter, setter, origin, end in zip(get_functions, set_functions, origins, ends):
        major_ticks = getter()
        minor_ticks = getter(minor=True)
        half_tick = (minor_ticks[1] - minor_ticks[0]) / 2
        new_minor_ticks = set(minor_ticks)
        new_major_ticks = set(major_ticks)
        # remove the tick if it is within half a minor tick of the origin value
        for tick in major_ticks:
            if (np.abs(origin - tick) < half_tick) | (tick < origin) | (tick > end):
                new_major_ticks.remove(tick)
        for tick in minor_ticks:
            if (np.abs(origin - tick) < half_tick) | (tick < origin) | (tick > end):
                new_minor_ticks.remove(tick)
        # set the new ticks
        setter(list(new_major_ticks))
        setter(list(new_minor_ticks), minor=True)

def get_spans(ax):
    # The base of those ticks should be read of the data
    xmin = np.inf
    xmax = -np.inf
    ymin = np.inf
    ymax = -np.inf
    for line in ax.get_lines():
        try:
            xmin = np.minimum(xmin, np.min(line.get_xdata()))
            xmax = np.maximum(xmax, np.max(line.get_xdata()))
            ymin = np.minimum(ymin, np.min(line.get_ydata()))
            ymax = np.maximum(ymax, np.max(line.get_ydata()))
        except ValueError:
            pass
    # check for already set custom limits or do we want to just set them with config?
    # TODO We probably should just have looked for the limits anyway?
    set_xmin = ax.get_xlim()[0]
    set_xmax = ax.get_xlim()[1]
    set_ymin = ax.get_ylim()[0]
    set_ymax = ax.get_ylim()[1]
    xmin = np.maximum(xmin, set_xmin)
    xmax = np.minimum(xmax, set_xmax)
    ymin = np.maximum(ymin, set_ymin)
    ymax = np.minimum(ymax, set_ymax)
    x_span = round_to_first_non_zero(xmax - xmin)
    y_span = round_to_first_non_zero(ymax - ymin)
    return x_span, y_span, (xmin, xmax, ymin, ymax)

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

def linear_corr(x, m, a, b):
    return a + m * x - b * x

def crit_poly_fit_corr(L, nu, A, B, omega):
    return A * L ** (1 / nu) * (1 + B * L **(-omega))

def ising_corr_poly_fit(L, A, B, omega):
    return A * L * (1 + B * L ** (-omega))

def crit_poly_fit(L, nu, A):
    return A * L ** (1/nu)

def find_nearest_value_and_index(x_arr, x_star):
    nearest_value = x_arr[0]
    nearest_index = 0
    min_diff = abs(x_star - x_arr[0])

    for i, x in enumerate(x_arr):
        diff = abs(x_star - x)
        if diff < min_diff:
            nearest_value = x
            nearest_index = i
            min_diff = diff

    return nearest_value, nearest_index

def get_intersection_index(y, z, x_y = [], x_z = []):
    def get_intersection_index(y, z):
        # assume two arrays y and z OF SAME SIZE and we find the index where the difference between them is the smallest
        dif = np.inf
        ind = 0
        for i in range(len(y)):
            difference = (np.abs(y[i] - z[i]))
            if difference < dif:
                dif = difference
                ind = i
        return ind
    if len(x_y) == 0:
        return get_intersection_index(y, z)
    else:
        dif = np.inf
        ind_i, ind_j = 0, 0
        for i,j in product(range(len(y)), range(len(z))):
            # we need a metric for the difference, lets take RS
            difference = np.sqrt((y[i]- z[j]) ** 2 + (x_y[i] - x_z[j]) ** 2)
            if difference < dif:
                dif =difference
                ind_i = i
                ind_j = j
        return (ind_i, ind_j)

def find_common_range(x1, x2, y1, y2):
    common_x_min = max(np.min(x1), np.min(x2))
    common_x_max = min(np.max(x1), np.max(x2))

    common_indices1 = (x1 >= common_x_min) & (x1 <= common_x_max)
    common_indices2 = (x2 >= common_x_min) & (x2 <= common_x_max)

    return x1[common_indices1], y1[common_indices1], y2[common_indices2]

def find_first_intersection(x1, x2, y1, y2, res=10000):
    # Interpolate the curves
    #print("diff_func:")
    x_range, _, _ = find_common_range(x1, x2, y1, y2)
    x_inter = np.linspace(x_range[0], x_range[-1], res)
    y1 = np.interp(x_inter, x1, y1)
    y2 = np.interp(x_inter, x2, y2)
    #print(x_range)
    #print(y1 - y2)
    diff_func = y1 - y2
    U_L_1_func = y1        # we need that for the U_L value at the intersection...
    #fig, ax = plt.subplots(1, 1)
    #ax.plot(x_inter, diff_func)
    #plt.show()
    #print("x_inter:", x_inter)
    #print("diff_func", diff_func)
    #print("np.argmin(diff_func(x_inter))", np.argmin(np.abs(diff_func)))
    #print("x_inter[np.argmin(diff_func(x_inter))]", x_inter[np.argmin(np.abs(diff_func))])
    #intersection_x = x_inter[np.argmin(np.abs(diff_func))]
    #print("diff_func(intersection_T)", diff_func(intersection_x))
    # Find the root (intersection) using numerical methods
    #print("diff_func.roots():", diff_func.roots())
    #intersection_x = diff_func.roots()[0]       # TODO does it work like that?

    # Okay I think we have to write a slightly better version of this function
    # We keep the diff_func, but we now check for sign changes
    diff_func_sign = np.sign(diff_func)
    signchange = ((np.roll(diff_func_sign, 1) - diff_func_sign) != 0).astype(int)
    signchange[0] = 0
    # we take the index of the first index change. This assumes that the lines dont cross for low temperatures
    # for high temperatures they might cross more often
    # This is still not very good?
    # We could do this sign thing for every diff curve and then choose always the
    # one that has the least total difference to the other sign arrays
    ind_signchange = np.where(signchange == 1)[0][0]
    intersection_x = x_inter[ind_signchange]
    U_L_intersection = U_L_1_func[ind_signchange]

    return intersection_x, U_L_intersection

def find_intersections(x_range, y1, y2, res=10000):
    # Interpolate the curves
    x_inter = np.linspace(x_range[0], x_range[-1], res)[1:]
    diff_func = np.interp(x_inter, x_range, y1 - y2)[1:]
    y_1_func = np.interp(x_inter, x_range, y1)[1:]        # we need that for the U_L value at the intersection...

    diff_func_sign = np.sign(diff_func)
    signchange = ((np.roll(diff_func_sign, 1) - diff_func_sign) != 0).astype(int)
    signchange[0] = 0

    ind_signchange = np.where(signchange == 1)[0]      # we exclude the first value because this means two U_L values are exactly the same and this means
    # either insane chance (I dont believe this could ever happen, or U_L = 1 if we are deep in the low temp region)
    intersections_x = []
    intersections_y = []
    for ind_sign in ind_signchange:
        intersections_x.append(x_inter[ind_sign])
        intersections_y.append(y_1_func[ind_sign])

    return intersections_x, intersections_y

def cut_data_around_peak(x_values, y_values, threshold_fraction=0.5, min_points_fraction=0.2):
    """
    reduces peaked data to the data around the peak
    :param threshold_fraction:  values below threshold_fraction * peak hight are cut
    :param min_points_fraction: we will be left at least with min_points_fraction of the datapoints we begun with
    :return: cut data
    """
    # we remove the offset of the y_values as they are not relevant for the peak position
    y_offset = np.min(y_values)
    y_values -= y_offset
    # Find the index and value of the peak
    peak_index = np.argmax(y_values)
    peak_value = y_values[peak_index]

    # Calculate the threshold based on the fraction of the peak value
    threshold = threshold_fraction * peak_value

    # Find the range of indices where the y-values are above the threshold
    above_threshold_indices = np.where(y_values > threshold)[0]
    # If the length of this is to small, we reduce the threshold fraction? Quick fix
    if len(above_threshold_indices) < min_points_fraction * len(y_values):
        # We always cut two values so if we only cut the two lowest values we cannot reach more values
        if len(above_threshold_indices) < len(x_values) - 2:
            new_threshold = 0.9 * threshold_fraction        # small steps?
            # We need to add back the offset so that it wont be zero next iteration
            y_values += y_offset
            return cut_data_around_peak(x_values, y_values, new_threshold, min_points_fraction)
    # I think if we
    # Extract the subset of data around the peak
    x_cut = np.array(x_values)[above_threshold_indices]
    y_cut = np.array(y_values)[above_threshold_indices]

    # we have to add the offset back again
    y_cut += y_offset

    return x_cut, y_cut

def cut_data_around_peak_order(x_values, y_values, threshold_fraction=0.5, min_points_fraction=0.2):
    """
    reduces peaked data to the data around the peak. Assumes the order that p = 0 is the first and then p > 0 and then p < 0
    :param threshold_fraction:  values below threshold_fraction * peak hight are cut
    :param min_points_fraction: we will be left at least with min_points_fraction of the datapoints we begun with
    :return: cut data
    """
    # we remove the offset of the y_values as they are not relevant for the peak position
    y_offset = np.min(y_values)
    y_values -= y_offset
    # Find the index and value of the peak
    peak_index = 0
    peak_value = y_values[peak_index]

    # Calculate the threshold based on the fraction of the peak value
    threshold = threshold_fraction * peak_value

    # Find the range of indices where the y-values are above the threshold
    cut_index = len(y_values)
    for ind, val in enumerate(y_values):
        if val < threshold:
            # means we found the index after which we want to cut.
            cut_index = ind
            break
    y_cut = np.concatenate((y_values[:cut_index], y_values[-cut_index:]))
    # If the length of this is to small, we reduce the threshold fraction? Quick fix
    if len(y_cut) < min_points_fraction * len(y_values):
        # We always cut two values so if we only cut the two lowest values we cannot reach more values
        if len(y_cut) < len(x_values) - 2:
            new_threshold = 0.9 * threshold_fraction        # small steps?
            # We need to add back the offset so that it wont be zero next iteration
            y_values += y_offset
            return cut_data_around_peak_order(x_values, y_values, new_threshold, min_points_fraction)
    # I think if we
    # Extract the subset of data around the peak
    x_cut = np.concatenate((x_values[:cut_index], x_values[-cut_index:]))

    # we have to add the offset back again
    y_cut += y_offset

    return x_cut, y_cut

def get_first_intersections(size_T_cum_dic, value_name):
    """
    returns the first intersection of every line pair
    :param size_T_cum_dic:
    :return:
    """
    intersections = []
    intersections_y = []
    sizes = list(size_T_cum_dic.keys())

    for (size1, size2) in combinations(sizes, 2):
        U_L_1 = size_T_cum_dic[size1][value_name]
        U_L_2 = size_T_cum_dic[size2][value_name]
        # We assume that the Temperatures are the same for the two curves...
        T_arr_1 = size_T_cum_dic[size1]["T"]
        T_arr_2 = size_T_cum_dic[size2]["T"]
        x_inter, y_inter = find_first_intersection(T_arr_1, T_arr_2, U_L_1, U_L_2)
        intersections.append(x_inter)
        intersections_y.append(y_inter)

    return intersections, intersections_y

def extract_cell(row_data, cell_nr, cell_L):
    L = int(np.sqrt(row_data.size))
    cells_per_row = L / cell_L
    col = cell_nr % cells_per_row       # number of the cell in its row
    row = cell_nr // cells_per_row

    cell = np.zeros(cell_L * cell_L)    # array for the values of the cell
    for j in range(cell_L):
        for i in range(cell_L):
            ind = int(L * (row * cell_L + j) + i + col * cell_L)
            cell[j * cell_L + i] = row_data[ind % (L * L)]

    return cell

def extract_cell_from_rectangle(row_data, cell_nr, cell_L, Lx, Ly):
    cells_per_row = Lx / cell_L
    col = cell_nr % cells_per_row       # number of the cell in its row
    row = cell_nr // cells_per_row

    cell = np.zeros(cell_L * cell_L)    # array for the values of the cell
    for j in range(cell_L):
        for i in range(cell_L):
            ind = int(
                # determines the left end of the lattice for the j-th row of the new cell
                Lx * (row * cell_L + j) +
                # iterates over everything to the right (left end of new cell to right end)
                i
                # displacement of the cell from the left end of the lattice
                + col * cell_L
            )
            cell[j * cell_L + i] = row_data[ind % int(Lx * Ly)]

    return cell

def extract_rectangle_from_rectangle(row_data, cell_nr, nr_cells, cell_Lx, cell_Ly, Lx, Ly):
    Lx = int(Lx)
    cell_Lx = int(cell_Lx)
    rows = np.zeros((int(cell_Ly), cell_Lx * nr_cells))

    # i can extract the rows of our combined subcells
    for row in range(int(cell_Ly)):
        rows[row] = row_data[row * Lx:row * Lx + cell_Lx * nr_cells]

    # we now somehow want to reshape this...
    cells_per_row = int(np.sqrt(nr_cells))
    row_length = cells_per_row * cell_Lx
    system = []

    for i in range(cells_per_row):
        for j in range(int(cell_Ly)):
            system.append(rows[j][i * row_length: (i + 1) * row_length])
    return np.array(system)

def extract_subsystem_from_matrix(matrix, Lx, cell_nr):
    """
    takes the matrix that is returned by reshape line and extracts the subsystem with cell_nr
    :param matrix:
    :param Lx:
    :param cell_nr:
    :return:
    """
    subsystem = matrix[:, cell_nr * Lx: (cell_nr + 1) * Lx]

    return subsystem

def reshape_line(line, dim_size_x, dim_size_y):
    """
    takes the line that is in your csv files and reshapes it to a 2D array with dimensions dim_size_y, dim_size_x
    :param line:
    :param dim_size_x:
    :param dim_size_y:
    :return:
    """

    line = np.array(line)
    mat = line.reshape((dim_size_y, dim_size_x))
    return mat

def list_directory_names(path):
    paths = os.listdir(path)
    dirs = []
    for dir in paths:
        if os.path.isdir(os.path.join(path, dir)):
            dirs.append(dir)
    return dirs

def second_order_num_diff(x, y):
    """
    takes two arrays and returns array with numerical diff dx/dy with second order accuracy if possible
    :param x: x-values
    :param y: depenant values
    :return: numerical diff array
    """
    dif = np.zeros(len(x))
    for i in range(len(x)):
        if (i >= 2) and (i < len(x) - 2):
            h = x[i+1] - x[i]
            # second order diff
            dy_dx = (- y[i + 2] + 8 * y[i + 1] - 8 * y[i - 1] + y[i - 2]) / (12 * h)
        elif i == 1 or i == len(x) - 2:
            # first oder central diff
            h = x[i+1] - x[i]
            dy_dx = 1 / (2 * h) * (y[i+1] - y[i-1])
        elif i == 0:
            h = x[i+1] - x[i]
            # vorwärtsdifferenz
            dy_dx = (y[i + 1] - y[i]) / h
        else:
            h = x[i] - x[i - 1]
            # rückwärtsdiffernz
            dy_dx = (y[i] - y[i - 1]) / h
        dif[i] = dy_dx
    return dif

def find_fwhm(x, y):
    """
    Find the full width at half maximum (FWHM) of a symmetric peak.

    Parameters:
    - x: Array of x-values
    - y: Array of y-values

    Returns:
    - fwhm: Full width at half maximum
    """
    # Find the index of the maximum y-value
    max_index = np.argmax(y)
    min_index = np.argmin(y)
    max_value = y[max_index]
    min_value = y[min_index]
    max_value = max_value - min_value
    # Define an interpolation function
    interpolate = interp1d(x, y - max_value / 2, kind='cubic')
    print(interpolate)
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, interpolate(x))
    ax.plot(x, y - max_value/2, ls="", marker="x", ms="1")

    # Find the roots (x-values) of the interpolation
    # Evaluate the interpolation at a set of points
    sample_points = np.linspace(min(x)/(1.5), max(x)/(1.5), 2000000)
    sample_values = interpolate(sample_points)

    # Find the indices where the sign changes (roots)
    sign_changes = np.where(np.diff(np.sign(sample_values)))[0]

    # Find the two roots closest to the maximum x-value
    closest_roots = sample_points[sign_changes]
    closest_roots = closest_roots[np.argsort(np.abs(closest_roots - x[max_index]))[:2]]

    ax.plot(closest_roots, interpolate(closest_roots), marker="o")
    plt.show()
    # Calculate FWHM
    fwhm = np.abs(closest_roots[1] - closest_roots[0])

    return fwhm

def critical_amplitude(eps, xi0):
    return xi0 / (eps ** 1)

def process_file_old(file_path, threshold, key='t', value='U_L', L=np.inf, cap=0.2):
    """
    Process a single file and calculate the average after the given threshold.
    Will now also return the error on the average
    """
    df = pd.read_csv(file_path)
    if 0 < threshold < 1:
        threshold = threshold * len(df[key])
    df = df[int(threshold):]
    nr_values = df.shape[0]
    if nr_values < 1:
        average_value = 0
        error = 0
    else:
        xi_values = np.array(df[value])
        # We set the values that are larger than 20% of the system size to be 20% of the system size
        xi_values = np.minimum(xi_values, L * cap)  # TODO okay haha you didnt think at all about that this is also supposed to work for U_L but it makes no difference
        average_value = np.mean(xi_values)
        stddev = np.std(xi_values)
        error = stddev / np.sqrt(nr_values)
    return average_value, error, nr_values

def process_file(file_path, threshold, key='t', value='U_L'):
    print("process file", file_path)
    file_path = Path(file_path)
    try:
        para_path = str(file_path.with_suffix(".txt"))
        parameters = read_parameters_txt(para_path)
        #print(parameters)
        f_avg = parameters[value]
        rel_error = parameters[f"{value}_error"]
        moving_factor = parameters["moving_factor"]
        total_nr_values = 0
    except:
        para_path = str(file_path.with_suffix(".txt"))
        print(f"couldnt read {para_path}")
        df = read_large_df(file_path, skiprows=1, sep="", cut_endline=1)
        if 0 < threshold < 1:
            threshold = threshold * len(df)
        total_nr_values = len(df)
        df = df[int(threshold):]
        nr_values = len(df)
        df = np.array(df)
        if value == "xiy":
            f = np.array(df[:, 2])
        else:
            f = np.array(df[:, 1])
        times = np.array(df[:, 0])


        f_avg = np.mean(f)
        f_dist_var = np.maximum(np.var(f), 1e-7)        # TODO quickfix because error of zero is unrealistic

        ds = np.maximum(times[1] - times[0], 1e-7)
        # try:
        #     # TODO this depends on the value, we probably should save this with the name, we currently not do this for U_L
        #     autocorr_time = paras["autocorrelation_time_" + value]
        # except KeyError:
        autocorr_time = integrated_autocorr_time(f, ds)
        print("autocorr_time ", autocorr_time, "  ds ", ds, " f_dist_var", f_dist_var)
        variance = autocorr_time / (nr_values * ds) * f_dist_var
        error = np.sqrt(variance)

        rel_error = error / f_avg
        # What are we doing with the moving factor? I do not really want
        # to return it here but if we already read the file...
        # I guess just return it, it is one of the useful values you want to extract
        # from a file

        # Moving factor
        #try:
        #    moving_factor = paras["moving_factor_" + value]
        #except KeyError:
        #    # We already have a function that calculates it?
        moving_factor = getMovingFactor(f, f_avg)

        parameters[value] = f_avg
        parameters[f"{value}_error"] = rel_error
        parameters["moving_factor"] = moving_factor

        # we would aslo have to write them to the same txt file?
        write_dict_to_file(parameters, para_path)

    return f_avg, rel_error, moving_factor, total_nr_values

def process_mag_file_to_U_L(file_path, threshold, key='t', value='m'):
    print(file_path)
    file_path = Path(file_path)
    try:
        parameters = read_parameters_txt(str(file_path.with_suffix(".txt")))
        U_L = parameters["U_L"]
        rel_error = parameters["U_L_error"]
        moving_factor = parameters["moving_factor"]
        total_nr_values = 0
    except:
        df = pd.read_csv(file_path, sep=";")
        print("fing couldnt read this stuff")
        if 0 < threshold < 1:
            threshold = threshold * len(df[key])
        total_nr_values = df.shape[0]
        df = df[int(threshold):]
        nr_values = df.shape[0]
        m = np.array(df[value])
        times = np.array(df[key])

        m_avg = np.mean(m)
        m2 = np.mean(m ** 2)
        m2_err = np.std(m ** 2) / np.sqrt(len(m))     # std is standard dev of dist
        m4 = np.mean(m ** 4)
        m4_err = np.std(m ** 4) / np.sqrt(len(m))

        U_L = m4 / m2 ** 2

        U_L_error = np.sqrt(pow(1 / m2 / m2 * m4_err, 2) + pow(2 * m4 / pow(m2, 3) * m2_err, 2))

        #f_dist_var = np.maximum(np.var(f), 1e-7)        # TODO quickfix because error of zero is unrealistic
    #
        #ds = times[1] - times[0]
        ## try:
        ##     # TODO this depends on the value, we probably should save this with the name, we currently not do this for U_L
        ##     autocorr_time = paras["autocorrelation_time_" + value]
        ## except KeyError:
        #autocorr_time = integrated_autocorr_time(f, ds)
    #
        #variance = autocorr_time / (nr_values * ds) * f_dist_var
        #error = np.sqrt(variance)

        rel_error = U_L_error / U_L
        # What are we doing with the moving factor? I do not really want
        # to return it here but if we already read the file...
        # I guess just return it, it is one of the useful values you want to extract
        # from a file

        # Moving factor
        #try:
        #    moving_factor = paras["moving_factor_" + value]
        #except KeyError:
        #    # We already have a function that calculates it?
        moving_factor = getMovingFactor(m, m_avg)

    return U_L, rel_error, moving_factor, total_nr_values

def process_new_mag_file_to_U_L(file_path, threshold, key='t', value='m'):
    print("This?", file_path)
    file_path = Path(file_path)
    # df = pd.read_csv(file_path, sep=";")
    try:
        para_path = str(file_path.with_suffix(".txt"))
        parameters = read_parameters_txt(para_path)
        #print(parameters)
        U_L = parameters["U_L"]
        rel_error = parameters["U_L_error"]
        moving_factor = parameters["moving_factor"]
        total_nr_values = 0
    except:
        print(f"couldnt read {para_path}")
        df = read_large_df(file_path, skiprows=1, sep=";")
        m = []
        for ind, line in enumerate(df):
            m += line
        m = np.array(m)

        parameters["equil_error"] = threshold

        if 0 < threshold < 1:
            threshold = int(threshold * len(m))
        total_nr_values = len(m)
        print("threshold ", int(threshold))
        if threshold:
            m = m[threshold:]
        #times = np.array(df[key])

        m_avg = np.mean(m)
        m2 = np.mean(m ** 2)
        m2_err = np.std(m ** 2) / np.sqrt(len(m))     # std is standard dev of dist
        m4 = np.mean(m ** 4)
        m4_err = np.std(m ** 4) / np.sqrt(len(m))

        U_L = m4 / m2 ** 2

        U_L_error = np.sqrt(pow(1 / m2 / m2 * m4_err, 2) + pow(2 * m4 / pow(m2, 3) * m2_err, 2))

        #f_dist_var = np.maximum(np.var(f), 1e-7)        # TODO quickfix because error of zero is unrealistic
    #
        #ds = times[1] - times[0]
        ## try:
        ##     # TODO this depends on the value, we probably should save this with the name, we currently not do this for U_L
        ##     autocorr_time = paras["autocorrelation_time_" + value]
        ## except KeyError:
        #autocorr_time = integrated_autocorr_time(f, ds)
    #
        #variance = autocorr_time / (nr_values * ds) * f_dist_var
        #error = np.sqrt(variance)

        rel_error = U_L_error / U_L
        # What are we doing with the moving factor? I do not really want
        # to return it here but if we already read the file...
        # I guess just return it, it is one of the useful values you want to extract
        # from a file

        # Moving factor
        #try:
        #    moving_factor = paras["moving_factor_" + value]
        #except KeyError:
        #    # We already have a function that calculates it?
        moving_factor = getMovingFactor(m, m_avg)

        parameters["m"] = m_avg
        parameters["U_L"] = U_L
        parameters[f"U_L_error"] = rel_error
        parameters["moving_factor"] = moving_factor

        # we would aslo have to write them to the same txt file?
        write_dict_to_file(parameters, para_path)

    return U_L, rel_error, moving_factor, total_nr_values

def recalculate_mag_file_to_U_L(file_path, threshold, key='t', value='m'):
    df = read_large_df(file_path, skiprows=1, sep=";")

    nr_values = len(df)
    m = []
    for ind, line in enumerate(df):
        m += line
    m = np.array(m)
    if 0 < threshold < 1:
        threshold = int(threshold * len(m))
    total_nr_values = len(m)
    if threshold:
        m = m[threshold:]
    #times = np.array(df[key])

    m_avg = np.mean(m)
    #print("m = ", m_avg)
    m2 = np.mean(m ** 2)
    #print("m2 = ", m2)
    m2_err = np.std(m ** 2) / np.sqrt(len(m))     # std is standard dev of dist
    m4 = np.mean(m ** 4)
    #print("m4 = ", m4)
    m4_err = np.std(m ** 4) / np.sqrt(len(m))

    U_L = m4 / m2 ** 2

    U_L_error = np.sqrt(pow(1 / m2 / m2 * m4_err, 2) + pow(2 * m4 / pow(m2, 3) * m2_err, 2))

    #f_dist_var = np.maximum(np.var(f), 1e-7)        # TODO quickfix because error of zero is unrealistic
#
    #ds = times[1] - times[0]
    ## try:
    ##     # TODO this depends on the value, we probably should save this with the name, we currently not do this for U_L
    ##     autocorr_time = paras["autocorrelation_time_" + value]
    ## except KeyError:
    #autocorr_time = integrated_autocorr_time(f, ds)
#
    #variance = autocorr_time / (nr_values * ds) * f_dist_var
    #error = np.sqrt(variance)

    rel_error = U_L_error / U_L
    # What are we doing with the moving factor? I do not really want
    # to return it here but if we already read the file...
    # I guess just return it, it is one of the useful values you want to extract
    # from a file

    # Moving factor
    #try:
    #    moving_factor = paras["moving_factor_" + value]
    #except KeyError:
    #    # We already have a function that calculates it?
    moving_factor = getMovingFactor(m, m_avg)
    print(f"U_L = {U_L}, error = {rel_error}")

    return U_L, rel_error, moving_factor, total_nr_values

def recalculate_vectorial_mag_file_to_U_L(file_path, threshold, key='t', value='m'):
    df = read_large_df(file_path, skiprows=1, sep=";")

    m = []
    for ind, line in enumerate(df):
        m += line
    m = np.array(m)
    mx = m[::2]
    my = m[1::2]
    if 0 < threshold < 1:
        threshold = int(threshold * len(mx))
    total_nr_values = len(mx)
    if threshold:
        mx = mx[threshold:]
        my = my[threshold:]


    m_avg = np.mean(m)

    m2 = np.mean(mx ** 2 + my ** 2)

    m2_err = np.std(mx ** 2 + my ** 2) / np.sqrt(len(mx))     # std is standard dev of dist
    m4 = np.mean((mx ** 2 + my ** 2) ** 2)

    m4_err = np.std((mx ** 2 + my ** 2) ** 2) / np.sqrt(len(mx))

    U_L = m4 / m2 ** 2

    U_L_error = np.sqrt(pow(1 / m2 / m2 * m4_err, 2) + pow(2 * m4 / pow(m2, 3) * m2_err, 2))

    rel_error = U_L_error / U_L
    moving_factor = getMovingFactor(m, m_avg)
    print(f"U_L = {U_L}, error = {rel_error}")

    return U_L, rel_error, moving_factor, total_nr_values

def process_size_folder(size_folder, threshold, key='T', value='U_L', file_ending='cum',
                        selected_temperatures=None, process_file_func=process_file):
    """
    Process all files in a size folder and return a dictionary with temperature and average values.
    """
    eps = 1e-4
    result = {key: [], value: [], "error": []}
    for temp_folder in os.listdir(size_folder):
        if (temp_folder != "plots") & (temp_folder[0] != "."):
            if selected_temperatures is not None:
                temp = float(temp_folder)
                selected = False
                # we have to look if approximately this temperature is present in selected_temperatures
                for sel_temp in selected_temperatures:
                    if sel_temp - eps < temp < sel_temp + eps:
                        selected = True
                if not selected:
                    continue

            temp_folder_path = os.path.join(size_folder, temp_folder)
            if os.path.isdir(temp_folder_path):
                val_avg, error, _, _ = process_temp_folder(temp_folder_path, threshold, file_ending, value, process_file_func)
                if val_avg:
                    result[key].append(float(temp_folder))
                    result[value].append(val_avg)
                    result["error"].append(error)

    result[value] = np.array(result[value])[np.argsort(result[key])]
    result["error"] = np.array(result["error"])[np.argsort(result[key])]
    result[key] = np.sort(result[key])
    # TODO if i need the error i can return it here but not for now
    return result

def process_temp_folder_old(temp_folder_path, threshold, file_ending, value, process_file_function=process_file_old):
    # So this averages the .cum or .corr files
    # We might need to calculate an error here, but how do we get the total error? -> gaussian error propagation
    temp_average = []
    temp_error = []
    nr_avg_values = []
    nr_subsystems = []
    for file_name in os.listdir(temp_folder_path):
        if file_name.endswith(file_ending):
            # it could be that some files have a cumulant averaged out of more subsystems, that shoudl be taken into consideration
            file_path = os.path.join(temp_folder_path, file_name)
            para_file_path = os.path.splitext(file_path)[0] + ".txt"
            parameters = read_parameters_txt(para_file_path)
            nr_subsys = parameters["nr_subsystems"]
            # fing sh... I guess It doesnt matter but it is so ugly
            if value== "U_L":
                L = np.inf
            else:
                direction = value[-1]
                L = parameters["subsystem_L" + direction]         # We need the subsystem sizes to see if the values we average make sense
            average_value, error, nr_values = process_file_function(file_path, threshold, 't', value, L=L)       # cap stays 20% for now I guess
            temp_average.append(average_value)
            temp_error.append(error)
            nr_avg_values.append(nr_values)
            nr_subsystems.append(nr_subsys)
    if temp_average:  # Why is here an if this should not be here?? If we are in a temperature folder we want to geht values back, am I right? If we dont have anything at all, we made something wrong before?
        N = np.sum(np.array(nr_avg_values) * np.array(nr_subsystems))
        val_avg = np.sum(np.array(temp_average) * np.array(nr_avg_values) * np.array(nr_subsystems)) / N
        error = 1 / N * np.sqrt(np.sum(np.array(nr_avg_values) * np.array(nr_subsystems) * np.array(temp_error)))
    else:  # Okay I dont care but this should not be like this I think. Is this even on the right plane? Shouldnt it be one bacK? OMG I dont understand my own code!!!
        val_avg = None
    return val_avg, error

def process_temp_folder(temp_folder_path, threshold, file_ending, value, process_file_func=process_file):
    # So this averages the .cum or .corr files
    # We might need to calculate an error here, but how do we get the total error? -> gaussian error propagation
    file_averages = []
    moving_factors = []
    file_errors = []
    nr_values = []
    for file_name in os.listdir(temp_folder_path):
        if file_name.endswith(file_ending):
            # it could be that some files have a cumulant averaged out of more subsystems, that shoudl be taken into consideration
            file_path = os.path.join(temp_folder_path, file_name)
            para_file_path = os.path.splitext(file_path)[0] + ".txt"
            # fing sh... I guess It doesnt matter but it is so ugly
            average_value, error, moving_factor, nr_vals = process_file_func(file_path, threshold, 't', value)
            file_averages.append(average_value)
            file_errors.append(error)
            moving_factors.append(moving_factor)
            nr_values.append(nr_vals)
    if file_averages:  # Why is here an if this should not be here?? If we are in a temperature folder we want to geht values back, am I right? If we dont have anything at all, we made something wrong before?
        # What are we doing now for the average, If one run had more subsystems this should just be reflecting in the error? So we just do
        # a weighted averagae here?
        # Frequently used weights seem to be the inverse variance
        file_errors = np.array(file_errors)
        weights = 1 / (file_errors ** 2)
        val_avg = np.average(file_averages, weights=weights)
        error = np.sqrt(1 / np.sum(weights))
        # Do we want to return an average moving factor? A weighted moving factor? Just the list?
    else:  # Okay I dont care but this should not be like this I think. Is this even on the right plane? Shouldnt it be one bacK? OMG I dont understand my own code!!!
        val_avg = None
        error = None
    return val_avg, error, moving_factors, nr_values

def average_ft(folderpath, ending=".ft"):
    files = os.listdir(folderpath)
    t_ft_k = {}
    t_ft_l= {}
    first_file = True
    nr_files = 0
    for file in (files):
        if file.endswith(ending):
            filepath = os.path.join(folderpath, file)
            df = pd.read_csv(filepath, sep=";")
            #print(filepath, "firstfile = ", first_file)
            for j, t in enumerate(df['t']):
                if first_file:
                    t_ft_k[t] = string_to_array(df["ft_k"][j])
                    t_ft_l[t] = string_to_array(df["ft_l"][j])
                else:
                    t_ft_k[t] += string_to_array(df["ft_k"][j])
                    t_ft_l[t] += string_to_array(df["ft_l"][j])
            first_file = False
            nr_files += 1
    # average
    #print(f"averaged {nr_files} files")
    for t in t_ft_k:
        t_ft_k[t] /= nr_files
        t_ft_l[t] /= nr_files
    return t_ft_k, t_ft_l

def average_ft_unequal_times(folderpath, ending=".ft"):
    files = os.listdir(folderpath)
    t_ft_k = {}
    t_ft_l= {}
    nr_files = 0
    para_file = find_first_txt_file(folderpath)
    parameters = read_parameters_txt(para_file)
    T_start = parameters["starting_temp"]
    T_end = parameters["end_temp"]
    tau = parameters["tau"]
    try:
        end_equil_time = parameters["equil_time_end"]
    except KeyError:
        end_equil_time = 0
    try:
        gamma = parameters["gamma"]
    except KeyError:
        gamma = 1
    quench_time = (T_start - T_end) ** (1 / gamma) * tau
    for file in (files):
        if file.endswith(ending):
            filepath = os.path.join(folderpath, file)
            df = pd.read_csv(filepath, sep=";")
            times = df['t']
            dt = times[1]- times[0]
            end_time = np.max(times)
            start_equil_time = end_time - quench_time - end_equil_time
            for j, t in enumerate(times):
                shift_t = t - start_equil_time
                try:
                    t_ft_k[shift_t] += string_to_array(df["ft_k"][j])
                    t_ft_l[shift_t] += string_to_array(df["ft_l"][j])
                except KeyError:
                    t_ft_k[shift_t] = string_to_array(df["ft_k"][j])
                    t_ft_l[shift_t] = string_to_array(df["ft_l"][j])
            nr_files += 1
    # average
    print(f"averaged {nr_files} files")
    for t in t_ft_k:
        t_ft_k[t] /= nr_files
        t_ft_l[t] /= nr_files

    t_ft_k = combine_close_values(t_ft_k, 1 / 2 * dt)
    t_ft_l = combine_close_values(t_ft_l, 1 / 2 * dt)
    return t_ft_k, t_ft_l

def xi_div(eps, xi0, nu):
    return xi0 / (np.abs(eps) ** nu)


def combine_close_values(dictionary, threshold):
    sorted_keys = sorted(dictionary.keys())
    combined_dict = {}
    current_key = sorted_keys[0]
    current_value_sum = dictionary[current_key]
    count = 1

    for key in sorted_keys[1:]:
        if key - current_key <= threshold:
            current_value_sum += dictionary[key]
            count += 1
        else:
            average_value = current_value_sum / count
            combined_dict[current_key] = average_value
            current_key = key
            current_value_sum = dictionary[key]
            count = 1

    # Add the last key-value pair
    average_value = current_value_sum / count
    combined_dict[current_key] = average_value

    return combined_dict

def average_lastlines_ft(folderpath, ending=".ft", nr_add_lines=10):
    # like average_ft but it only returns the ft for the latest t
    files = os.listdir(folderpath)
    t_ft_k = np.array([])
    t_ft_l = np.array([])
    first_file = True
    nr_files = 0
    for file in (files):
        if file.endswith(ending):
            print(file)
            filepath = os.path.join(folderpath, file)
            df = pd.read_csv(filepath, sep=";")
            #print(filepath, "firstfile = ", first_file)
            if first_file:
                # print("df['ft_k'].iloc[-1]")
                # print(df["ft_k"].iloc[-1])
                t_ft_k = string_to_array(df["ft_k"].iloc[-1])
                t_ft_l = string_to_array(df["ft_l"].iloc[-1])
            else:
                t_ft_k += string_to_array(df["ft_k"].iloc[-1])
                t_ft_l += string_to_array(df["ft_l"].iloc[-1])
            first_file = False

            for i in np.arange(0, nr_add_lines):
                t_ft_k += string_to_array(df["ft_k"].iloc[-(1 + nr_add_lines - i)])
                t_ft_l += string_to_array(df["ft_l"].iloc[-(1 + nr_add_lines - i)])

            nr_files += 1
    # average
    #print(f"averaged {nr_files} files")
    t_ft_k /= nr_files * (nr_add_lines + 1)
    t_ft_l /= nr_files * (nr_add_lines + 1)
    return t_ft_k, t_ft_l

def get_avail_tau_dic(root_dir):
    directory_dict = {}
    sizes = find_size_folders(root_dir)
    for size_folder in sizes:
        size_path = os.path.join(root_dir, size_folder)
        if os.path.isdir(size_path) and size_folder[0] != ".":
            for tau_folder in os.listdir(size_path):
                tau_path = os.path.join(size_path, tau_folder)
                if os.path.isdir(tau_path) and tau_folder != "plots" and tau_folder[0] != ".":
                    if float(tau_folder) not in directory_dict:
                        directory_dict[float(tau_folder)] = []
                    directory_dict[float(tau_folder)].append(int(size_folder))
    return directory_dict

def get_frequencies_fftw_order(nr_times):
    K = nr_times // 2
    freqs = []
    for i in range(nr_times):
        if i + K < nr_times:
            ind = i + K
        else:
            if nr_times % 2 != 0:
                ind = i - K - 1
            else:
                ind = i - K
        freqs.append(2 * np.pi * (ind - K) / nr_times)
    return freqs

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))
def lorentz_ft(x, xi, a, b):
    return b + a * xi ** 2 / (1 + (x) ** 2 * xi ** 2)

def MF_lorentz(x, xi, a):
    return a * xi / (1 + (x) ** 2 * xi ** 2)

def lorentz_pi(x, xi, a):
    return a * 2 * xi / (1 + (4 * np.pi ** 2) * (xi ** 2) * (x ** 2))

def lorentz_offset(x, xi, a, b):
    return b + a * 2 * xi / (1 + (xi ** 2) * (x ** 2))
def fit_lorentz(p, ft, fitfunc=MF_lorentz, errors=None):
    try:
        popt, pcov = curve_fit(fitfunc, p, ft, sigma=errors, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))

    except RuntimeError:
        print("RuntimeError in your function fit_lorentz. Maybe no fit possible because"
              " of to many fluctuations or because it is just flat. Returning 0")
        return (0, 0), None
        exit()
        # function has form of a delta peak
        # we delete the largest value
        ind = np.argmax(ft)

        ft = np.insert(ft, ind + 1, 1/2 * ft[ind])
        ft = np.insert(ft, ind, 1 / 2 * ft[ind])
        p = np.insert(p, ind + 1, p[ind] + 1/2 * p[ind + 1])
        p = np.insert(p, ind, (p[ind] + 1 / 2 * p[ind-1]))

        #ft = np.delete(ft, ind)
        #p = np.delete(p, ind)
        print("had to insert values")
        print(p)

        # maybe not cut off the maximum but insert values?

        return fit_lorentz(p, ft)
    return popt, perr

def pi_formatter(value, ticknumber):
    fraction = value / np.pi
    fraction = fraction.as_integer_ratio()
    if fraction[0] == 0:
        return "0"
    elif fraction[0] < 0:
        return r"$-\frac{" + str(abs(fraction[0])) +  "}{" + str(fraction[1]) + "} \pi$"
    else:
        return r"$\frac{" + str(fraction[0]) +  "}{" + str(fraction[1]) + "} \pi$"

def pi_formatter_fraction(value, ticknumber):
    fraction = value / np.pi
    fraction = fraction.as_integer_ratio()
    zähler = abs(fraction[0])
    if zähler == 1:
        zähler = ""
    else:
        zähler = str(zähler)

    if fraction[0] == 0:
        return "0"
    elif fraction[0] < 0:
        return r"$-\frac{" + zähler +  "\,\pi\,}{" + str(fraction[1]) + "}$"
    else:
        return r"$\frac{" + zähler +  "\,\pi\,}{" + str(fraction[1]) + "}$"

def cut_zero_imp(p, ft, ft_err=None):
    # nr cuts is the number of values to cut around zero
    p, ft = np.array(p), np.array(ft)
    # be aware this only works if we have an even lenght and the zero is included...
    if len(p) % 2 == 0:
        ft = ft[p != 0]
        if ft_err:
            ft_err = np.array(ft_err)
            ft_err = ft_err[p != 0]
            p = p[p != 0]
            return p, ft, ft_err
        else:
            p = p[p != 0]
            return p, ft
    else:
        # TODO if not we will cut the first value. Quickfix?!?!
        ft = ft[1:]
        p = p[1:]
        return p, ft

def cut_around_zero(p, ft, nr_additional_cuts = 0):
    # assumes order
    if nr_additional_cuts == 0:
        p = p[1:]
        ft = ft[1:]
    else:
        p = p[1 + nr_additional_cuts // 2:-nr_additional_cuts // 2]
        ft = ft[1 + nr_additional_cuts // 2:-nr_additional_cuts // 2]
    return p, ft

def rescale_t(t, tau, t_eq, zoom = 1):
    total_time = np.max(t)
    t = np.array(t) - t_eq
    new_t = []
    # the values that sit between t_eq and total_time - t_eq shall be scaled
    new_t += list(t[t < 0])
    new_t += list(t[(t >= 0) & (
            t < total_time - 2 * t_eq)] / tau * zoom)
    new_t += list(t[t >= total_time - 2 * t_eq] - (
            1 - zoom / tau) * (
                          total_time - 2 * t_eq))
    t_q_s = (total_time - 2 * t_eq) / tau * zoom
    return np.array(new_t), t_q_s

def plot_struct_func(px, py, fx, fy, error_x=None, error_y=None):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axx = axes[0]
    axy = axes[1]

    axx.errorbar(px, fx, yerr=error_x, ls=" ", marker=".", label="Structure Func", color="C1", ecolor="black", capsize=3)
    axx.set_xlabel(r"$p_x$")
    axx.set_ylabel(r"$S(p_x)$")
    axy.errorbar(py, fy, yerr=error_y, ls=" ", marker=".", label="Structure Func", color="C1", ecolor="black", capsize=3)
    axy.set_xlabel(r"$p_y$")
    axy.set_ylabel(r"$S(p_y)$")
    axx.legend()
    axy.legend()
    return fig, axes

def derivative_F_BJ_Ising(BJ):
    num = 1 / np.tanh(2 * BJ) * (1 / np.sinh(2 * BJ))
    dom = np.sqrt((1 / np.sinh(2 * BJ)) ** 2 + 1)
    return num / dom
def T_c_Ising(T, J1, J2):
    return np.sinh(2 * np.abs(J1) / T) * np.sinh(2 * np.abs(J2) / T) - 1

def T_c_Ising_eff(BJ1, BJ2):
    return np.sinh(2 * np.abs(BJ1)) * np.sinh(2 * np.abs(BJ2)) - 1

def T_c_Ising_eff_ratio(BJ, f):
    return np.sinh(2 * np.abs(BJ * f)) * np.sinh(2 * np.abs(BJ)) - 1

def T_c_est_Ising_eff_ratio(f):
    BJ2 = fsolve(T_c_Ising_eff_ratio, 1 / f, args=(f), maxfev=10000)
    return BJ2

def T_c_eff_ratio(BJ, f):
    return np.sinh(2 * np.abs(BJ * f)) * np.sinh(2 * np.abs(BJ)) - 1

def T_c_est_Ising_eff_ratio(f):
    BJ2 = fsolve(T_c_Ising_eff_ratio, 1 / f, args=(f), maxfev=10000)
    return BJ2

def T_c_est_Ising_eff(BJ):
    BJ2 = fsolve(T_c_Ising_eff, 1 / BJ, args=(BJ), maxfev=10000)
    return BJ2

def T_c_XY(T, J_large, J_small):
    return 2 * T / J_large * np.log(2 * T / J_small) - 1

def T_c_est(J_para, J_perp, h=0):
    # we need to solve the transcendet equation of the XY model
    # scipy f_solve needs a guess of the critical temperature which is difficult but lets try with the mean of the Js
    T_c_est = (J_para + J_perp) / 2
    if J_para > J_perp:
        T_c_est = fsolve(T_c_XY, T_c_est, args=(J_para, J_perp))
    else:
        T_c_est = fsolve(T_c_XY, T_c_est, args=(J_perp, J_para))
    return T_c_est

def BJ_x(BJ_y):
    # calculates the critical effective coupling beta * Jx in dependence of the
    # given coupling beta * Jy
    BJ_x = (-2) * np.log(BJ_y / 2)
    if BJ_x < BJ_y:
        # if thats the case we used the wrong formula
        # is that all legal? it feels so illegal
        BJ_x = 2 * np.exp(-BJ_y / 2)
    return BJ_x

def BJ_assert(BJ_x, BJ_y):
    # has to return zero otherwise something is wrong
    if np.abs(BJ_x) > np.abs(BJ_y):
        return np.abs(BJ_x) / 2 + np.log(BJ_y / 2)
    else:
        return np.abs(BJ_y) / 2 + np.log(BJ_x / 2)

def BJ_large(BJ_small):
    return (-2) * np.log(np.abs(BJ_small) / 2)

def BJ_small(BJ_large):
    return 2 * np.exp(-np.abs(BJ_large) / 2)

def BJ_est(BJ_para, BJ_perp, h):
    BJ_est = (BJ_para + BJ_perp) / 2        # woah I dont have any clue on how to estimate this

def find_closest_key(dictionary, target_value):
    closest_key = None
    min_difference = float('inf')

    for key in dictionary:
        difference = abs(key - target_value)
        if difference < min_difference:
            min_difference = difference
            closest_key = key

    return closest_key

def p_approximation(p_guess, theta_equil, J_para, J_perp, h):
    """
    Approximates the p we have to use to achieve the desired theta_equil
    :param p_guess: A guess of the p we have to use
    :param theta_equil: the supposed theta in equilibrium. 7/18 pi for 70 degrees
    :param J_para: coupling constant along dimer rows
    :param J_perp: coupling constant across dimer rows
    :param h: strenght of the symmetry breaking field
    :return: approximated p
    """
    p = 1 / theta_equil * np.tan(theta_equil * p_guess) - p_guess
    q = 4 / theta_equil * np.sin(4 * theta_equil) / np.cos(theta_equil * p_guess) * (J_para + J_perp) / h
    p_plus = - (1/2) * p  +  np.sqrt((p/2) ** 2 - q)
    p_minus = - (1/2) * p  -  np.sqrt((p/2) ** 2 - q)
    return p_plus, p_minus

def accept_xi_inv_fit(reg, xmin, xmax, Tc_est, tolerance):
    xi_ampl = 1 / reg.slope
    Tc = - reg.intercept * xi_ampl

    return (Tc_est * (1-tolerance) < Tc < Tc_est * (1 + tolerance))

def best_fit_inv_old(T_arr, xi_inv_arr, Tc_est, tolerance, min_r_squared=0, min_points=3):
    """
    Searches the best linear fit through the data. The T_arr, xi_inv_arr should be linear
    around the critical temperature
    :param T_arr: The temperature values
    :param xi_inv_arr: The corresponding inverse correlation lengths
    :param Tc_est: the estimated critical temperature
    :param tolerance: the tolerance for the Tc that we get out of the fit. It
                        should lie inside [Tc_est - tolerance * Tc_est, Tc_est + tolerance Tc_est]
    :param min points: minimum number of points used in the regression
    :return:
    """
    # okay we want to fit the area above Tc so we just start to with all
    # values, then omit the first one and so on
    # But should we not also exclude values that are to far away and also not linear
    # anymore?

    # We try it with the linear regression stuff that chatgpt suggested because
    # If I use the naive mse, the fit will probably favor less points
    # I would have to install some kind of benefit for including more points
    best_starting_pos = 0
    best_ending_pos = 0
    # The minimum r_squared value that we need is the starting value for the best_r_sqaured
    best_r_squared = min_r_squared      # This value ranges to 1 I think and 1 is a really good fit?
    best_reg = None
    for starting_pos in range(len(T_arr) - min_points + 1):
        for ending_pos in range(starting_pos + min_points, len(T_arr)+1):
            # We should have at least 4 points i would say.
            T_fit = T_arr[starting_pos:ending_pos]
            xi_inv_fit = xi_inv_arr[starting_pos:ending_pos]

            reg = linregress(T_fit, xi_inv_fit)
            xi_ampl = 1 / reg.slope
            Tc = - reg.intercept * xi_ampl
            # We only accept the outcome if we get the expected result (haha)
            # The critical temperature that we get out of the fit has to bein +-10% range of the one that we calculated with the binder cumulant
            if ((reg.rvalue ** 2 > best_r_squared) and
                    (Tc_est * (1-tolerance) < Tc < Tc_est * (1 + tolerance))):
                best_starting_pos = starting_pos
                best_ending_pos = ending_pos
                best_r_squared = reg.rvalue ** 2
                best_reg = reg

    return best_reg, best_starting_pos, best_ending_pos


def min_x_accept_function(reg, starting_tau, ending_tau, min_tau):
    #print("starting_tau = ", starting_tau, " vs min tau = ", min_tau)
    return starting_tau >= min_tau

def best_lin_reg(x, y, min_r_squared, min_points=4, accept_function=None,
                 accept_function_args=None, more_points=True, require_end=False):
    x = np.array(x)
    y = np.array(y)
    best_starting_pos = 0
    best_ending_pos = 0
    # The minimum r_squared value that we need is the starting value for the best_r_sqaured
    best_r_squared = min_r_squared     # This value ranges to 1 I think and 1 is a really good fit?
    if more_points:
        best_r_squared *=1
    best_reg = None

    if min_points == 0:
        min_points = len(x)

    for starting_pos in range(len(x) - min_points + 1):
        for ending_pos in range(starting_pos + min_points, len(x)+1):
            # We should have at least 4 points i would say.
            if require_end:
                if ending_pos != len(x):
                    continue
            x_fit = x[starting_pos:ending_pos]
            y_fit = y[starting_pos:ending_pos]

            reg = linregress(x_fit, y_fit)

            if accept_function:
                accepted = accept_function(reg, x[starting_pos], x[ending_pos-1], *accept_function_args)
            else:
                accepted = True
            # We only accept the outcome if we get the expected result (haha)
            # The critical temperature that we get out of the fit has to bein +-10% range of the one that we calculated with the binder cumulant

            # reg.rvalue ** 2 seems to be a good measure of the linearity, but I somehow want another method of judging
            # which can incorporate the number of points used, so rather takes values that have more points
            # how would we do that? if we just divide reg.rvalue by the number of points, so basically a rvalue per point
            # we will always tend to using more points and loose linearity. Maybe by the root of the number of points?
            if more_points:
                quality_value = reg.rvalue ** 2 + (0.001 * (ending_pos - starting_pos - min_points)) # small bonus for fits that have more than 4 points
            else:
                quality_value = reg.rvalue ** 2
            if ((quality_value > best_r_squared) and accepted):
                best_starting_pos = starting_pos
                best_ending_pos = ending_pos
                best_r_squared = quality_value
                best_reg = reg
            #else:
                # if not accepted:
                #     print("starting_pos = ", starting_pos, "ending_pos = ", ending_pos, "not accepted")

    return best_reg, best_starting_pos, best_ending_pos

def best_fit_inv(T_arr, xi_inv_arr, Tc_est, tolerance, min_r_squared=0, min_points=6, more_points=False):
    return best_lin_reg(T_arr, xi_inv_arr, min_r_squared, min_points, accept_function=accept_xi_inv_fit,
                        accept_function_args=(Tc_est, tolerance), more_points=more_points)

def fold(x, y, fold=3):
    x = np.array(x)
    y = np.array(y)
    x_fold = []
    y_fold = []
    nr_points = len(x)
    fold = int(fold)
    nr_folded_points = nr_points // fold
    for point_nr in range(nr_folded_points):
        x_avg = 0
        y_avg = 0
        for nr_in_fold in range(fold):
            ind = point_nr * fold + nr_in_fold
            x_avg += x[ind]
            y_avg += y[ind]
        x_avg /= fold
        y_avg /= fold

        x_fold.append(x_avg)
        y_fold.append(y_avg)

    return np.array(x_fold), np.array(y_fold)

def equilibrium_angle_equation(theta, J_para, J_perp, h, p):
    zero = p * np.sin(theta * p) + 4 * np.sin(4 * theta) * (J_para + J_perp) / h
    return zero

def get_equilibrium_position(J_para, J_perp, h, p, theta_guess=0.8):
    """
    Calculates the equilibrium position of my XY model
    :param J_para: coupling parameter in parallel direction
    :param J_perp: coupling parameter in perpendicular direction
    :param h: strenght of the external field
    :param p: multiplicity of the external field
    :return: theta*
    """
    theta_equil = fsolve(equilibrium_angle_equation, theta_guess, args=(J_para, J_perp, h, p))
    return theta_equil

def find_time_of_temperature_change(csv_file):
    previous_temperature = None
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) < 2:  # Skip empty lines or lines with fewer than 2 columns
                continue
            time, temperature = float(row[0]), float(row[1])
            if previous_temperature is not None and temperature < previous_temperature * 0.99:
                return previous_time  # Return time of prior line
            previous_time, previous_temperature = time, temperature
    return None  # Return None if no temperature change is found

def meanAbsDiff(f):
    diffs = np.ediff1d(f)
    return np.mean(np.abs(diffs))

def getMovingFactor(f, avg_f, fraction=0.8):
    f_end = avg_f
    recent_ind = int(fraction * len(f))
    f_start = np.mean(f[:recent_ind])

    mean_abs_delta = meanAbsDiff(f)

    if mean_abs_delta == 0:
        # This should usually not happen, only if we are in low temp and
        # really nothing happens at all.
        moving_factor = 0
    else:
        moving_factor = np.abs(f_end - f_start) / ((len(f) - recent_ind) * mean_abs_delta)

    return moving_factor

def integrated_autocorr_time(f, ds):
    # We are looking at a file at which we did not write donw the autocorrelation time yet
    autocorr_function = st.acf(f, nlags=len(f) - 1)
    # I guess I have to rebuild the windowing algorithm here because otherwise
    # This acf thing only returns 40 values somehow, is this a bad sign?
    autocorr_time = 0
    for M in range(10, len(f)):
        # M is the cut
        cut_autocorr_function = autocorr_function[:M]
        lags = ds * np.arange(len(cut_autocorr_function))
        autocorr_time_M = np.trapz(cut_autocorr_function, lags)
        if (M * ds >= 5 * autocorr_time_M):
            autocorr_time = autocorr_time_M
            break
        autocorr_time = autocorr_time_M

    if np.isnan(autocorr_time):
        # if we again have only ones, wait but we dont for the files that I see here
        # Okay but if is really constant, the autocorrtime should be 0 I guess
        # or 1 I guess, I dont know
        autocorr_time = 1

    return autocorr_time

def get_largest_value_with_linestyle(ax, linestyle):
    largest_value = None

    for line in ax.lines:
        line_linestyle = line.get_linestyle()
        if line_linestyle == linestyle:
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            max_y_value = max(y_data)

            if largest_value is None or max_y_value > largest_value:
                largest_value = max_y_value

    return largest_value

def get_smallest_value_with_linestyle(ax, linestyle):
    largest_value = None

    for line in ax.lines:
        line_linestyle = line.get_linestyle()
        if line_linestyle == linestyle:
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            max_y_value = min(y_data)

            if largest_value is None or max_y_value < largest_value:
                largest_value = max_y_value

    return largest_value

def prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_k_fit, peak_cut_threshold, set_fts_to_zero,
                     min_points_fraction, nr_additional_cuts=0):
    p_k = get_frequencies_fftw_order(len(ft_k_fit))
    if cut_zero_impuls:
        p_k, ft_k_fit = cut_around_zero(p_k, ft_k_fit, nr_additional_cuts)
    if set_fts_to_zero:
        ft_k_fit -= np.min(ft_k_fit)
    if cut_around_peak:
        p_k, ft_k_fit = cut_data_around_peak(p_k, ft_k_fit, threshold_fraction=peak_cut_threshold,
                                             min_points_fraction=min_points_fraction)
    return ft_k_fit, p_k

def calc_structure_factor(filepath, zero_padding=False, zero_padding_center=False):
    #df = read_csv(filepath)
    df = read_large_df(filepath)
    filepath = Path(filepath)
    parameters = read_parameters_txt(str(filepath.with_suffix(".txt")))
    Lx = int(parameters["subsystem_Lx"])
    Ly = int(parameters["subsystem_Ly"])
    dim_size_x = int(parameters["dim_size_x"])
    if zero_padding:
        result_dim_y = 8 * Ly
        result_dim_x = 8 * Lx
    else:
        result_dim_y = Ly
        result_dim_x = Lx

    nr_subsystems = int(parameters["nr_subsystems"])
    struct_fact_k = {}
    struct_fact_l = {}
    for line in df:
        t = line[0]
        sin_ft = np.zeros((nr_subsystems, result_dim_y, result_dim_x), dtype=np.complex128)
        cos_ft = np.zeros((nr_subsystems, result_dim_y, result_dim_x), dtype=np.complex128)
        lattice = line[2:]      # or is this again to easy for f*ing pandas
        lattice = reshape_line(lattice, dim_size_x, Ly)
        for cell_nr in range(nr_subsystems):
            cell = extract_subsystem_from_matrix(lattice, Lx, cell_nr)
            cell = chess_board_trafo_rectangular_subsystems(cell, Lx, Ly)
            if zero_padding:
                new_cell_zeros = np.zeros((result_dim_y, result_dim_x), dtype=np.complex128)
                if zero_padding_center:
                    start_row = result_dim_y // 2 - Ly // 2
                    start_col = result_dim_x // 2 - Lx // 2
                    end_row = start_row + Ly
                    end_col = start_col + Lx
                    new_cell_zeros[start_row:end_row, start_col:end_col] = cell
                else:
                    new_cell_zeros[:Ly, :Lx] = cell
                cell = new_cell_zeros
            #print("cell:\n", cell)
            sin_cell = np.sin(cell) + 0j
            cos_cell = np.cos(cell) + 0j
            sin_trafo = (np.fft.fft2(sin_cell)) #np.fft.fftshift
            #print("sin output: \n", sin_trafo)
            cos_trafo = (np.fft.fft2(cos_cell)) #np.fft.fftshift
            #print("cos output: \n", cos_trafo)
            sin_ft[cell_nr] = sin_trafo
            cos_ft[cell_nr] = cos_trafo

        sin_ft_squared = np.zeros((result_dim_y, result_dim_x))
        cos_ft_squared = np.zeros((result_dim_y, result_dim_x))
        for cell_nr in range(nr_subsystems):
            #print("abs sin ft")
            #print(np.abs(sin_ft[cell_nr]))
            sin_ft_squared += np.abs(sin_ft[cell_nr]) ** 2
            #print("abs cos ft")
            #print(np.abs(cos_ft[cell_nr]))
            cos_ft_squared += np.abs(cos_ft[cell_nr]) ** 2
        # Averaging?
        # sin_ft_squared /= nr_subsystems
        # cos_ft_squared /= nr_subsystems

        st_fact_k = np.zeros(result_dim_y)
        st_fact_l = np.zeros(result_dim_x)

        # now for the two directions
        for l in range(result_dim_x):
            for k in range(result_dim_y):
                st_fact_k[k] += sin_ft_squared[k, l] + cos_ft_squared[k, l]
                st_fact_l[l] += sin_ft_squared[k, l] + cos_ft_squared[k, l]

        #print("ft_k_squared of this csv file")
        #print(st_fact_k)

        # averaging
        st_fact_k /= (Lx ** 4 * nr_subsystems)
        st_fact_l /= (Ly ** 4 * nr_subsystems)

        # reordering or something
        #st_fact_k = np.concatenate((np.array([st_fact_k[0]]), st_fact_k[1: len(st_fact_k) // 2 + 1][::-1], st_fact_k[len(st_fact_k) // 2 + 1:][::-1]))
        #st_fact_l = np.concatenate((np.array([st_fact_l[0]]), st_fact_l[1: len(st_fact_l) // 2][::-1],
        #                            st_fact_l[len(st_fact_l) // 2:][::-1]))

        struct_fact_k[t] = st_fact_k
        struct_fact_l[t] = st_fact_l
    return struct_fact_k, struct_fact_l


def dft2d(mat):
    """
    computes the classic 2D dft
    :param mat:
    :return:
    """
    transform = np.zeros_like(mat, dtype=np.complex128)

    M = mat.shape[0]
    N = mat.shape[1]
    for k in range(mat.shape[0]):
        for l in range(mat.shape[1]):
            for m in range(mat.shape[0]):       # shape 0 is in y direction, or lets say in m direction
                for n in range(mat.shape[1]):
                    transform[k, l] += mat[m, n] * np.exp(-2j * np.pi * m * k / M) * np.exp(-2j * np.pi * n * l / N)

    return transform


def get_2nd_moment_corr_length(ft):
    """
    calculates the 2nd moment correlation lenght from the fourier transform (haha it doesnt) Implementation from
    Stackoverflow, cannot be correct. Assumes order
    :param ft:
    :return:
    """
    L = len(ft)
    ft_zero = ft[0]
    ft_smallest_wave = ft[1]

    return  L / ( 2 * np.pi) * np.sqrt(ft_zero / ft_smallest_wave - 1)


def zoom_plot(ax, new_xlim):
    # Get the data from the original plot
    lines = ax.lines
    labels = [line.get_label() for line in lines]
    xdata = lines[0].get_xdata()
    ydata = lines[0].get_ydata()

    # Create a new figure and axis
    fig, ax_zoomed = plt.subplots()

    # Plot the data in the zoomed-in axis
    ax_zoomed.plot(xdata, ydata)

    # Set new x-limits
    ax_zoomed.set_xlim(new_xlim)

    # Set labels
    ax_zoomed.set_xlabel(ax.get_xlabel())
    ax_zoomed.set_ylabel(ax.get_ylabel())
    ax_zoomed.set_title(ax.get_title())

    # Return the new axis
    return fig, ax_zoomed

def write_dict_to_file(dictionary, filepath):
    if not filepath.endswith('.txt'):
        print("Filepath must end with '.txt'")
        return

    try:
        with open(filepath, 'w') as file:
            for key, value in dictionary.items():
                file.write(f"{key},{value}\n")
        print("Dictionary has been written to the file successfully.")
    except IOError as e:
        print(f"An error occurred while writing to the file: {e}")

def find_size_folders(directory):
    integer_folders = []
    if not os.path.isdir(directory):
        print("Invalid directory path.")
        return integer_folders

    for folder in os.listdir(directory):
        if os.path.isdir(os.path.join(directory, folder)):
            try:
                int(folder)  # Try converting folder name to an integer
                integer_folders.append(folder)
            except ValueError:
                continue
    return sorted(integer_folders, reverse=True)

def central_difference_wrong(y, x):
    """
    Calculate the central difference of order n.

    Parameters:
    y (function): The function to differentiate.
    x (np.array): The x values of the function.
    n (int): The order of the derivative.

    Returns:
    np.array: The central difference of order n.
    """
    n = len(y) // 2
    h = x[1] - x[0]  # step size
    diff = np.zeros_like(y)

    for i in range(n, len(x) - n):
        for j in range(1, n+1):
            diff[i] += ((-1)**j / (2**j)) * (y[i+j] - y[i-j])
        diff[i] /= h**n
    return diff

def central_difference(y, x):
    """
    Calculate the central difference of order n.

    Parameters:
    y (function): The function to differentiate.
    x (np.array): The x values of the function.
    n (int): The order of the derivative.

    Returns:
    np.array: The central difference of order n.
    """
    if len(y) % 2 == 0:
        print("EVEN NUMBERS NOT IMPLEMENTED")
        print("y ", y)
        print("x ", x)
        pass
    else:
        n = len(y) // 2     # ''accuracy''
        h = x[1] - x[0]  # step size

        coeff = np.zeros(len(y))
        values = np.zeros(len(y))

        for i in np.arange(-n, n + 1):
            ind = i + n
            if i == 0:
                coeff[ind] = 0
            else:
                coeff[ind] = central_difference_coefficients(2 * n, i)
            values[ind] = y[ind]
        summe = np.sum(coeff * values)
        diff = summe / h
        return diff

def central_difference_coefficients(n, p):
    num = (-1.0) ** (p + 1) * (math.factorial(n) ** 2)
    dom = p * math.factorial((n - p)) * math.factorial(n + p)
    return num / dom


def calculate_angle_from_slopes(slope1, slope2):
    # Calculate angle between lines using arctan
    # angle_radians = math.atan(abs((slope2 - slope1) / (1 + slope1 * slope2)))
    angle_radians = np.arctan(slope1) - np.arctan(slope2)
    # Convert angle to degrees
    # angle_degrees = math.degrees(angle_radians)
    print("angle in radians ", angle_radians)
    angle_degrees = np.abs(angle_radians / np.pi * 180)
    return angle_degrees
def angle_between_curves(x_values, y1_values, y2_values, x_star):
    """
    Calculate the angle between two curves y1(x) and y2(x) at their intersection point x_star.

    Parameters:
        x_values (array-like): Array of x values.
        y1_values (array-like): Array of y values corresponding to y1(x).
        y2_values (array-like): Array of y values corresponding to y2(x).
        x_star (float): Intersection point of the two curves.

    Returns:
        float: Angle between the curves at the intersection point in radians.
    """
    # Find the index corresponding to the intersection point x*
    index_x_star = np.argmin(np.abs(x_values - x_star))

    # Calculate derivatives using finite differences
    dy1_dx = np.gradient(y1_values, x_values)
    dy2_dx = np.gradient(y2_values, x_values)

    slope1 = dy1_dx[index_x_star]
    slope2 = dy2_dx[index_x_star]
    # Calculate angle in degrees from previous routine
    angle_degrees = calculate_angle_from_slopes(slope1, slope2)

    return angle_degrees


def main():
    print("This file is made to import, not to execute")

if __name__ == "__main__":
    main()
