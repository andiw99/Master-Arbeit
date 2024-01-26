import csv
import os
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib
from itertools import product, combinations
from IPython import embed
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import fsolve


# import matplotlib; matplotlib.use("TkAgg")

colors = ["#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356", "#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356"]
markers = [".", "x", "+", "*", "D", "1", "2", "v", "^"]

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

def string_to_array(input_string):
    # Split the input string by commas and convert each element to a float
    # input_string = input_string.strip()
    values = [float(x) for x in input_string.split(',')]

    # Create a NumPy array from the list of float values
    array = np.array(values)

    return array


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

def read_large_df(filepath, rows=None):
    df = []
    with open(filepath) as f:
        for i, line in enumerate(f):
            if rows:
                if i in rows:
                    df.append(string_to_array(line[:-2]))
            else:
                df.append(string_to_array(line[:-2]))
    return df

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
    if nr_rows < nr_of_meshs:
        nr_of_meshs = nr_rows
        print(f"{nr_of_meshs} snapshots not available, using {nr_rows} instead.")
    # equidistant row numbers to use
    # Select the rows with the row equidistant row numbers
    rows = np.linspace(0, nr_rows - 1, nr_of_meshs, endpoint=True)
    rows = [int(row) for row in rows]
    df = read_large_df(filepath, rows)
    stretch = Lx/Ly
    if config["subsystem"]:
        stretch = np.minimum(parameters["subsystem_Lx"], cell_L) / parameters["subsystem_Ly"]
        print("stretch = " , stretch)
    elif cell_L != 0:
        stretch = 1
    if nr_of_meshs == 2:
        fig, axes = plt.subplots(2, 1, figsize=[2 + 10 * stretch, 10])
    else:
        fig, axes = plt.subplots(int(np.sqrt(nr_of_meshs)), int(np.sqrt(nr_of_meshs)), figsize=[2 + 10 * stretch, 10])
    i = 0
    for axs in (axes):
        if nr_of_meshs == 2:
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
                    subsystem_Lx = parameters["subsystem_Lx"]
                    subsystem_Ly = parameters["subsystem_Ly"]
                    print(np.sqrt(parameters["nr_subsystems"] * parameters["x_y_factor"]))
                    max_nr_subsystems_per_row = int(np.sqrt(parameters["nr_subsystems"] * parameters["x_y_factor"]) + 0.5)
                    max_nr_subsystems_per_col = int(parameters["nr_subsystems"] / max_nr_subsystems_per_row)
                    nr_subsystems_per_row = np.minimum(int(config["cell_L"] / subsystem_Lx), max_nr_subsystems_per_row)
                    if nr_subsystems_per_row == 0:
                        nr_subsystems_per_row = 1
                        nr_subsystems_per_col = 1
                    else:
                        nr_subsystems_per_col = np.minimum(int(config["cell_L"] / subsystem_Ly), max_nr_subsystems_per_col)
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
                    row = extract_cell_from_rectangle(row, config["cell_nr"], config["cell_L"], Lx, Ly)
            im = plot_rectangular_colormesh(ax, row, parameters, config)
            config["cell_L"] = actual_cell_L
            i += 1
    plt.tight_layout()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.04, 0.02, 0.7])
    well_pos = np.sqrt(parameters["beta"] / 2)
    nr_ticks = config["nr_colorbar_ticks"]
    if config["angle"] == 1:
        ticks = np.linspace(0, 2 * np.pi, nr_ticks, endpoint=True)
        tick_labels = np.linspace(0, 2 * np.pi, nr_ticks, endpoint=True)
        tick_labels = [f"{tick_label:.2f}" for tick_label in tick_labels]
    elif config["angle"] == 2:
        ticks = np.linspace(- np.pi / 2, np.pi / 2, nr_ticks, endpoint=True)
        tick_labels = np.linspace(-np.pi / 2, np.pi / 2, nr_ticks, endpoint=True)
        tick_labels = [f"{tick_label:.2f}" for tick_label in tick_labels]
    else:
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
    subsystem_Lx = int(parameters["subsystem_Lx"])
    subsystem_Ly = int(parameters["subsystem_Ly"])
    cell_L = int(config["cell_L"])
    if config["subsystem"]:
        pass
    elif cell_L != 0:
        row = row.reshape((cell_L, cell_L))
    else:
        row = row.reshape((Ly, Lx))
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
            if cell_L < Lx:
                row = chess_board_trafo_rectangular_subsystems(row, cell_L, Ly)
            else:
                row = chess_board_trafo_rectangular_subsystems(row, Lx, Ly)
        else:
            row = chess_board_trafo(row)
    if config["angle"] == 1:
        cf = ax.pcolormesh(row, cmap="viridis_r", vmax=2 * np.pi, vmin=0)
    elif config["angle"] == 1:
        cf = ax.pcolormesh(row, cmap="viridis_r", vmax=np.pi/2, vmin=-np.pi/2)
    else:
        well_pos = np.sqrt(parameters["beta"] / 2)
        cf = ax.pcolormesh(row, cmap="viridis_r", vmax=2 * well_pos, vmin=-2 * well_pos)

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
    df = pd.read_csv(filepath, delimiter=",", header=None, index_col=0)
    para_set = {}
    for label in df.index:
        try:
            para_set[label] = float(df.loc[label, 1])
        except ValueError:
            pass
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
    print(subsystems_per_row, "  ", subsystems_per_col)
    for row in range(subsystems_per_row):
        for col in range(subsystems_per_col):
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
    for i, digit in enumerate(str(number)):
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
        return config

PLOT_DEFAULT_CONFIG = {
    "ylabelsize": 11,
    "xlabelsize": 11,
    "titlesize": 14,
    "xtickfontsize": 9,
    "ytickfontsize": 9,
    "legendfontsize": 11,
    "ticklength": 6,
    "tickwidth": 2,
}

def delete_last_line(filename):
    # Read all lines from the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Remove the last line
    if len(lines) > 0:
        lines.pop()

    # Write the remaining lines back to the file
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

    x_span, y_span, (xmin, xmax, ymin, ymax) = get_spans(ax)

    if ax.get_yscale() != "log":
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=y_span/5))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=y_span/5 / 5))
    if ax.get_xscale() != "log":
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_span/5))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_span/5 / 5))

    # We want to have inline ticks
    ax.tick_params(direction='in', which='both', length=config["ticklength"], width=config["tickwidth"], labelsize=9)
    ax.tick_params(direction='in', which='minor', length=int(config["ticklength"] * 0.75), width=int(config["tickwidth"] * 0.75), labelsize=9)

    if ax.get_xscale() != "log":
        remove_origin_ticks(ax)

    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)

    # title, achsenbeschriftungen, legend
    get_functions = [ax.get_title, ax.get_xlabel, ax.get_ylabel]        # haha functions in a list, just python things
    set_functions = [ax.set_title, ax.set_xlabel, ax.set_ylabel]
    default_texts = [os.path.splitext(os.path.basename(sys.argv[0]))[0], "x", "y"]
    for i, get in enumerate(get_functions):
        if get() == "":
            # If we have empty string as title or stuff
            set_functions[i](default_texts[i])

    # rotate the y label
    ax.set_ylabel(ax.get_ylabel(), rotation=0, ha="right", fontsize=config["ylabelsize"])
    ax.set_xlabel(ax.get_xlabel(), fontsize=config["xlabelsize"])
    ax.set_title(ax.get_title(), fontsize=config["titlesize"])
    #legend
    ax.legend(fontsize=config["legendfontsize"])
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
    xmin = np.infty
    xmax = -np.infty
    ymin = np.infty
    ymax = -np.infty
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
        dif = np.infty
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
        dif = np.infty
        ind_i, ind_j = 0, 0
        for i,j in product(range(len(y)), range(len(z))):
            # we need a metric for the difference, lets take RS
            difference = np.sqrt((y[i]- z[j]) ** 2 + (x_y[i] - x_z[j]) ** 2)
            if difference < dif:
                dif =difference
                ind_i = i
                ind_j = j
        return (ind_i, ind_j)


def find_intersection(x_range, y1, y2, res=1000):
    # Interpolate the curves
    #print("diff_func:")
    x_inter = np.linspace(x_range[0], x_range[-1], res)
    diff_func = np.interp(x_inter, x_range, y1 - y2)
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
    # We could do this sign thing for every diff curve and then choose always the one that has the least total difference to the other sign arrays
    ind_signchange = np.where(signchange == 1)[0][0]
    intersection_x = x_inter[ind_signchange]

    return intersection_x

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
        new_threshold = 0.9 * threshold_fraction        # small steps?
        return cut_data_around_peak(x_values, y_values, new_threshold, min_points_fraction)
    # I think if we
    # Extract the subset of data around the peak
    x_cut = x_values[above_threshold_indices]
    y_cut = y_values[above_threshold_indices]

    # we have to add the offset back again
    y_cut += y_offset

    return x_cut, y_cut

def get_intersections(size_T_cum_dic):
    intersections = []
    sizes = list(size_T_cum_dic.keys())

    for (size1, size2) in combinations(sizes, 2):
        U_L_1 = size_T_cum_dic[size1]["U_L"]
        U_L_2 = size_T_cum_dic[size2]["U_L"]
        # We assume that the Temperatures are the same for the two curves...
        T_arr = size_T_cum_dic[size1]["T"]
        intersection = find_intersection(T_arr, U_L_1, U_L_2)
        intersections.append(intersection)

    return intersections

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
            system.append(rows[j][i *  row_length: (i + 1) * row_length])
    return np.array(system)
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

def process_file(file_path, threshold, key='t', value='U_L'):
    """
    Process a single file and calculate the average after the given threshold.
    """
    df = pd.read_csv(file_path)
    if 0 < threshold < 1:
        threshold = threshold * len(df[key])
    df = df[df[key] >= threshold]
    nr_values = df.shape[0]
    if nr_values < 1:
        average_value = 0
    else:
        average_value = df[value].mean()
    return average_value, nr_values


def process_size_folder(size_folder, threshold, key='T', value='U_L', file_ending='cum', selected_temperatures=None):
    """
    Process all files in a size folder and return a dictionary with temperature and average values.
    """
    eps = 1e-4
    result = {key: [], value: []}
    for temp_folder in os.listdir(size_folder):
        if selected_temperatures is not None:
            temp = float(temp_folder)
            selected = False
            # we have to look if approximately this temperature is present in selected_temperatures
            for sel_temp in selected_temperatures:
                if sel_temp - eps < temp < sel_temp + eps:
                    selected = True
            if not selected:
                continue
        if (temp_folder != "plots") & (temp_folder[0] != "."):
            temp_folder_path = os.path.join(size_folder, temp_folder)
            if os.path.isdir(temp_folder_path):
                temp_average = []
                nr_avg_values = []
                nr_subsystems = []
                for file_name in os.listdir(temp_folder_path):
                    if file_name.endswith(file_ending):
                        # it could be that some files have a cumulant averaged out of more subsystems, that shoudl be taken into consideration
                        file_path = os.path.join(temp_folder_path, file_name)
                        average_value, nr_values = process_file(file_path, threshold, 't', value)
                        para_file_path = os.path.splitext(file_path)[0] + ".txt"
                        parameters = read_parameters_txt(para_file_path)
                        nr_subsys = parameters["nr_subsystems"]
                        temp_average.append(average_value)
                        nr_avg_values.append(nr_values)
                        nr_subsystems.append(nr_subsys)
                if temp_average:
                    val_avg = (np.sum(np.array(temp_average) * np.array(nr_avg_values) * np.array(nr_subsystems)) /
                               np.sum(np.array(nr_avg_values) * np.array(nr_subsystems)))
                    result[key].append(float(temp_folder))
                    result[value].append(val_avg)

    result[value] = np.array(result[value])[np.argsort(result[key])]
    result[key] = np.sort(result[key])

    return result

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
            print(filepath, "firstfile = ", first_file)
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
    print(f"averaged {nr_files} files")
    for t in t_ft_k:
        t_ft_k[t] /= nr_files
        t_ft_l[t] /= nr_files
    return t_ft_k, t_ft_l

def average_lastline_ft(folderpath, ending=".ft"):
    # like average_ft but it only returns the ft for the latest t
    files = os.listdir(folderpath)
    t_ft_k = np.array([])
    t_ft_l = np.array([])
    first_file = True
    nr_files = 0
    for file in (files):
        if file.endswith(ending):
            filepath = os.path.join(folderpath, file)
            df = pd.read_csv(filepath, sep=";")
            print(filepath, "firstfile = ", first_file)
            if first_file:
                print("df['ft_k'].iloc[-1]")
                print(df["ft_k"].iloc[-1])
                t_ft_k = string_to_array(df["ft_k"].iloc[-1])
                t_ft_l = string_to_array(df["ft_l"].iloc[-1])
            else:
                t_ft_k += string_to_array(df["ft_k"].iloc[-1])
                t_ft_l += string_to_array(df["ft_l"].iloc[-1])
            first_file = False
            nr_files += 1
    # average
    print(f"averaged {nr_files} files")
    t_ft_k /= nr_files
    t_ft_l /= nr_files
    return t_ft_k, t_ft_l


def get_frequencies_fftw_order(nr_times):
    K = nr_times // 2
    freqs = []
    for i in range(nr_times):
        if i + K < nr_times:
            ind = i + K
        else:
            ind = i - K
        freqs.append(2 * np.pi * (ind - nr_times/2) / nr_times)
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
        popt, pcov = curve_fit(fitfunc, p, ft, sigma=errors)
        perr = np.sqrt(np.diag(pcov))

    except RuntimeError:
        print("RuntimeError in your function fit_lorentz. Maybe no fit possible because"
              " of to many fluctuations or because it is just flat. Returning 0")
        return (0, 0, 0), None
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

def cut_zero_imp(p, ft, ft_err=None):
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

def plot_struct_func(px, py, fx, fy, error_x=0, error_y=0):
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


def main():
    print("This file is made to import, not to execute")

if __name__ == "__main__":
    main()
