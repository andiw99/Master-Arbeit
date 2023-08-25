import csv
import os
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib
from itertools import product

# import matplotlib; matplotlib.use("TkAgg")

colors = ["#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356", "#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356"]


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

def plot_multiple_times(filepath, config={"nr_of_meshs": 16, "cell_L": 128, "cell_nr": 1, "chess_trafo": 1, "nr_colorbar_ticks": 5}):

    para_filepath = os.path.splitext(filepath)[0] + ".txt"
    parameters = read_parameters_txt(para_filepath)

    # extract parameters out of config
    nr_of_meshs = config["nr_of_meshs"]
    cell_L = config["cell_L"]
    parameters["cell_L"] = cell_L
    nr_rows = parameters["nr_save_values"]
    # equidistant row numbers to use
    # Select the rows with the row equidistant row numbers
    df = []
    rows = np.linspace(0, nr_rows-1, nr_of_meshs, endpoint=True)
    rows = [int(row) for row in rows]
    print(rows)
    with open(filepath) as f:
        for i, line in enumerate(f):
            if i in rows:
                print("appending...")
                df.append(string_to_array(line[:-2]))

    print(df)

    fig, axes = plt.subplots(int(np.sqrt(nr_of_meshs)), int(np.sqrt(nr_of_meshs)), figsize=[12, 10])
    i = 0
    for axs in (axes):
        for ax in (axs):
            ind = i

            # read time and temperature and extract the values
            parameters["t"] = df[ind][0]
            parameters["T"] = df[ind][1]
            row = np.array(df[ind][2:])

            if config["cell_L"]:
                row = extract_cell(row, config["cell_nr"], config["cell_L"])
            im = plot_colormesh(ax, row, parameters, config)
            i += 1
    plt.tight_layout()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.83, 0.04, 0.02, 0.7])
    well_pos = np.sqrt(parameters["beta"] / 2)
    nr_ticks = config["nr_colorbar_ticks"]
    ticks = np.linspace(- well_pos, well_pos, nr_ticks, endpoint=True)
    tick_labels = np.linspace(-1, 1, nr_ticks, endpoint=True)
    tick_labels = [str(tick_label) for tick_label in tick_labels]
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=ticks)
    cbar.ax.set_yticklabels(tick_labels)

    textstr = ''
    wanted_paras=["dt", "J", "alpha", "beta", "eta", "tau", "Jy", "lat_dim"]
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
        chess_trafo = parameters["J"] < 0
        print("auto chess trafo yields: chess_trafo=", chess_trafo)
    if chess_trafo:
        print("doing chess trafo")
        row = chess_board_trafo(row)
    print(row)
    well_pos = np.sqrt(parameters["beta"] / 2)
    cf = ax.pcolormesh(row, cmap="copper", vmax=well_pos, vmin=-well_pos)

    return cf


def find_first_csv_file(folder_path):
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".csv"):
            return os.path.join(folder_path, file_name)

    return None

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
    wanted_paras=["cell_L", "T", "J", "alpha", "beta", "eta", "tau", "Jy", "lat_dim"]
    for key in wanted_paras:
        try:
            fname += key + "=" + str(paras[key])
        except:
            pass
    return fname


def list_folders_and_subfolders(directory_path):
    folder_list = []

    for root, dirs, files in os.walk(directory_path):
        for dir_name in dirs:
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


def configure_ax(fig, ax, config=None):
    """
    Takes a fig and an axes and configures the axes and stuff. If no config map is provided standard config is used
    :param fig:
    :param ax:
    :param config:
    :return: void
    """

    x_span, y_span = get_spans(ax)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=x_span/5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=x_span/5 / 5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=y_span/5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(base=y_span/5 / 5))
    # We want to have inline ticks
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

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
    ax.set_ylabel(ax.get_ylabel(), rotation=0, ha="right")
    #legend
    ax.legend()
    plt.tight_layout()


def mark_point(ax, x, y, c="C0"):
    ax.scatter(x, y, c=c)
    ax.plot([x, x], [0, y], c=c, ls='--', lw=2, alpha=0.75)
    ax.plot([0, x], [y, y], c=c, ls='--', lw=2, alpha=0.75)

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
    return x_span, y_span


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



def main():
    print("This file is made to import, not to execute")

if __name__ == "__main__":
    main()
