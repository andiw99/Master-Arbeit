import os
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt


markers = ["o", "s", "^", "v", "D", "p", "1", "2","*", "x", "+", "v", "^"]
colors = ["#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356", "#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356"]
colors += colors + colors + colors + colors
colors += colors + colors + colors + colors
def string_to_list(input_string):
    return [float(x) for x in input_string.split(',')]
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

def find_first_txt_file(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".txt"):
                return os.path.join(root, file)
    return None  # Return None if no .txt file is found

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

def write_lists_to_file(x, y, filename, header):
    with open(filename, 'w') as file:
        file.write(header)
        for x_val, y_val in zip(x, y):
            file.write(f'{x_val},{y_val}\n')


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


def get_point_kwargs_color(color, markeredgewidth=1):
    return {"linestyle": "None", "markerfacecolor": "none", "markeredgecolor": color, "markeredgewidth": markeredgewidth}

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
                  title_fontsize=config["legendfontsize"], alignment="left")
    if config["tight_layout"]:
        plt.tight_layout()


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
