import os
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from decimal import Decimal
from glob import glob
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



markers = ["o", "s", "^", "v", "D", "p", "1", "2","*", "x", "+", "v", "^"]
colors = ["#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356", "#00305d", "#006ab2", "#009de0", "#00893a", "#65b32e", "#94C356"]
colors += colors + colors + colors + colors
colors += colors + colors + colors + colors

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


def process_folder_avg_fast(folderpath, file_ending):
    all_m = []
    cum = []
    value_files = glob(f"{folderpath}/*.{file_ending}")

    for i, file_path in enumerate(value_files):
        df, t = read_large_df_with_times(file_path, skiprows=1, sep=";")
        # this df is a list of lists, we dont have the times at all
        # I think I dont need the times, I just want to calculate the difference variation
        if i == 0:
            for m_t in df:
                all_m.append(m_t)
        else:
            for j, m_t in enumerate(df):
                all_m[j] += m_t
    for m_t in all_m:
        m2 = np.mean(np.array(m_t) ** 2)
        m4 = np.mean(np.array(m_t) ** 4)
        U_L = m4 / m2 ** 2
        cum.append(U_L)
    times = np.array(t)
    return cum, times

def process_folder_avg_balanced(folderpath, file_ending):
    cum = []
    times = []
    value_files = glob(f"{folderpath}/*.{file_ending}")

    para_path = find_first_txt_file(folderpath)
    parameters = read_parameters_txt(para_path)
    max_nr_m_values = 2e8   # maximum number of m values that can be processed in one batch

    try:
        nr_m_values = parameters["nr_mag_values"]
        subsystems_per_file = parameters["nr_subsystems"]
        total_nr_mag_values = nr_m_values * subsystems_per_file * len(value_files)
        nr_splits = int(math.ceil(total_nr_mag_values / max_nr_m_values))      # number of splits we have to spilt the run into because otherwise ram is to full
        t_vals_per_split = int(nr_m_values / nr_splits)

        for split in range(nr_splits):
            all_m = []
            for i, file_path in enumerate(value_files):
                df, t = read_large_df_with_times(file_path, skiprows=1 + split * t_vals_per_split,
                                                 sep=";", max_row=(split+1) * t_vals_per_split + 1)
                # this df is a list of lists, we dont have the times at all
                # I think I dont need the times, I just want to calculate the difference variation
                if i == 0:
                    for m_t in df:
                        all_m.append(m_t)
                else:
                    for j, m_t in enumerate(df):
                        all_m[j] += m_t
            # if all m is so large we need to do the calculations and reset it
            for m_t in all_m:
                m2 = np.mean(np.array(m_t) ** 2)
                m4 = np.mean(np.array(m_t) ** 4)
                U_L = m4 / m2 ** 2
                cum.append(U_L)
            times += t
        times = np.array(times)
    except KeyError:
        cum, times = process_folder_avg_fast(folderpath, file_ending)

    average_path = f"{folderpath}/this.cum_avg"
    header ="t,U_L\n"
    print("trying to write")
    write_lists_to_file(times, cum, filename=average_path, header=header)

    return cum, times

def read_folder_avg(folderpath, file_ending):
    # this one just wants to read already calculated files
    averaged_files = glob(f"{folderpath}/*.cum_avg")
    if averaged_files:
        file_path = averaged_files[0]
        df = pd.read_csv(file_path, delimiter=",", index_col=False)
        cum = np.array(df["U_L"])
        times = df['t']
    else:
        cum, times = process_folder_avg_balanced(folderpath, file_ending)
    return cum, times

def get_folder_average(folderpath, file_ending="mag", value="U_L", process_folder_avg_func=process_folder_avg_balanced):
    """
    This is supposed to look at a file that has observables like corr length and average them if we have multiple files
    this one only works for same times
    :param folderpath: path to the simulation
    :param file_ending: file ending of the file that holds the observables
    :param value: name of the value in the file
    :return: averaged value array, should it maybe also return the time?
    """
    value_files = glob(f"{folderpath}/*.{file_ending}")
    times = np.array([])
    if file_ending == "cum":
        cum = np.array([])
        for i, file_path in enumerate(value_files):
            df = pd.read_csv(file_path, delimiter=",", index_col=False)
            this_cum = np.array(df[value])
            # I think I dont need the times, I just want to calculate the difference variation
            if i == 0:
                cum = this_cum
                times = df['t']
            else:
                cum += this_cum
        # averaging
        cum = cum / len(value_files)
        times = np.array(times)
        return cum, times
    else:
        cum, times = process_folder_avg_func(folderpath, file_ending)
        # averaging
        times = np.array(times)
        return cum, times
def get_results_time_resolved(sizes, simulation_path, Tc=None, file_ending="mag", value_name="U_L",
                              process_folder_avg_func=process_folder_avg_balanced):
    size_cum_dic = {}
    size_times_dic = {}
    for size in sizes:
        size_path = os.path.join(simulation_path, str(size))
        if os.path.isdir(size_path):
            if not Tc:
                for temppath in os.listdir(size_path):
                    if temppath != "plots" and temppath[0] != ".":
                        Tc = float(temppath)
                        break
            temp_path = os.path.join(size_path,
                                     f"{Tc:.6f}")  # we only want to look at the current used critical temperature
            if os.path.isdir(temp_path):
                # okay so here we extract the cumulant for the a certain size temp combination
                # I guess we use the same tactic as in the ccumovertime script since I think there is no simpler
                # method to extract both the time and the values
                # we can use get folder avage here
                cum, times = get_folder_average(temp_path, file_ending=file_ending, value=value_name,
                                                process_folder_avg_func=process_folder_avg_func)
                size_cum_dic[size] = cum
                size_times_dic[size] = times
    return size_cum_dic, size_times_dic

def fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=3, xlim=0.5, config=None):
        z_list = np.linspace(1, 3, 201)
        sizes = list(size_cum_dic.keys())
        for i, size in enumerate(sizes):
            # We prbably want a 'base fold' for the smallest size
            cum = size_cum_dic[size]
            times = size_times_dic[size]

            t_fold, cum_fold = fold(times, cum, fold=fold_nr)
            # plot the points
            ax.plot(t_fold, cum_fold, **get_point_kwargs_color(colors[2 * i], markeredgewidth=1.5), marker=markers[i],
                    markersize=8, label=f"$L_\parallel$={size}", alpha=0.5)
            # the new interploation works with np.interp
            # the number of points of interp should be... i dont know at least a bit larger than the number of folded values?
            # what even happens when you try to use less values haha?
            # we could just use base_fold * len(t_fold) to get the same number of values we had before but folded?
            # The thing is for the error calculation we need the same interval for both sizes, I guess here it doesnt matter to much
            t_inter_plot = np.linspace(np.min(t_fold), np.max(t_fold), fold_nr * len(t_fold))
            #cum_inter_plot = np.interp(t_inter_plot, t_fold, cum_fold)

            #ax.plot(t_inter_plot, cum_inter_plot, color=colors[2 * i], alpha=0.5)
            # okay so much for the plotting, now the refitting
            # if there is a larger size we want to map this size onto it
            if i < len(sizes) - 1:
                next_size = sizes[i + 1]
                # cumulant and time values of the next size:
                cum_next = size_cum_dic[next_size]
                times_next = size_times_dic[next_size]
                # the folding shall be calculated based on the z that we use? so that t_fold and t_fold_next have
                # approximately the same number of points in the shared interval
                b = next_size / size  # next size is the larger size so b > 0
                # values to keep track wich z is the best
                best_z = 0
                best_msd = np.infty
                best_t_compare_arr = []
                best_cum_compare_arr = []
                # we need to try out every z
                for z in z_list:
                    # first we resacle the time of the small size?
                    times_next_rescaled = times_next / (b ** z)
                    # now we can look vor the shared interval by times_next and times_rescaled
                    t_lower_limit = np.maximum(np.min(times_next_rescaled), np.min(
                        times_next))  # the lower limit is the greater minimum of the times
                    t_upper_limit = np.minimum(np.max(times_next_rescaled), np.max(times_next))

                    # between those bounds we need an array with dense t values to compare the two cum curves
                    # the number of values should be... i dont know how large, just large enough? The arrays we will use
                    # to calculate the interpolation will use folded values, so just len(times_next) or something like this?
                    t_compare_arr = np.linspace(t_lower_limit, t_upper_limit, len(times_next))
                    # now we need the interpolations, but we want to use the interpolation of folded values, so we fold?
                    # so now the folding of the previous will stay the same, but the folding of the next will change
                    # if we for example rescale by a factor of 4, the folding should be 4 times as large?
                    t_next_rescaled_fold, cum_next_fold = fold(times_next_rescaled, cum_next,
                                                               b ** z * fold_nr)  # one with just the interpolation of the next curve without rescaling (Constant and independent of z)
                    # oh maaaaaaan you did it wrong you wanted to rescale the larger size because it has more values which can be folded to get fewer fluctuations
                    cum_next_compare_arr = np.interp(t_compare_arr, t_next_rescaled_fold, cum_next_fold)
                    cum_compare_arr = np.interp(t_compare_arr, t_fold, cum_fold)
                    # Okay from them we can now calculate an error?
                    err_arr = cum_next_compare_arr - cum_compare_arr
                    msd = np.mean(err_arr ** 2)
                    if msd < best_msd:
                        best_z = z
                        best_msd = msd
                        best_t_compare_arr = t_compare_arr
                        best_cum_compare_arr = cum_next_compare_arr
                # If we did this we want to plot it

                ax.plot(best_t_compare_arr[1:],
                        best_cum_compare_arr[1:],
                        linestyle="-", label=rf"${next_size} \rightarrow {size}$,"
                                             f"\tz = {best_z:.3f}", color=colors[2 * (i)], linewidth=2)
        ax.set_ylabel(r"$U_L$")
        ax.set_xlabel("$t \,/\, I$")
        #ax.set_xscale("log")
        ax.set_xlim(0, ax.get_xlim()[1] * xlim)

        if config is None:
            config = {
                "increasefontsize": 1.25,
                "labelhorizontalalignment": "right",
            }

        configure_ax(fig, ax, config)
        #ax.set_title("z extraction")
