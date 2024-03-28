import numpy as np

import matplotlib.pyplot as plt
from FunctionsAndClassesLight import *
from glob import glob

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

def fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=3, xlim=0.5):
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
        ax.set_xlabel("t / I")
        #ax.set_xscale("log")
        ax.set_xlim(0, ax.get_xlim()[1] * xlim)

        config = {
            "increasefontsize": 0.75,
            "labelhorizontalalignment": "right",
        }

        configure_ax(fig, ax, config)
        #ax.set_title("z extraction")

def main():

    simpath = "../../Generated content/Final/z-measurement-small/h=5/5/z"

    Tc = 2.85
    fold = 40
    folder_avg_function=process_folder_avg_balanced

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_5_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=5/1/z"

    Tc = 1.975000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 7, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_02_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=0.2-2/1/z"

    Tc = 1.970000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(80, 208, 3, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_05_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=0.5/1/z"

    Tc = 1.970000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_1_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/1/z"

    Tc = 1.975000
    fold = 120
    folder_avg_function=process_folder_avg_balanced

    sizes = np.linspace(72, 120, 3, dtype=np.int64)
    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 2, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()


def plot_eta_1_h_05():
    simpath = "../../Generated content/Final/z-measurement-small/0.5/z"

    Tc = 1.731000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 2, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()



if __name__ == "__main__":
    main()