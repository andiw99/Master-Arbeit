from FunctionsAndClasses import *
import glob
from matplotlib import pyplot as plt


def fold(t, U_L, fold=3):
    t = np.array(t)
    U_L = np.array(U_L)
    t_fold = []
    U_L_fold = []
    nr_points = len(t)
    nr_folded_points = nr_points // fold
    for point_nr in range(nr_folded_points):
        t_avg = 0
        U_L_avg = 0
        for nr_in_fold in range(fold):
            ind = point_nr * fold + nr_in_fold
            t_avg += t[ind]
            U_L_avg += U_L[ind]
        t_avg /= fold
        U_L_avg /= fold

        t_fold.append(t_avg)
        U_L_fold.append(U_L_avg)

    return np.array(t_fold), np.array(U_L_fold)
def main():
    root = "../../Generated content/Silicon/Subsystems/z extraction/high temperature/small eta/"
    ending = "cum"
    value = "U_L"
    title = ""
    z_min_fit = 1.9
    z_max_fit = 2.2
    res = 200
    fold_points = 20

    size_dirs =sorted(os.listdir(root))
    print(size_dirs)

    t_upper_limit = 1600
    t_lower_limit = 0

    cum_map = {}
    t_map = {}
    for size_dir in size_dirs:
        size_path = os.path.join(root, size_dir)
        if os.path.isdir(size_path):
            temp_dirs = os.listdir(size_path)
            for temp_dir in temp_dirs:
                temp_path = os.path.join(size_path, temp_dir)
                print(temp_path)
                cum_files = glob.glob(f"{temp_path}/*.{ending}")
                para_path = glob.glob(f"{temp_path}/*.txt")[0]
                parameters = read_parameters_txt(para_path)
                nr_cum_values = int(parameters[f"nr_{ending}_values"])
                print("nr of values: ", nr_cum_values)
                Lx = int(parameters["subsystem_Lx"])
                T = parameters["T"]
                print(cum_files)
                # averageing over the cumulant files
                cumulant = np.zeros(nr_cum_values+1)
                t_arr = []
                nr_cum_values_arr = np.zeros(nr_cum_values+1)
                for cum_path in cum_files:
                    print(cum_path)
                    df = pd.read_csv(cum_path, delimiter=",", index_col=False)
                    cum = df[value]
                    print("len df value:", len(df[value]))
                    print(len(cumulant))
                    print("current cum")
                    print(cum)
                    # okay we do something wild now, if one cum file isnt finished yet or something, we add the finished vals
                    if len(cum) < len(cumulant):
                        print("Attention, following file is corrupted or not yet finished")
                        print(cum_path)
                        cumulant[0:len(cum)] += cum
                        nr_cum_values_arr[0:len(cum)] += 1
                    else:
                        cumulant += cum
                        nr_cum_values_arr += 1
                    if len(df["t"]) > len(t_arr):
                        print(f"len(df[t]) larger than len(t_arr): {len(df['t'])} > {len(t_arr)}")
                        t_arr = df["t"]
                # averaging
                cumulant = cumulant / nr_cum_values_arr
                cumulant = cumulant[cumulant > 0]           # cum is never zero so we can assume that the spaces with zero have not been reached
                print(cumulant)
                cum_map[(Lx, T)] = cumulant
                t_map[(Lx, T)] = t_arr

    fig, ax = plt.subplots(1, 1)

    if t_upper_limit:
        ax.set_xlim(0, t_upper_limit)
    else:
        ax.set_xlim(0, np.max(t_arr))

    z_list = [2.1, 2.2, 2.0]
    prev_size_temp = 0
    for size_temp in cum_map.keys():
        t_arr_plt = np.array(t_map[size_temp])
        cum_arr = np.array(cum_map[size_temp])
        print(len(cum_arr))
        print(len(t_arr_plt))
        ax.plot(t_arr_plt[t_arr_plt > t_lower_limit], cum_arr[t_arr_plt > t_lower_limit], label=f"Lx={size_temp[0]},  T = {size_temp[1]}")
        if(prev_size_temp):
            actual_size = size_temp[0]
            b = actual_size / prev_size_temp[0]
            # prev is the smaller L so b > 0
            for z in z_list:
                t_arr = np.array(b ** z * t_map[prev_size_temp])
                ax.plot(t_arr[t_arr > t_lower_limit], np.array(cum_map[prev_size_temp])[t_arr > t_lower_limit], linestyle="dotted", label=f"Lx={prev_size_temp[0]},  T = {prev_size_temp[1]}  rescaled z = {z}")
        prev_size_temp = size_temp
    if t_lower_limit:
        ax.set_xlim(t_lower_limit, ax.get_xlim()[1])
    configure_ax(fig, ax)
    plt.savefig(root + f"/{ending}-over-time.png", format="png")
    plt.show()

    # okay we now will try to evaluate the error between the two sets

    fig, ax = plt.subplots(1, 1)

    if t_upper_limit:
        ax.set_xlim(0, t_upper_limit)
    else:
        ax.set_xlim(0, np.max(t_arr))
    z_list = np.linspace(z_min_fit, z_max_fit, res)
    prev_size_temp = 0
    res = 1000
    print(cum_map.keys())
    sizes = [pair[0] for pair in cum_map.keys()]
    keys = np.array(list(cum_map.keys()))[np.argsort(sizes)]
    keys = [(key[0], key[1]) for key in keys]
    print(keys)
    for i, size_temp in enumerate(keys):
        t_arr_plt = np.array(t_map[size_temp])
        cum_arr = np.array(cum_map[size_temp])

        t_fold, U_L_fold = fold(t_arr_plt, cum_arr, fold=fold_points)
        ax.plot(t_fold[t_fold > t_lower_limit], U_L_fold[t_fold > t_lower_limit], linestyle='', marker="x", markersize=5,
                color=f"C{i}", label=f"$L_x$={size_temp[0]},  T = {size_temp[1]}")
        cum_inter = interp1d(t_fold,  U_L_fold)
        #ax.plot(t_fold[t_fold > t_lower_limit],
        #        cum_inter(t_fold[t_fold > t_lower_limit]),
        #        label=f"$L_x$={size_temp[0]},  T = {size_temp[1]}",
        #        color=f"C{i}")
        max_t = np.max(np.array(t_fold))
        min_t = np.min(np.array(t_fold))
        min_squared_error = np.infty
        min_z = 1
        if(prev_size_temp):
            actual_size = size_temp[0]
            b = actual_size / prev_size_temp[0]
            # prev is the smaller L so b > 0
            for z in z_list:
                pre_t_arr, cum_arr = fold(t_map[prev_size_temp], cum_map[prev_size_temp], fold=fold_points)
                pre_t_arr = b ** z * pre_t_arr
                pre_cum_inter = interp1d(pre_t_arr, cum_arr, kind="cubic")
                upper_t = np.minimum(np.max(pre_t_arr), max_t)
                lower_t = np.maximum(np.min(pre_t_arr), min_t)
                err_t_arr = np.linspace(lower_t, upper_t, res)
                err_arr = (pre_cum_inter(err_t_arr) - cum_inter(err_t_arr)) ** 2
                total_squared_error = np.sum(err_arr)
                if(total_squared_error < min_squared_error):
                    min_squared_error = total_squared_error
                    min_z = z
            plot_t_arr = b ** min_z * t_map[prev_size_temp]
            t_arr = np.array(plot_t_arr)
            ax.plot(t_arr[t_arr > t_lower_limit],
                    np.array(cum_map[prev_size_temp])[t_arr > t_lower_limit],
                    linestyle="-", label=f"$L_x$={prev_size_temp[0]},"
                                              f"  T = {prev_size_temp[1]}  rescaled z = {min_z:.3f}", color=f"C{i}")
        prev_size_temp = size_temp
    ax.set_ylabel(r"$U_L$")
    ax.set_xlabel("t")
    configure_ax(fig, ax)
    ax.set_title(title)
    plt.savefig(root + f"/{ending}-over-time-scan.png", format="png")
    plt.show()


if __name__ == "__main__":
    main()