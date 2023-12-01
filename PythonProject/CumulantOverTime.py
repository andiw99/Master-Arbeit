from FunctionsAndClasses import *
import glob
from matplotlib import pyplot as plt

def main():
    root = "../../Generated content/Subsystems/Cumulant Test"
    ending = "cum"

    size_dirs =sorted(os.listdir(root))
    print(size_dirs)

    cum_map = {}
    t_map = {}
    for size_dir in size_dirs:
        size_path = os.path.join(root, size_dir)
        temp_dirs = os.listdir(size_path)
        for temp_dir in temp_dirs:
            temp_path = os.path.join(size_path, temp_dir)

            cum_files = glob.glob(f"{temp_path}/*.cum")
            para_path = glob.glob(f"{temp_path}/*.txt")[0]
            parameters = read_parameters_txt(para_path)
            nr_cum_values = int(parameters["nr_cum_values"])
            Lx = int(parameters["subsystem_Lx"])
            T = parameters["T"]
            print(cum_files)
            # averageing over the cumulant files
            cumulant = np.zeros(nr_cum_values+1)
            t_arr = np.zeros(nr_cum_values+1)
            for cum_path in cum_files:
                df = pd.read_csv(cum_path, delimiter=",", index_col=False)
                cum = df["U_L"]
                cumulant += cum
                t_arr = df["t"]
            # averaging
            cumulant /= len(cum_files)
            cum_map[(Lx, T)] = cumulant
            t_map[(Lx, T)] = t_arr

    fig, ax = plt.subplots(1, 1)

    ax.set_xlim(0, np.max(t_arr))
    z = 3
    prev_size_temp = 0
    for size_temp in cum_map.keys():
        ax.plot(t_map[size_temp], cum_map[size_temp], label=f"Lx={size_temp[0]},  T = {size_temp[1]}")
        if(prev_size_temp):
            actual_size = size_temp[0]
            b = actual_size / prev_size_temp[0]
            # prev is the smaller L so b > 0
            t_map[prev_size_temp] = b ** z * t_map[prev_size_temp]

            ax.plot(t_map[prev_size_temp], cum_map[prev_size_temp], linestyle="dotted", label=f"Lx={prev_size_temp[0]},  T = {prev_size_temp[1]}  rescaled")
        prev_size_temp = size_temp

    configure_ax(fig, ax)
    plt.show()


if __name__ == "__main__":
    main()