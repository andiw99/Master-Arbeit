from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
import matplotlib.ticker as ticker



def main():
    # parameters
    root = "../../Generated content/Performance Benchmark/"
    name = "runtimes"
    gpu_dirs = os.listdir(root)
    runtimes_map = {}
    system_sizes_map = {}
    sim_times_map = {}
    min_size = 64
    max_size = 2048
    size_factor = 0.85
    ms = 10
    colors = ["#94c356", "#00305D", "#009FE3"]
    fig, ax = plt.subplots(1, 1, figsize=(7 * size_factor, 4.8 * size_factor))
    fig2, ax2 = plt.subplots(1, 1, figsize=(7* size_factor, 4.8 * size_factor))

    i = 0
    for j, gpu_dir in enumerate(gpu_dirs):
        gpu_path = os.path.join(root, gpu_dir)
        if os.path.isdir(gpu_path):
            gpu = os.path.basename(gpu_dir)
            runtimes = []
            system_sizes = []
            sim_times = []
            # Loop through the directory contents and print the directories
            size_dirs = os.listdir(gpu_path)

            for item in size_dirs:
                dir_path = os.path.join(gpu_path, item)
                print(os.path.isdir(dir_path))
                if (item != "plots") & os.path.isdir(dir_path):
                    print(item)
                    # Create the full path to the item
                    # Check if the item is a directory
                    if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                        filename = dir_path + "/" + name
                        print("reading: ", filename)
                        df = pd.read_csv(filename, header=None)
                        runtime = int(df[0].iloc[0])
                        files = os.listdir(dir_path)
                        parameters = {}
                        for f in files:
                            # we take the first file to be the parameters
                            if(os.path.splitext(f)[1] == ".txt"):
                                parameters = read_parameters_txt(os.path.join(dir_path, f))
                        n = parameters["n"]
                        if n >= min_size ** 2:
                            runtimes.append(runtime)
                            system_sizes.append(parameters["n"])
                            sim_times.append(parameters["end_t"])

            runtime_per_site_per_sim_time = np.array(runtimes) / np.array(system_sizes) / np.array(sim_times)
            runtime_per_site_per_sim_time = runtime_per_site_per_sim_time[np.argsort(system_sizes)]
            runtimes_per_sim_time = np.array(runtimes) / np.array(sim_times)
            runtimes_per_sim_time = runtimes_per_sim_time[np.argsort(system_sizes)]
            system_sizes = sorted(system_sizes)

            print(system_sizes)
            print(runtimes_per_sim_time)

            ax.plot(system_sizes, runtime_per_site_per_sim_time, linestyle="", marker="x", label=gpu, ms=ms, c=colors[i])
            ax2.plot(system_sizes, np.array(runtimes_per_sim_time) / 1000 / 60, linestyle="", marker="x", ms=ms, c=colors[i],
                     label=gpu)
            i += 1
            runtimes_map[gpu] = runtimes
            system_sizes_map[gpu] = system_sizes
            sim_times_map[gpu] = sim_times

    #x_ticks = sorted(system_sizes)[1:]
    #x_tick_labels = [str(int(np.sqrt(size))) + " x " + str(int(np.sqrt(size))) for size in x_ticks]
    # I think the x ticks will be rather constant so I will hardcode them here
    x_ticks = [64 * 64, 128 * 128, 256 * 256, 512 * 512, 1024 * 1024, 2048 * 2048]
    x_tick_labels = ["64 x 64", "128 x 128", "256 x 256", "512 x 512", "1024 x 1024", "2048 x 2048"]
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("System size n")
    ax.set_ylabel(r"$\frac{T}{t}/ms$", rotation=0, ha="right")
    ax.set_title("Time per Site")
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.tick_params(which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=0, width=0, labelsize=9)
    ax.legend()


    ax2.set_xscale("log")
    ax2.set_xlabel("System size n")
    ax2.set_ylabel("T/m")
    ax2.set_title("Computation time $T$ for a fixed Simulation time $t=5000$")
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_tick_labels)
    ax2.grid(which='major', linestyle='--', alpha=0.5)
    ax2.tick_params(which='both', length=6, width=2, labelsize=9)
    ax2.tick_params(direction='in', which='minor', length=0, width=0, labelsize=9)
    ax2.legend()
    fig.tight_layout()
    fig2.tight_layout()
    fig.savefig(root + "/time-per-site.pdf", format="pdf", transparent=True)
    fig2.savefig(root + "/total-time.pdf", format="pdf", transparent=True)
    plt.show()


if __name__ == "__main__":
    main()