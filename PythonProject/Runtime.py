from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
import matplotlib.ticker as ticker



def main():
    # parameters
    root = "../../Generated content/Performance Benchmark"
    name = "runtimes"
    root_dirs = os.listdir(root)

    print(root_dirs)
    runtimes = []
    system_sizes = []
    # Loop through the directory contents and print the directories
    for item in root_dirs[::-1]:
        if (item != "plots"):
            # Create the full path to the item
            dir_path = os.path.join(root, item)
            # Check if the item is a directory
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                print("reading: ", filename)
                df = pd.read_csv(filename, header=None)
                runtime = int(df[0].iloc[0])
                runtimes.append(runtime)
                files = os.listdir(dir_path)
                parameters = {}
                for f in files:
                    # we take the first file to be the parameters
                    if(os.path.splitext(f)[1] == ".txt"):
                        parameters = read_parameters_txt(os.path.join(dir_path, f))
                system_sizes.append(parameters["n"])


    runtime_per_site = np.array(runtimes) / np.array(system_sizes)

    plt.plot(system_sizes, runtimes, linestyle="", marker="x")
    fig, ax = plt.subplots(1, 1)
    ax.plot(system_sizes, runtime_per_site, linestyle="", marker="x", label="AMD RX 6800")
    x_ticks = sorted(system_sizes)[1:]
    x_tick_labels = [str(int(np.sqrt(size))) + " x " + str(int(np.sqrt(size))) for size in x_ticks]
    ax.set_xscale("log")
    #ax.set_yscale("log")
    ax.set_xlabel("System size n")
    ax.set_ylabel("T/ms")
    ax.set_title("Computation time $T$ per site for a fixed Simulation time $t=5000$")
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.set_ylim((0, np.max(runtime_per_site) * 1.05))
    ax.tick_params(which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=0, width=0, labelsize=9)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    main()