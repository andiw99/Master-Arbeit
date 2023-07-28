from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker


def main():
    root = "../../Generated content/Coulomb/Binder"
    name = "binder.cumulants"
    root_dirs = os.listdir(root)
    print(root_dirs)
    errors_for_fit = False
    # arrays to save the xi corrsponding to T

    T_dic = {}
    cum_dic = {}
    # Loop through the directory contents and print the directories
    for size_dir in root_dirs:
        size_path = os.path.join(root, size_dir)
        if os.path.isdir(size_path):
            cum_path = size_path + "/"  + name
            n = get_size(size_path, temp_dirs = os.listdir(size_path))

            df = pd.read_csv(cum_path, delimiter=",", index_col=False)
            print(df)
            T = df["T"]
            U = df["U"]
            U = np.array(U)[np.argsort(T)]
            T = np.sort(T)
            T_dic[n] = T
            cum_dic[n] = U


    fig, ax = plt.subplots(1, 1)

    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
    # TODO minor locator muss
    ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    for size in T_dic.keys():
        print(cum_dic[size])
        ax.plot(T_dic[size], cum_dic[size], ls="", marker="+")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_title("Corr Length depending on T")
    save_plot(root, "/cum.png")



    plt.show()


def get_size(size_path, temp_dirs):
    parameters = {}
    for dir in temp_dirs:
        if (os.path.isdir(os.path.join(size_path, dir))):
            dir_path = os.path.join(size_path, dir)
            file_paths = os.listdir(dir_path)
            for f in file_paths:
                if (os.path.splitext(f)[1] == ".txt"):
                    parameters = read_parameters_txt(os.path.join(dir_path, f))
    return np.sqrt(parameters["n"])


if __name__ == "__main__":
    main()