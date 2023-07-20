import numpy as np
from FunctionsAndClasses import *
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker


def main():
    root = "../../Generated content/Fit testing/"
    name = "quench.process"
    png_name = "quench.png"
    root_dirs = os.listdir(root)

    for item in root_dirs:
        if (item != "plots"):
            dir_path = os.path.join(root, item)
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                print("reading: ", filename)
                df = read_struct_func(filename)
                # We assume that we have one that is called 0.txt
                paras = read_parameters_txt(dir_path + "/0.txt")
                tau = paras["tau"]
                t_eq = paras["t_eq"]

                fig, axes = plt.subplots(1, 1)
                t = np.array(df["t"][1:])
                # normalizing so the quench begins at zero:
                t = t - t_eq
                t_tau = t / tau
                xi_x = np.array(df["xi_x"][1:])
                axes.plot(t_tau, xi_x, ls="", marker=".", label=r"$\xi_x$")
                xi_y = np.array(df["xi_y"][1:]) + 1/4 * np.max(xi_x)
                axes.plot(t_tau, xi_y, ls="", marker=".", label=r"$\xi_y$")
                # Setze Tickmarken und Labels
                axes.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
                axes.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
                axes.xaxis.set_major_locator(ticker.MultipleLocator(base=0.25))
                axes.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.05))
                # TODO minor locator muss
                axes.yaxis.set_minor_locator((plt.MultipleLocator(0.01)))
                # Füge Gitterlinien hinzu
                axes.grid(which='major', linestyle='--', alpha=0.5)

                axes.set_ylabel(r"$\xi$")
                axes.set_xlabel(r"t$/ \tau_Q$")
                axes.legend()
                plt.show()

if __name__ == "__main__":
    main()