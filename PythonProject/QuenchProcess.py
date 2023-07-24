import numpy as np
from FunctionsAndClasses import *
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker


def main():
    root = "../../Generated content/Quenching/Testing"
    name = "quench.process"
    png_name = "quench.png"
    root_dirs = os.listdir(root)
    print(root_dirs)
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
                T_start = paras["starting_T"]
                T_end = paras["end_T"]
                print(paras)
                t = np.array(df["t"])
                # normalizing so the quench begins at zero:
                t = t - t_eq
                t_tau = t / tau
                T = np.zeros(len(t_tau))
                # calcing T
                for i, t_eq_tau in enumerate(t_tau):
                    if T_start - t_eq_tau < T_end:
                        T[i] = T_end
                    elif t_eq_tau > 0:
                        T[i] = T_start - t_eq_tau
                    else:
                        T[i] = T_start

                xi_x = np.array(df["xi_x"][1:])
                xi_y = np.array(df["xi_y"][1:]) + 1/4 * np.max(xi_x)


                fig, axes = plt.subplots(1, 1)
                axes.plot(t_tau[1:], xi_x, ls="", marker=".", label=r"$\xi_x$", ms=2.5)
                axes.plot(t_tau[1:], xi_y, ls="", marker=".", label=r"$\xi_y$", ms=2.5)

                # Setze Tickmarken und Labels
                axes.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
                axes.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

                maxt_tau = int(np.max(t_tau))
                axes.xaxis.set_major_locator(ticker.MultipleLocator(base=maxt_tau / 4))
                axes.xaxis.set_minor_locator(ticker.MultipleLocator(base=maxt_tau / 4 / 5))
                # TODO minor locator muss
                axes.yaxis.set_minor_locator((plt.MultipleLocator(0.01)))
                # Füge Gitterlinien hinzu
                axes.grid(which='major', linestyle='--', alpha=0.5)
                # second x axis for the temperature
                axes2 = axes.twinx()
                axes2.plot(t_tau, T, c="red", label="T", alpha=0.3)
                # Setze Tickmarken und Labels
                axes2.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
                axes2.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
                axes2.yaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
                axes2.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
                axes2.set_ylabel("T")

                axes.set_ylabel(r"$\xi$")
                axes.set_xlabel(r"t$/ \tau_Q$")
                axes.set_title("Quench protocol")
                fig.legend()
                plt.show()

if __name__ == "__main__":
    main()