import numpy as np
from FunctionsAndClasses import *
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as animation
from functools import partial
import matplotlib; matplotlib.use("TkAgg")


def main():
    root = "../../Generated content/New/Overdamped Quenching/"
    name = "quench.process"
    png_name = "quench.png"
    root_dirs = os.listdir(root)
    print(root_dirs)

    xi_avg_dic = {}
    t_tau_dic = {}

    # The problem could be that we cannot read in the df if we save like 1000 values.
    # We might have to consider building something that splits up the file in 100 MB portions?


    plt.show()

    for item in root_dirs:
        if (item != "plots"):
            dir_path = os.path.join(root, item)
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                print("reading: ", filename)
                df = read_struct_func(filename)
                sys_df = read_csv(dir_path + "/0.csv")
                print(sys_df.shape[0])
                nr_times = sys_df.shape[0]
                # We assume that we have one that is called 0.txt
                paras = read_parameters_txt(dir_path + "/0.txt")
                tau = paras["tau"]
                t_eq = paras["t_eq"]
                T_start = paras["starting_T"]
                T_end = paras["end_T"]
                t = np.array(df["t"])
                beta = paras["beta"]
                # normalizing so the quench begins at zero:
                t = t - t_eq
                t_tau = t / tau
                T = np.zeros(len(t_tau))
                # calcing T
                xi_x = np.array(df["xi_x"][1:])
                xi_y = np.array(df["xi_y"][1:])
                lat_dim = int(np.sqrt(sys_df.iloc[0].size))

                fig, axes = plt.subplots(1, 2, figsize=(16, 9))
                corr_ax = axes[1]


                # Setze Tickmarken und Labels
                corr_ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
                corr_ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

                maxt_tau = np.maximum(int(np.max(t_tau)), 1)
                corr_ax.xaxis.set_major_locator(ticker.MultipleLocator(base=maxt_tau / 4))
                corr_ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=maxt_tau / 4 / 5))
                # TODO minor locator muss
                corr_ax.yaxis.set_minor_locator((plt.MultipleLocator(0.01)))
                # Füge Gitterlinien hinzu
                corr_ax.grid(which='major', linestyle='--', alpha=0.5)
                # second x axis for the temperature
                axes2 = corr_ax.twinx()
                # axes2.plot(t_tau, T, c="red", label="T", alpha=0.3)
                # Setze Tickmarken und Labels
                axes2.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
                axes2.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
                axes2.yaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
                axes2.yaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
                axes2.set_ylabel("T")

                corr_ax.set_ylabel(r"$\xi$")
                corr_ax.set_xlabel(r"t$/ \tau_Q$")
                corr_ax.set_title(rf"Quench protocol $\tau_Q = {int(tau)}$")
                save_plot(root + "plots/", str(tau) + ".png")
                ln, = corr_ax.plot([], [], ls="", marker=".")
                sys_ax = axes[0]
                x = np.arange(0.5, lat_dim + 1, 1)
                y = np.arange(0.5, lat_dim + 1, 1)
                z_values = np.array(sys_df.iloc[0][2:-1], dtype=float).reshape((lat_dim, lat_dim))
                vmin = -1.5 * np.sqrt(beta/2)
                vmax = 1.5 * np.sqrt(beta/2)
                pcm = sys_ax.pcolormesh(z_values, vmin=vmin, vmax=vmax)

                def animate_colormesh(pcm, df_row):
                    # takes one row and plots it
                    lat_dim = int(np.sqrt(df_row.size))
                    z_values = np.array(df_row, dtype=float).reshape((lat_dim, lat_dim))
                    pcm.set_data(z_values)
                    plt.draw()


                def update(frame, pcm, ln,  xi, sys_df, t):
                    t.append(frame)
                    i = len(t) % 127
                    df_row = sys_df.iloc[i][2:-1]
                    lat_dim = int(np.sqrt(df_row.size))
                    z_values = np.array(df_row, dtype=float).reshape((lat_dim, lat_dim))
                    pcm.set_array(z_values)
                    ln.set_data(t, xi[:i])
                    return ln, pcm


                corr_ax.set_ylim(np.min(xi_x) - 0.1 * np.min(xi_x), np.max(xi_x) + 0.1 * np.min(xi_x))
                corr_ax.set_xlim(np.min(t_tau), np.max(t_tau))
                for i, t_eq_tau in enumerate(t_tau[1:]):
                    if T_start - t_eq_tau < T_end:
                        T[i] = T_end
                    elif t_eq_tau > 0:
                        T[i] = T_start - t_eq_tau
                    else:
                        T[i] = T_start


                print(len(t_tau[1:-1]))
                ani = animation.FuncAnimation(fig, partial(update, pcm=pcm, ln=ln, xi=xi_x, sys_df=sys_df, t=[]), frames=t_tau[1:-1], interval=50, blit=True)
                plt.show()
                # corr_ax.scatter(t_tau[i], xi_x[i], c="C0")
                # corr_ax.scatter(t_tau[i], T[i], c="red")
                # save_plot(root + "Animation plots/", "frame" + str(i))
                # print(i, t_tau[i], xi_x[i])


                exit()
                corr_ax.plot(t_tau[1:], xi_x, ls="", marker=".", label=r"$\xi_x$", ms=2.5)
                corr_ax.plot(t_tau[1:], xi_y, ls="", marker=".", label=r"$\xi_y$", ms=2.5)

                xi = 1/2 * (xi_x + xi_y)
                xi_y +=  1/4 * np.max(xi_x)     # for plotting

                fig.legend()


                # save the stuff in the dics
                t_tau_dic[tau] = t_tau
                xi_avg_dic[tau] = xi

    fig, axes = plt.subplots(1, 1)
    # plot everything

    # Setze Tickmarken und Labels
    axes.tick_params(direction='in', which='both', length=6, width=2,
                     labelsize=9)
    axes.tick_params(direction='in', which='minor', length=3, width=1,
                     labelsize=9)
    axes.set_yscale("log")
    maxt_tau = 1
    max_xi = 0.5
    for tau in t_tau_dic.keys():
        axes.plot(t_tau_dic[tau][1:], xi_avg_dic[tau], ls="", marker=".", label=rf"{tau}", ms=2.5)
        maxt_tau = np.maximum(int(np.max(t_tau_dic[tau])), maxt_tau)
        axes.xaxis.set_major_locator(ticker.MultipleLocator(base=maxt_tau / 4))
        axes.xaxis.set_minor_locator(ticker.MultipleLocator(base=maxt_tau / 4 / 5))
        # TODO minor locator muss
        max_xi = round(np.maximum(np.max(xi_avg_dic[tau]), max_xi), 1)
        axes.yaxis.set_major_locator(ticker.MultipleLocator(base=max_xi / 4))
        axes.yaxis.set_minor_locator(ticker.MultipleLocator(base=max_xi / 4 /5 ))
        # Füge Gitterlinien hinzu
        axes.grid(which='major', linestyle='--', alpha=0.5)
        # second x axis for the temperature
        axes.set_ylabel(r"$\xi$")
        axes.set_xlabel(r"t$/ \tau_Q$")
        axes.set_title(rf"Quench protocol")
    axes.set_xlim(-0.5, 1.25)
    fig.legend()
    save_plot(root + "plots/", "together" + ".png")
    plt.show()

if __name__ == "__main__":
    main()