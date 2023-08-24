import numpy as np
from FunctionsAndClasses import *
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as animation
from functools import partial
import matplotlib; matplotlib.use("TkAgg")


class Anim():
    def __init__(self, fig, axes, root):
        self.fig = fig
        self.sys_ax = axes[0]
        self.corr_ax = axes[1]
        self.twin_corr_ax = axes[2]
        self.init = True
        self.root = root
        self.tau_dirs = []
        self.get_tau_dirs()
        self.tau_nr = 0         # nr of current tau
        # read first df
        self.sys_values = 0
        self.sys_per_corr = 0
        self.read_folder(0)

    def update(self, frame):
        # index has to be modulo sinze we animate multipel taus
        ind = frame % ((self.sys_values))
        print(frame, ind)
        if (ind == (self.sys_values) -1):
            #self.ani.pause()
            if self.tau_nr < len(self.tau_dirs):
                print("what is even this tau nr? ", self.tau_nr)
                self.read_folder(tau_nr=self.tau_nr)
                self.next_quench()
                #self.ani.resume()
                return self.ln, self.pcm
            else:
                print("Exhausting all tau folders")
                self.ani.pause()
                return
        # updating system

        df_row = self.sys_df.iloc[ind]
        lat_dim = int(np.sqrt(df_row.size))
        z_values = np.array(df_row, dtype=float).reshape((lat_dim, lat_dim))
        print("Here?")
        self.pcm.set_array(z_values)
        # updating corr length
        if (ind % self.sys_per_corr == 0):
            # only every like 10th run
            xi_ind = int(ind / self.sys_per_corr) + 1
            self.ln.set_data(self.t_tau[:xi_ind], self.xi[:xi_ind])
            self.ln_temp.set_data(self.t_tau[:xi_ind], self.T[:xi_ind])
            lower_bound = np.minimum(np.min(self.xi[:xi_ind]) - 0.1 * np.min(self.xi[:xi_ind]), self.corr_ax.get_ylim()[0])
            upper_bound = np.maximum(np.max(self.xi[:xi_ind]) + 0.1 * np.min(self.xi[:xi_ind]), self.corr_ax.get_ylim()[1])
            self.corr_ax.set_ylim(lower_bound, upper_bound)

        print("why?")
        return self.ln, self.pcm


    def reshape(self, df_row):
        lat_dim = int(np.sqrt(df_row.size))
        return np.array(df_row, dtype=float).reshape((lat_dim, lat_dim))


    def next_quench(self):
        # assuming t_tau stays the same
        v = np.sqrt(self.beta/2)
        # plot the first value of the time and the correlation length
        print("Before plotting first xi")
        self.ln,  = self.corr_ax.plot(self.t_tau[0], self.xi[0], ls="",
                                      marker=".", label=rf"$\tau = {self.tau}$", c=colors[self.tau_nr])
        print("Before plotting first temp")
        self.ln_temp, = self.twin_corr_ax.plot(self.t_tau[0], self.T[0], alpha=0.5, c="red")
        #self.twin_corr_ax.set_xlim(np.min(self.t_tau), np.max(self.t_tau))
        self.twin_corr_ax.set_ylim(np.min(self.T), np.max(self.T))
        # plot first system frame
        # Frames of the animation are all ts excluded the first one
        # but we just use the count
        # plot the first value of the time and the correlation length
        # plot first system frame
        # pcm will be overwritten so we use set_array
        self.corr_ax.legend()
        if self.init:
            self.init_plots()
            self.pcm = self.sys_ax.pcolormesh(self.reshape(self.sys_df.iloc[0]), vmin=-v, vmax=v)
            self.ani = animation.FuncAnimation(self.fig, self.update, repeat=False, interval=50)
            self.init = False
        else:
            self.pcm.set_array(self.reshape(self.sys_df.iloc[0]))
            #self.ani = animation.FuncAnimation(self.fig, self.update, repeat=False, interval=10)

    def init_plots(self):
        tau_span = np.max(self.t_tau[:-2]) - np.min(self.t_tau[2:])
        self.corr_ax.set_xlim(np.min(self.t_tau[2:] - 0.05 * tau_span), np.max(self.t_tau[:-2]) + 0.05 * tau_span)
        self.corr_ax.set_ylim(np.min(self.xi) - 0.1 * np.min(self.xi), np.max(self.xi) + 0.1 * np.min(self.xi))

    def get_tau_dirs(self):
        root_dirs = os.listdir(self.root)
        tau_vals = []
        for dir in root_dirs:
            if(dir != "plots") & (os.path.isdir(self.root + dir)):
                self.tau_dirs.append(self.root + dir)
                tau_vals.append(float(dir))
        # sort for tau size
        self.tau_dirs = np.array(self.tau_dirs)[np.argsort(tau_vals)]


    def read_folder(self, tau_nr):
        dir_path = self.tau_dirs[tau_nr]
        filename = dir_path + "/quench.process"
        sys_df = read_csv(dir_path + "/0.csv")
        df = read_struct_func(filename)
        paras = read_parameters_txt(dir_path + "/0.txt")
        self.sys_df = sys_df.iloc[:, 2:-1]
        self.sys_values = sys_df.shape[0]
        self.tau = paras["tau"]
        t_eq = paras["t_eq"]
        T_start = paras["starting_T"]
        T_end = paras["end_T"]
        t = np.array(df["t"])
        self.beta = paras["beta"]
        # normalizing so the quench begins at zero:
        t = t - t_eq
        self.t_tau = t / self.tau
        self.T = np.zeros(len(self.t_tau))
        # calcing T
        xi_x = np.array(df["xi_x"][1:])
        xi_y = np.array(df["xi_y"][1:])
        for i, t_eq_tau in enumerate(self.t_tau[1:]):
            if T_start - t_eq_tau < T_end:
                self.T[i] = T_end
            elif t_eq_tau > 0:
                self.T[i] = T_start - t_eq_tau
            else:
                self.T[i] = T_start
        self.sys_per_corr = int(self.sys_values / len(self.t_tau))
        self.xi = xi_x
        self.tau_nr += 1

    def safe(self, path="foo"):
        FFwriter = animation.FFMpegWriter(fps=10)
        self.ani.save('../../Generated content/animation.mp4', writer=FFwriter)


def main():
    root = "../../Generated content/Overdamped Quenching/"
    name = "quench.process"
    png_name = "quench.png"
    root_dirs = os.listdir(root)
    print(root_dirs)
    xi_avg_dic = {}
    t_tau_dic = {}

    # The problem could be that we cannot read in the df if we save like 1000 values.
    # We might have to consider building something that splits up the file in 100 MB portions?

    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    corr_ax = axes[1]

    # Setze Tickmarken und Labels
    corr_ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    corr_ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

    #maxt_tau = np.maximum(int(np.max(t_tau)), 1)
    #corr_ax.xaxis.set_major_locator(ticker.MultipleLocator(base=maxt_tau / 4))
    #corr_ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=maxt_tau / 4 / 5))
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
    #corr_ax.set_title(rf"Quench protocol $\tau_Q = {int(tau)}$")
    axes = np.append(axes, axes2)
    ani_class = Anim(fig, axes, root)
    ani_class.next_quench()
    plt.title(r"Quench Process for different Values of $\tau$")
    ani_class.safe()
    exit()
    for item in root_dirs:
        if (item != "plots"):
            dir_path = os.path.join(root, item)
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                df = read_struct_func(filename)
                sys_df = read_csv(dir_path + "/0.csv")
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






                sys_df_trunc = sys_df.iloc[:, 2:-1]

                ani_class.next_quench(sys_df_trunc, xi_x, t_tau)
                # corr_ax.scatter(t_tau[i], xi_x[i], c="C0")
                # corr_ax.scatter(t_tau[i], T[i], c="red")
                # save_plot(root + "Animation plots/", "frame" + str(i))
                # print(i, t_tau[i], xi_x[i])
    plt.show()

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