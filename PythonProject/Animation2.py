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
        self.cell_L = 128
        self.lines_at_a_time = 100

    def update(self, frame):
        # index has to be modulo sinze we animate multipel taus
        ind = frame % ((self.sys_values))
        print(frame, ind, self.sys_values)

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
                self.ani.event_source.stop()
                exit()
        # updating system
        if ind % self.lines_at_a_time == 0:
            self.ts, self.Ts, self.sys_rows = self.read_multiple_lines(ind, self.lines_at_a_time)
        self.t_tau.append((self.ts[ind % self.lines_at_a_time]-self.t_eq) / self.tau)
        self.T.append(self.Ts[ind % self.lines_at_a_time])
        z_values = np.array(self.sys_rows[ind % self.lines_at_a_time], dtype=float).reshape((self.cell_L, self.cell_L))
        if self.J < 0:
            z_values = chess_board_trafo(z_values)
        self.pcm.set_array(z_values)
        # updating corr length

        # only every like 10th run
        xi_ind = int(ind) + 1
        print(xi_ind)
        self.ln.set_data(self.t_tau, self.xi[:xi_ind])
        self.ln_temp.set_data(self.t_tau, self.T)
        lower_bound = np.minimum(np.min(self.xi[:xi_ind]) - 0.1 * np.min(self.xi[:xi_ind]), self.corr_ax.get_ylim()[0])
        upper_bound = np.maximum(np.max(self.xi[:xi_ind]) + 0.1 * np.min(self.xi[:xi_ind]), self.corr_ax.get_ylim()[1])
        self.corr_ax.set_ylim(lower_bound, upper_bound)

        return self.ln, self.pcm


    def reshape(self, df_row):
        lat_dim = int(np.sqrt(df_row.size))
        return np.array(df_row, dtype=float).reshape((lat_dim, lat_dim))


    def next_quench(self):
        # assuming t_tau stays the same
        v = np.sqrt(self.beta/2)
        # plot the first value of the time and the correlation length
        print("Before plotting first xi")
        self.ln,  = self.corr_ax.plot([], [], ls="",
                                      marker=".", label=rf"$\tau = {self.tau}$", c=colors[self.tau_nr - 1])
        print("Before plotting first temp")
        self.ln_temp, = self.twin_corr_ax.plot([], [], alpha=0.5, c="red")
        #self.twin_corr_ax.set_xlim(np.min(self.t_tau), np.max(self.t_tau))
        self.twin_corr_ax.set_ylim(self.paras["end_T"] - 0.05, self.paras["starting_T"] + 0.05)
        # plot first system frame
        # Frames of the animation are all ts excluded the first one
        # but we just use the count
        # plot the first value of the time and the correlation length
        # plot first system frame
        # pcm will be overwritten so we use set_array
        self.corr_ax.legend(bbox_to_anchor=(0.07, 0.3, 0.2, 0.5))
        if self.init:
            self.ani = animation.FuncAnimation(self.fig, self.update, repeat=False, interval=50, init_func=self.init_plots)
            self.pcm = self.sys_ax.pcolormesh(self.reshape(np.zeros(self.cell_L ** 2)), vmin=-1.5 * v, vmax=1.5 * v)
            print("Before starting the animation")
            ticks = np.linspace(- 1.5 *  v, 1.5 * v, 7, endpoint=True)
            tick_labels = np.linspace(-1.5, 1.5, 7, endpoint=True)
            tick_labels = [str(tick_label) for tick_label in tick_labels]
            self.cbar = self.fig.colorbar(self.pcm, ticks=ticks)
            self.cbar.ax.set_yticklabels(tick_labels)
            plt.tight_layout()
            self.init = False
        else:
            self.pcm.set_array(self.reshape(np.zeros(self.cell_L ** 2)))
            #self.ani = animation.FuncAnimation(self.fig, self.update, repeat=False, interval=10)

    def init_plots(self):
        tau_span = (self.paras["end_t"]) / self.tau
        print("tau_span:", tau_span)
        self.corr_ax.set_xlim(- self.t_eq/ self.tau - 0.05 * tau_span, tau_span - self.t_eq / self.tau + 0.05 * tau_span)
        self.corr_ax.set_ylim(np.min(self.xi) - 0.1 * np.min(self.xi), np.max(self.xi) + 0.1 * np.min(self.xi))


    def read_line(self, ind):
        """
        returns t, T and the system
        :param ind:
        :return:
        """
        with open(self.dir_path + "/0.csv") as f:
            for i, line in enumerate(f):
                if i == ind:
                    data = string_to_array(line[:-2])
                    t = data[0]
                    T = data[1]
                    sys_row = extract_cell(data[2:], 0, self.cell_L)
                    return t, T, sys_row

    def read_multiple_lines(self, ind, n):
        """
        returns t, T and the system values for the next n lines
        :param ind:
        :return:
        """
        with open(self.dir_path + "/0.csv") as f:
            ts = []
            Ts = []
            sys_rows = []
            for i, line in enumerate(f):
                if (i > ind):
                    if (i > ind + n):
                        return ts, Ts, sys_rows
                    data = string_to_array(line[:-2])
                    t = data[0]
                    T = data[1]
                    sys_row = extract_cell(data[2:], 0, self.cell_L)
                    ts.append(t)
                    Ts.append(T)
                    sys_rows.append(sys_row)
            return ts, Ts, sys_rows


    def get_tau_dirs(self):
        root_dirs = os.listdir(self.root)
        tau_vals = []
        for dir in root_dirs:
            print(dir)
            print(os.path.isdir(self.root + dir))
            if(dir != "plots") & (os.path.isdir(self.root + dir)):
                self.tau_dirs.append(self.root + dir)
                tau_vals.append(float(dir))
        # sort for tau size
        self.tau_dirs = np.array(self.tau_dirs)[np.argsort(tau_vals)]


    def read_folder(self, tau_nr):
        self.dir_path = self.tau_dirs[tau_nr]
        dir_path = self.dir_path
        filename = dir_path + "/quench.process"
        df = read_struct_func(filename)
        paras = read_parameters_txt(dir_path + "/0.txt")
        self.paras = paras
        self.sys_values = int(paras["nr_save_values"])
        self.tau = paras["tau"]
        self.t_eq = paras["t_eq"]
        T_start = paras["starting_T"]
        T_end = paras["end_T"]
        self.beta = paras["beta"]
        self.J = paras["J"]
        # normalizing so the quench begins at zero:
        self.t_tau = []
        self.T = []
        # calcing T
        xi_x = np.array(df["xi_x"][1:])
        xi_y = np.array(df["xi_y"][1:])
        self.xi = 1/2 * (xi_x + xi_y)
        self.tau_nr += 1

    def safe(self, path="foo"):
        FFwriter = animation.FFMpegWriter(fps=40)
        self.ani.save('../../Generated content/LargeAnimationLegend.mp4', writer=FFwriter, dpi=300)


def main():
    root = "../../Generated content/Defense2/Large Quench Video/"
    name = "quench.process"
    png_name = "quench.png"
    root_dirs = os.listdir(root)
    print(root_dirs)
    xi_avg_dic = {}
    t_tau_dic = {}

    # The problem could be that we cannot read in the df if we save like 1000 values.
    # We might have to consider building something that splits up the file in 100 MB portions?

    fig, axes = plt.subplots(1, 2, figsize=(16, 9),  gridspec_kw={'width_ratios': [19, 10]})
    corr_ax = axes[1]

    # Setze Tickmarken und Labels
    corr_ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    corr_ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

    #maxt_tau = np.maximum(int(np.max(t_tau)), 1)
    #corr_ax.xaxis.set_major_locator(ticker.MultipleLocator(base=maxt_tau / 4))
    #corr_ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=maxt_tau / 4 / 5))
    # TODO minor locator muss
    #corr_ax.yaxis.set_minor_locator((plt.MultipleLocator(0.01)))
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
    # plt.show()
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