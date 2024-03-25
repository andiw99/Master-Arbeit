import functools

from FunctionsAndClasses import *
import matplotlib; matplotlib.use("TkAgg")
from functools import partial
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.integrate import dblquad
from matplotlib.ticker import MultipleLocator

def boltzmann_dist(x, p, beta, potential):
    return np.exp(- beta * ((p ** 2) / 2 + potential(x, p)))


def boltzmann_dist_x(x, beta, potential):
    return np.sqrt(2 * np.pi / beta) * np.exp(- beta * potential(x))

def int_boltzmann_dist_x(x, beta, potential, x2_range):
    # We want to integrate the boltzmann dist over the x2 range
    # to basically average over the cut point
    # the potential should be a function only in x2?
    # or we make it one with functools partial
    # so potential is now a function in x1, x2 and yields the total? potential
    # energy
    W = np.zeros(len(x))
    for i, x_point in enumerate(x):
        pot = functools.partial(potential, x2=x_point)
        # we now want to integrate e^-beta pot over the x2 range
        W_x2 = np.sqrt(2 * np.pi / beta) * np.exp(- beta * pot(x2_range))
        print(W_x2)
        W_avg = np.trapz(W_x2, x2_range) / (x2_range[-1] - x2_range[0])
        W[i] = W_avg
    print(W)
    return W

def boltzmann_dist_2d(x1, x2, beta, potential):
    return np.exp(- beta * potential(x1, x2))

def start_ani(event, ax, bars, x, x_cut, dx, nr_bins=100):
    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        for line in x:
            [b.remove() for b in bars]
            line_filter = filter_cut(line[2:], x_cut, dx)
            count, bins, bars = ax.hist(line_filter, nr_bins, density=True, color=colors[1])
            plt.draw()
            ax.set_title(f"t = {line[0]}")
            plt.pause(0.1)


def cos_pair_potential(x1, x2, h, J):
    return h * np.cos(2 * x1) - J * np.cos(x1 - x2)


def cos_combined_pair_potential(x1, x2, h, J, p):
    return h * np.cos(p * x1) + h * np.cos(p * x2) - J * np.cos(2 * (x1 - x2))

def cos_potential_x(x, p=2, h=20):
    return  h * np.cos(p * x)

def filter_cut(x, x_cut, dx):
    # so in x there are always (x1, x2) pairs and we now want a cut
    # of W(x1) around a specific x2
    # Therefore we filter only the pairs (x1, x2) with x2 in x2 - dx/2, x2 + dx/2
    x1 = x[0::2]        # all left lattice sites
    x2 = x[1::2]        # all right lattice sites
    x1 = x1[(x2 < x_cut + dx) & (x2 > x_cut - dx)]

    return x1


def main():
    root = "../../Generated content/Silicon/Benchmarks/Pairs/1e-5/longest/2"
    #root = "../../Generated content/Final/Benchmarks/2e-2/eta=1/h=100-J=40/T=50/longer/2"
    root = "../../Generated content/Final/Benchmarks/1e-5/Final/2"
    root_dirs = list_directory_names(root)
    file_extension = ".csv"
    potential = cos_potential_x
    potential_total = cos_combined_pair_potential
    # x_range = np.linspace(0, 2 * np.pi, 200)
    x_range = np.linspace(-np.pi / 2, np.pi / 2, 200)
    T = 20
    x = []
    dx = 0.15    # we look at a cut of the distribution that is 0.1 thick
    x2 = np.pi / 3   # the cut is located around x2
    nr_bins = 100
    x2_range = np.linspace(x2 - dx, x2 + dx, 100)


    for temp_dir in root_dirs:
        if (temp_dir != "NER plots") and (temp_dir != "plots"):
            temp_path = os.path.join(root, temp_dir)

            # now we need all ner files in this folder
            files = os.listdir(temp_path)
            for file in files:
                if os.path.splitext(file)[1] == file_extension:
                    file_path = os.path.join(temp_path, file)
                    para_filepath = os.path.splitext(file_path)[0] + ".txt"
                    parameters = read_parameters_txt(para_filepath)
                    with open(file_path, "rb") as f:
                        nr_rows = sum(1 for _ in f)
                    row = range(nr_rows)
                    x = read_large_df_array(file_path, row)
                    print(len(x))
                    T = parameters["T"]
                    dt = parameters["dt"]
                    alpha = parameters["alpha"]
                    J = parameters["J"]
                    p = parameters["p"]
                    break

        # first lets plot it in 3d to get a feeling for the distribution
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X1, X2 = np.meshgrid(x_range, x_range)
        pot_2d = functools.partial(potential_total, h=alpha, J=J, p=p)
        pot_vals = pot_2d(X1, X2)
        W = np.exp(-1 / T * pot_vals)
        boltz_2d = functools.partial(boltzmann_dist_2d, beta=1 / T, potential=pot_2d)
        Z = dblquad(boltz_2d, 0, 2 * np.pi, 0, 2 * np.pi)
        W /= Z[0]

        surf = ax.plot_surface(X1, X2, W, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        # Customize the z axis.
        ax.zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax.zaxis.set_major_formatter('{x:.02f}')

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()

        beta = 1 / T
        pot_single = functools.partial(potential, h=alpha, p=p)
        pot_total = functools.partial(potential_total, h=alpha, J=J, p=p)
        W_x = int_boltzmann_dist_x(x_range, beta, pot_total, x2_range)
        W_x_single = boltzmann_dist_x(x_range, beta, pot_single)
        Z = np.trapz(W_x, x_range)
        Z_single = np.trapz(W_x_single, x_range)
        W_x /= Z
        W_x_single /= Z_single
        fig, ax = plt.subplots(1, 1, figsize=(10, 10 * 4.8/6.4))
        x_start = filter_cut(x[0][2:], x2, dx)
        count, bins, bars = ax.hist(x_start, nr_bins, density=True, label=f"d$\sigma= {dt}$", color=colors[1])
        #ax.plot(x_range, W_x, label=f"T = {T:.2f}\n"+rf"$\vartheta_2 = ${x2:.2f}", color=colors[0], linewidth=3)
        ax.plot(x_range, W_x, label=r"$p(\vartheta_1) \propto $exp$({-\beta H(\vartheta_1, \vartheta_2)})$" + "\n" + rf"$\vartheta_2 = ${x2:.2f}",
                color=colors[4], linewidth=3)
        #ax.plot(x_range, W_x_single, label=f"single particle dist")
        ax.set_title(f"t = {x[0][0]}, dt = {dt}")
        ax.set_ylim(1e-7, 1.2 * np.max(W_x))
        ax.set_xlim(np.min(x_range), np.max(x_range))
        ax.set_ylabel(r"$p(\vartheta_1)$")
        ax.set_xlabel(r"$\vartheta_1$")

        #ax.set_yscale("log")
        config = {
            "increasefontsize": 1.5,
            "labelhorizontalalignment": "right",
        }
        configure_ax(fig, ax, config)
        ax.xaxis.set_major_locator(MultipleLocator(base=np.pi / 8))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(pi_formatter))
        ax.xaxis.set_minor_locator(MultipleLocator(base=np.pi / (5 * 4)))
        animation = partial(start_ani, ax=ax, bars=bars, x=x, x_cut=x2, dx=dx, nr_bins=nr_bins)
        fig.canvas.mpl_connect("button_press_event", animation)
        create_directory_if_not_exists(root + "/plots/")
        plt.show()
        ax.set_title(f"")
        fig.savefig(root + f"/plots/p-dt-{dt}-theta2-{x2:.2f}.png", dpi=300)





if __name__ == "__main__":
    main()