import functools

from FunctionsAndClasses import *
import matplotlib; matplotlib.use("TkAgg")
from functools import partial
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.integrate import dblquad

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
    for i,x_point in enumerate(x):
        pot = functools.partial(potential, x2=x_point)
        # we now want to integrate e^-beta pot over the x2 range
        W_x2 = np.exp(- beta * pot(x2_range))
        print(W_x2)
        exit()
        W_avg = np.trapz(W_x2, x2_range) / (x2_range[-1] - x2_range[0])
        W[i] = W_avg
    return W

def boltzmann_dist_2d(x1, x2, beta, potential):
    return np.exp(- beta * potential(x1, x2))

def start_ani(event, ax, bars, x, x_cut, dx, nr_bins=100):
    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        for line in x:
            [b.remove() for b in bars]
            line_filter = filter_cut(line[2:], x_cut, dx)
            count, bins, bars = ax.hist(line_filter, nr_bins, density=True, color="C0")
            plt.draw()
            ax.set_title(f"t = {line[0]}")
            plt.pause(0.1)



def cos_pair_potential(x1, x2, h, J):
    return h * np.cos(2 * x1) - J * np.cos(x1 - x2)


def cos_combined_pair_potential(x1, x2, h, J):
    return h * np.cos(2 * x1)  + h * np.cos(2 * x2) - J * np.cos(x1 - x2)



def filter_cut(x, x_cut, dx):
    # so in x there are always (x1, x2) pairs and we now want a cut
    # of W(x1) around a specific x2
    # Therefore we filter only the pairs (x1, x2) with x2 in x2 - dx/2, x2 + dx/2
    x1 = x[0::2]        # all left lattice sites
    x2 = x[1::2]        # all right lattice sites
    print(len(x1))
    x1 = x1[(x2 < x_cut + dx) & (x2 > x_cut - dx)]

    print(len(x1), x1)
    return x1


def main():
    root = "../../Generated content/XY Pairs/Test/"
    #root = "../../Generated content/Testing Convergence/0.01/"
    root_dirs = list_directory_names(root)
    file_extension = ".csv"
    potential = cos_pair_potential
    potential_total = cos_combined_pair_potential
    x_range = np.linspace(0, 2 * np.pi, 200)
    T = 20
    x = []
    dx = 0.1        # we look at a cut of the distribution that is 0.1 thick
    x2 = np.pi/2    # the cut is located around x2
    nr_bins = 100
    x2_range = np.linspace(x2 - dx, x2 + dx, 10)

    # first lets plot it in 3d to get a feeling for the distribution
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X1, X2 = np.meshgrid(x_range, x_range)
    pot_2d = functools.partial(potential_total, h=20, J=20)
    pot_vals = pot_2d(X1, X2)
    W = np.exp(-1/T * pot_vals)
    print(pot_2d(np.pi/2, np.pi/2))
    boltz_2d = functools.partial(boltzmann_dist_2d, beta=1/T, potential=pot_2d)
    print(boltz_2d(np.pi/2, np.pi/2))
    Z = dblquad(boltz_2d, 0, 2 * np.pi, 0, 2*np.pi)
    print(Z)
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
                    row = range(int(parameters["nr_saves"]))
                    x = read_large_df(file_path, row)[1:]
                    T = parameters["T"]
                    dt = parameters["dt"]
                    alpha = parameters["alpha"]
                    J = parameters["J"]
                    break

        beta = 1 / T
        pot = functools.partial(potential, x2=x2, h=alpha, J=J)
        pot_total = functools.partial(potential_total, h=alpha, J=J)
        W_x = int_boltzmann_dist_x(x_range, beta, pot, x2_range)
        Z = np.trapz(W_x, x_range)
        W_x /= Z
        fig, ax = plt.subplots(1, 1)
        x_start = filter_cut(x[0][2:], x2, dx)
        count, bins, bars = ax.hist(x_start, nr_bins, density=True, label=f"dt = {dt}")
        ax.plot(x_range, W_x, label=f"T = {T:.2f}")
        ax.set_title(f"t = {x[0][0]}, dt = {dt}")
        ax.set_ylim(0, 1.5 * np.max(W_x))
        configure_ax(fig, ax)
        animation = partial(start_ani, ax=ax, bars=bars, x=x, x_cut=x2, dx=dx, nr_bins=nr_bins)
        fig.canvas.mpl_connect("button_press_event", animation)
        plt.show()





if __name__ == "__main__":
    main()