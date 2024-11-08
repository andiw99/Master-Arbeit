import functools

from FunctionsAndClasses import *
import matplotlib; matplotlib.use("TkAgg")
from functools import partial

def boltzmann_dist(x, p, beta, potential):
    return np.exp(- beta * ((p ** 2) / 2 + potential(x, p)))


def boltzmann_dist_x(x, beta, potential):
    return np.sqrt(2 * np.pi / beta) * np.exp(- beta * potential(x))


def cos_potential_x(x, p_multi=2, alpha=20):
    return alpha * np.cos(p_multi * x)


def harmonic_potential(x, alpha=20):
    return alpha * x**2

def start_ani(event, ax, bars, x, nr_bins=100):
    mode = event.canvas.toolbar.mode
    if event.button == 1 and event.inaxes == ax and mode == "":
        for line in x:
            [b.remove() for b in bars]
            count, bins, bars = ax.hist(line[2:], nr_bins, density=True, color="C0")
            plt.draw()
            ax.set_title(f"t = {line[0]}")
            plt.pause(0.1)



def main():
    root = "../../Generated content/BBK/Dist Convergence/BBK/0.05"
    #root = "../../Generated content/Testing Convergence/0.01/"
    root_dirs = list_directory_names(root)
    file_extension = ".csv"
    potential = cos_potential_x
    x_range = np.linspace(0, 2 * np.pi, 200)
    T = 1
    x = []
    nr_bins = 100

    for temp_dir in root_dirs[::-1]:
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
                    break

        beta = 1 / T
        pot = functools.partial(potential, alpha=alpha)
        W_x = boltzmann_dist_x(x_range, beta, pot)
        Z = np.trapz(W_x, x_range)
        W_x /= Z
        fig, ax = plt.subplots(1, 1)
        count, bins, bars = ax.hist(x[0][2:], nr_bins, density=True, label=f"dt = {dt}")
        ax.plot(x_range, W_x, label=f"T = {T:.2f}")
        ax.set_title(f"t = {x[0][0]}, dt = {dt}")
        ax.set_ylim(0, 1.5 * np.max(W_x))
        configure_ax(fig, ax)
        animation = partial(start_ani, ax=ax, bars=bars, x=x, nr_bins=nr_bins)
        fig.canvas.mpl_connect("button_press_event", animation)
        plt.show()





if __name__ == "__main__":
    main()