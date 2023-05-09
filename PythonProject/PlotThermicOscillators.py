from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np

def plot_trajectories(df, parameters):
    times = df.iloc[:, 0]
    x_values = df.iloc[:, 1:-1]
    n = x_values.shape[1]
    x_avg = x_values.mean(axis=1)

    # extract values for theoretical curve
    eta = parameters["eta"]
    alpha = parameters["alpha"]
    N = times.size
    a = times.iloc[0]
    b = times.iloc[-1]
    x0 = x_avg[0]

    print(a, "  ", b)

    t, osc = theoretical_trajectory(eta, alpha, x0, a, b, N)

    T = parameters["T"]
    dt = parameters["dt"]
    fig, axes = plt.subplots(1, 1)
    axes.plot(times, x_avg, label=f" {n} Averaged Oscialltors in Bath ")
    axes.plot(t, osc, linestyle="dashed", label="Theoretical Trajectory")
    axes.set_xlabel("t")
    axes.set_ylabel("x(t)")
    axes.set_ylim((-x0, 1.1 * x0))
    axes.set_xlim((a, b))
    axes.set_title(f"Uncoupled harmonic oscillators in T={T} \n dt = {dt}")
    axes.legend()

    plt.show()


def theo_msd(eta, alpha, T, a, b, N):
    t = np.linspace(a, b, N)
    return t, T / alpha * (1 - np.exp(-2 * alpha / eta * t))
def plot_theo_msd(df, parameters, savepath):
    eta = parameters["eta"]
    alpha = parameters["alpha"]
    T = parameters["T"]
    times = df.iloc[:, 0]
    N = times.size
    a = times.iloc[0]
    b = times.iloc[-1]

    x_values = df.iloc[:, 1:-1]
    n = x_values.shape[1]


    t, msd_theo = theo_msd(eta, alpha, T, a, b, N)
    msd_avg = x_values ** 2
    msd_avg = msd_avg.mean(axis=1)
    dt = parameters["dt"]
    fig, axes = plt.subplots(1, 1)
    axes.plot(times, msd_avg, label=f" {n} Oscialltors")
    axes.plot(t, msd_theo, linestyle="dashed", label="Theoretical MSD")
    axes.set_xlabel("t")
    axes.set_ylabel("MSD(t)")
    axes.set_ylim((0.9 * np.min(msd_avg), 1.1 * np.max(msd_avg)))
    axes.set_xlim((a, b))
    axes.set_title(f"Uncoupled harmonic oscillators in T={T} \n dt = {dt}")
    axes.legend()
    name = str(dt + T)
    save_plot(savepath, name)
    plt.show()


def theoretical_trajectory(eta, alpha, x0, a, b, N=200):
    xi = eta/ (2 * np.sqrt(alpha))
    w0 = np.sqrt(alpha)
    w = w0 * np.sqrt(1 - xi ** 2)
    y = w * xi / np.sqrt(1 - xi ** 2)

    # stützstellen
    t = np.linspace(a, b, N)
    osc = x0 * np.exp(-y * t) * np.cos(w * t)

    return t, osc


def main():
    root = "../../../Generated content/Convergence Check MSD/"
    # /home/andi/Documents/Master-Arbeit Code/Generated content/GPU Oscillators/eta=0.20/T=500.00/dt=0.0010
    filepaths = new_files_in_dir(root, root, plot_all=False)
    print(filepaths)
    for filepath in filepaths:
        df = read_csv(filepath)
        txt_filepath = os.path.splitext(filepath)[0] + ".txt"
        parameters = read_parameters_txt(txt_filepath)

        # average and plot, actually no big deal
        plot_theo_msd(df, parameters, root + "plots/")
        #plot_trajectories(df, parameters)

    plotted_files = open(root + "/plotted_files.txt", "a")
    for f in filepaths:
        plotted_files.write(f + "\n")
    # newline for visibility
    plotted_files.write("\n")


if __name__ == "__main__":
    main()