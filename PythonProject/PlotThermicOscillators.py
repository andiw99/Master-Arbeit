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

def theo_sigma_xx(eta, alpha, T, a, b, N = 100):
    lambda_1 = 0.5 * (eta + np.sqrt(eta ** 2 - 8 * alpha + 0j))
    lambda_2 = 0.5 * (eta - np.sqrt(eta ** 2 - 8 * alpha + 0j))
    t = np.linspace(a, b, N)
    pref = eta * T / ((lambda_1 - lambda_2) ** 2)
    sum_1 = (lambda_1 + lambda_2) / (lambda_1 * lambda_2)
    sum_2 = 4 / (lambda_1 + lambda_2) * (np.exp(- (lambda_1 + lambda_2) * t) - 1)
    sum_3 = 1 / lambda_1 * np.exp(-2 * lambda_1 * t)
    sum_4 = 1 / lambda_2 * np.exp(-2 * lambda_2 * t)
    return t, np.real(pref * (sum_1 + sum_2 - sum_3 - sum_4))

def plot_theo_msd(fig, axes, df, parameters, savepath, title="", color=colors[5]):
    eta = parameters["eta"]
    alpha = parameters["alpha"]
    T = parameters["T"]
    dt = parameters["dt"]
    eta = parameters["eta"]

    times = df.iloc[:, 0]
    N = times.size
    a = times.iloc[0]
    b = times.iloc[-1]

    x_values = df.iloc[:, 2:-1]
    n = x_values.shape[1]
    print(x_values)

    t, msd_theo = theo_sigma_xx(eta, alpha, T, a, b, N)
    msd_theo = msd_theo
    msd_avg = x_values ** 2
    msd_avg = msd_avg.mean(axis=1)
    times = np.array(times)
    msd_avg = np.array(msd_avg)
    dt = parameters["dt"]

    axes.plot(t, msd_theo, label=rf"$\eta = {eta:.1f}$", color=color)
    axes.plot(times, msd_avg, linestyle="dotted", color=colors[0])
    axes.set_xlabel("t")
    axes.set_ylabel(r"$\left\langle x^2(t)\right \rangle$")
    axes.set_ylim((0.9 * np.min(msd_avg), 1.05 * np.max(msd_avg)))
    axes.set_xlim((a, b))
    #axes.set_title(f"Uncoupled harmonic oscillators in T={T} \n dt = {dt}")
    axes.legend()



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
    root = "../../Generated content/BBK/MSD/BBK/0.025"
    root = "../../Generated content/Final/Benchmarks/HO/Euler/T=20/dt = 0.0005/eta comp"
    title = "Euler method"
    # /home/andi/Documents/Master-Arbeit Code/Generated content/GPU Oscillators/eta=0.20/T=500.00/dt=0.0010

    a = 0
    b = 3
    T = 20
    eta = 0.5
    alpha = 20

    fig, ax = plt.subplots(1, 1)
    t, sigma = theo_sigma_xx(eta, alpha, T, a, b, 200)
    t_msd, msd = theo_msd(eta, alpha, T, a, b, 200)
    ax.plot(t, sigma)
    # ax.plot(t_msd, msd)
    configure_ax(fig, ax)
    plt.show()
    filepaths = list_directory_names(root)
    print(filepaths)
    for fpath in filepaths:
        temppath = os.path.join(root, fpath)
        files = os.listdir(temppath)
        if fpath != "plots":
            fig, axes = plt.subplots(1, 1)
            ind = 0
            for i, file in enumerate(files):
                if os.path.splitext(file)[1] == ".csv":

                    filepath = os.path.join(temppath, file)
                    df = read_csv(filepath)
                    print(df)
                    txt_filepath = os.path.splitext(filepath)[0] + ".txt"
                    parameters = read_parameters_txt(txt_filepath)
                    dt = parameters["dt"]

                    # average and plot, actually no big deal
                    plot_theo_msd(fig, axes, df, parameters, os.path.join(root, "plots/"), title=title, color=colors[3 + 2 * ind])
                    ind += 1

            axes.plot([], [], label=f"dt = {dt}", linestyle="dotted", color=colors[0])

            name = f"{T}-{dt}-{alpha}-{eta}.png"

            config = {
                "increasefontsize": 0.75,
                "labelhorizontalalignment": "center",
                "labelverticalalignment": "bottom",
                "labelrotation": 90,
                "legendtitle": title,
            }
            savepath =  os.path.join(root, "plots/")
            configure_ax(fig, axes, config)
            create_directory_if_not_exists(savepath)
            plt.savefig(os.path.join(savepath, name), format="png")
            plt.show()

        #plot_trajectories(df, parameters)

    plotted_files = open(root + "/plotted_files.txt", "a")
    for f in filepaths:
        plotted_files.write(f + "\n")
    # newline for visibility
    plotted_files.write("\n")


if __name__ == "__main__":
    main()