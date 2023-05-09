from FunctionsAndClasses import *
from matplotlib import pyplot as plt

def main():
    root = "../../Generated content/Convergence Check Diffusion more time/"
    name = "results.csv"

    df = read_csv(root + name)

    print(df)
    # listen anlegen
    print(df.iloc[:, 1])
    dt = df.index
    mu_em = df.iloc[:, 0]
    msd_em = df.iloc[:, 1]
    err_mu_em = df.iloc[:, 2]
    err_msd_em = df.iloc[:, 3]

    mu_lm  = df.iloc[:, 4]
    msd_lm  = df.iloc[:, 5]
    err_mu_lm  = df.iloc[:, 6]
    err_msd_lm  = df.iloc[:, 7]

    fig, axes = plt.subplots(2, 1)
    # plot versus dt
    axes[0].set_title("Distribution parameters of simple Diffusion \n Werte ermittelt aus 40000 Realisierungen \n Mean Value")
    axes[0].plot(dt, err_mu_em, label = "Error Mittelwert")
    axes[0].set_xlabel("dt")
    axes[0].set_ylabel("abs error")
    axes[1].plot(dt, err_msd_em, label="Error MSD")
    axes[0].plot(dt, err_mu_lm, label="Error Mittelwert LM")
    axes[1].plot(dt, err_msd_lm, label="Error MSD LM")
    axes[1].set_xlabel("dt")
    axes[1].set_ylabel("abs error")
    axes[1].set_title("MSD")
    axes[0].legend()
    axes[1].legend()
    axes[0].set_xscale("log")
    axes[1].set_xscale("log")
    plt.tight_layout()

    # plot vs 1/dt
    fig, axes = plt.subplots(2, 1)
    steps = 1/np.array(dt)
    axes[0].set_title(
        "Distribution parameters of simple Diffusion \n Werte ermittelt aus 40000 Realisierungen \n Mean Value")
    axes[0].plot(steps, err_mu_em, label="Error Mittelwert")
    axes[0].set_xlabel("steps")
    axes[0].set_ylabel("abs error")
    axes[1].plot(steps, err_msd_em, label="Error MSD")
    axes[0].plot(steps, err_mu_lm, label="Error Mittelwert LM")
    axes[1].plot(steps, err_msd_lm, label="Error MSD LM")
    axes[1].set_xlabel("steps")
    axes[1].set_ylabel("abs error")
    axes[1].set_title("MSD")
    axes[0].legend()
    axes[1].legend()

    axes[0].set_xscale("log")
    axes[1].set_xscale("log")

    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()