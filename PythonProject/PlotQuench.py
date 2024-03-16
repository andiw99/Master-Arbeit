from Suite import quench_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    simpath = "../../Generated content/Silicon/Subsystems/Suite/h/1.7320508075688776/Quench/"

    taus = [512, 2048, 8192]
    xi_ampl = 1.2
    Tc = 1.2
    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc)
    plt.show()

    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc)
    new_xlim = (0.3, 0.5)
    axes[1].set_xlim(new_xlim)
    axes[1].set_ylim(5, 20)
    configure_ax(fig, axes[1])
    #ax = zoom_plot(axes[1], new_xlim)
    ax = axes[1]
    ax.set_xlim(new_xlim)
    plt.show()

if __name__ == "__main__":
    main()