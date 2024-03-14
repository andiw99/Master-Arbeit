from Suite import amplitude_measurement

def main():
    simpath = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=1/0.4161791450287818/Amplitude"

    amplitude_measurement.plot_divergence(simpath, Tc=1.7477, direction="xiy")

if __name__ == "__main__":
    main()