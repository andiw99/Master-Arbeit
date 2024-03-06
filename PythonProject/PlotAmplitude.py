from Suite import amplitude_measurement

def main():
    simpath = "../../Generated content/Silicon/Subsystems/Suite/h/0.4161791450287818/Amplitude/"

    amplitude_measurement.plot_divergence(simpath, Tc=0.895)

if __name__ == "__main__":
    main()