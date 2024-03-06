from FunctionsAndClasses import *

def main():
    # I think it is easiest here just to copy the values, I will only need
    # this 3 times or so

    h = [0.1, 0.4161791450287818, 1.7320508075688776, 7.208434242404265, 10]
    Tc = [0.7714, 0.8951, 1.1813, 1.3423, 1.2741]

    fig, ax = plt.subplots(1, 1)

    ax.plot(Tc, h)
    configure_ax(fig, ax)
    plt.show()

if __name__ == "__main__":
    main()