from FunctionsAndClasses import *


def crit_h(T, A, B, C, D):
    tmp1 = np.exp(B / T)
    tmp2 = - A * T ** 2
    return D + C * np.exp(tmp2 * tmp1)


def main():
    # I think it is easiest here just to copy the values, I will only need
    # this 3 times or so
    J_para = 3.11

    h = [0.01, 0.1, 0.2, 0.4161791450287818, 0.8, 1.7320508075688776, 3, 7.208434242404265]#, 10]
    Tc = [0.704, 0.7714, 0.8318, 0.8951, 1.03, 1.1813, 1.3154, 1.3399] #1.2741]
    Tc = np.array(Tc)
    Tc_fit = np.array(Tc) #- 0.618

    Tc = Tc / J_para

    A, B, C = 0.1, 1, 5

    popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    h_est = crit_h(Tc_fit, *popt)

    fig, ax = plt.subplots(1, 1)

    ax.plot(Tc, h, marker='s', **blue_point_kwargs, label="crit. points")
    ax.plot(Tc, h_est, label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h$ / meV")
    ax.set_xlabel("$T_c / J_\parallel$")

    config = {
        "labelrotation": 90,
        "legendfontsize": 15,
    }

    configure_ax(fig, ax, config)
    plt.show()

if __name__ == "__main__":
    main()