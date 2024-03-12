from FunctionsAndClasses import *


def crit_h(T, A, B, C, D):
    tmp1 = np.exp(B / T)
    tmp2 = - A * T ** 2
    return D + C * np.exp(tmp2 * tmp1)

def find_x_star(x_values, y_values, y_star):
    # Interpolate the curve
    f = interp1d(y_values, x_values, kind='linear')

    # Find x where y(x) = y*
    x_star = f(y_star)

    return x_star
def main():
    # I think it is easiest here just to copy the values, I will only need
    # this 3 times or so
    J_para = 3.11
    Tc_XY = 0.618
    J_para2 = 120000
    Tc_XY2 = 20020
    h = [0.01, 0.1, 0.2, 0.4161791450287818, 0.8, 1.7320508075688776, 3, 7.208434242404265]#, 10]
    Tc = [0.704, 0.7714, 0.8318, 0.8951, 1.03, 1.1813, 1.3154, 1.3399] #1.2741]
    h2 = [10000]
    Tc2 = [27200]

    h2 = np.array(h2) /     Tc_XY2
    Tc2 = np.array(Tc2) /   Tc_XY2

    Tc_fit = np.array(Tc) #- 0.618

    h = np.array(h)     / J_para
    Tc = np.array(Tc)   / Tc_XY

    A, B, C = 0.1, 1, 5

    popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    h_est = crit_h(Tc_fit, *popt)

    fig, ax = plt.subplots(1, 1)

    ax.plot(Tc, h, marker='s', **blue_point_kwargs, label="crit. points")
    ax.plot(Tc, h_est, label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h / J_\parallel$")
    ax.set_xlabel("$T_c / J_\parallel$")

    config = {
        "labelrotation": 90,
        "legendfontsize": 15,
    }

    configure_ax(fig, ax, config)
    plt.show()

    Tc_desired = 1.59

    h_to_use = find_x_star(h, Tc, Tc_desired)

    fig, ax = plt.subplots(1, 1)

    ax.plot(h, Tc, marker='s', **blue_point_kwargs, label="crit. points")
    ax.plot(h2, Tc2, marker='s', **get_point_kwargs_color(colors[4]), label="$J_{exp}$")
    mark_point(ax, h_to_use, Tc_desired, label=f"$h/ J_\parallel = {h_to_use:.3f}$")
    ax.set_xscale("log")
    ax.set_xlabel("$h / J_\parallel$")
    ax.set_ylabel("$T_c / J_\parallel$")

    config = {
        "labelrotation": 90,
        "legendfontsize": 15,
    }

    configure_ax(fig, ax, config)
    plt.show()

if __name__ == "__main__":
    main()