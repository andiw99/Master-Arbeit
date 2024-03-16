from FunctionsAndClasses import *
from scipy.special import lambertw

def crit_h(T, A, B, C, D):
    tmp1 = np.exp(B / T)
    tmp2 = - A * T ** 2
    return C * np.exp(tmp2 * tmp1)

def find_x_star(x_values, y_values, y_star):
    # Interpolate the curve
    f = interp1d(y_values, x_values, kind='linear')

    # Find x where y(x) = y*
    x_star = f(y_star)

    return x_star
def main():
    # I think it is easiest here just to copy the values, I will only need
    # this 3 times or so
    """
    J_para = 3.11
    Tc_XY = 0.618
    J_para2 = 120000
    J_para3 = 93000
    Tc_XY2 = 20020
    Tc_XY3 = 18500
    h = [0.01, 0.1, 0.2, 0.4161791450287818, 0.8, 1.7320508075688776, 3, 7.208434242404265]#, 10]
    Tc = [0.704, 0.7714, 0.85, 0.91, 1.03, 1.1813, 1.3154, 1.3399] #1.2741]
    h2 = [2800, 10000]
    Tc2 = [22750, 27200]

    h3 = [930, 5500]
    Tc3 = [21500, 24695]

    h2 = np.array(h2) /     Tc_XY2
    Tc2 = np.array(Tc2) /   Tc_XY2

    h3 = np.array(h3) / Tc_XY3
    Tc3 = np.array(Tc3) / Tc_XY3

    Tc_fit = np.array(Tc) / Tc_XY

    h = np.array(h)     / Tc_XY
    Tc = np.array(Tc)   / Tc_XY

    A, B, C = 1, 1, 100

   # popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    print(popt)
    #A, B, C, D = popt
    #print(f"A = {A}, B = {B}, C = {C}, D = {D}")
    h_est = crit_h(Tc_fit * Tc_XY, *popt)

    fig, ax = plt.subplots(1, 1)

    ax.plot(Tc[:-1], h[:-1], marker='s', **blue_point_kwargs, label="crit. points")
    ax.plot(Tc[:-1], h_est[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.plot(Tc3, h3, marker='s', **get_point_kwargs_color(colors[5]), label="$J_{exp} 31$")
    ax.set_ylabel("$h / J_\parallel$")
    ax.set_xlabel("$T_c / J_\parallel$")

    config = {
        "labelrotation": 90,
        "legendfontsize": 15,
    }
    #ax.set_yscale("log")

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-h(T).png", dpi=500, format="png")
    plt.show()

    Tc_desired = 1.232

    h_to_use = find_x_star(h, Tc, Tc_desired)

    fig, ax = plt.subplots(1, 1)

    ax.plot(h, Tc, marker='s', **blue_point_kwargs, label="crit. points")
    ax.plot(h2, Tc2, marker='s', **get_point_kwargs_color(colors[4]), label="$J_{exp} 60$")
    ax.plot(h3, Tc3, marker='s', **get_point_kwargs_color(colors[5]), label="$J_{exp} 31$")
    mark_point(ax, h_to_use, Tc_desired, label=f"$h/ J_\parallel = {h_to_use:.3f}$")
    ax.set_xscale("log")
    ax.set_xlabel("$h / J_\parallel$")
    ax.set_ylabel("$T_c / J_\parallel$")

    config = {
        "labelrotation": 90,
        "legendfontsize": 15,
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-T(h).png", dpi=500, format="png")
    plt.show()

"""


    J_para = 3.11
    Tc_XY = 0.618
    J_para2 = 120000
    J_para3 = 93000
    Tc_XY2 = 20020
    Tc_XY3 = 18500
    h =  [0.01, 0.1, 0.2, 0.4161791450287818, 0.8, 1.7320508075688776, 3, 7.208434242404265]#, 10]
    Tc = [0.704, 0.7714, 0.84, 0.91, 1.03, 1.1813, 1.3154, 1.3399] #1.2741]
    h2 = [2800, 10000]
    Tc2 = [22750, 27200]

    h3 = [930, 5500, 17000]
    Tc3 = [21500, 24695, 28800]

    h2 = np.array(h2)
    Tc2 = np.array(Tc2, dtype=np.float64)
    h3 = np.array(h3)
    Tc3 = np.array(Tc3, dtype=np.float64)
    Tc_fit = np.array(Tc)
    h = np.array(h)
    Tc = np.array(Tc)


    h /= Tc

    h3 =  h3 / Tc3
    h2 =  h2 / Tc2


    Tc /= J_para
    Tc2 /= J_para2
    Tc3 /= J_para3

    A, B, C, D = 1, 1, 100, 0

    # popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    TcXY = 1 / (2 * lambertw(31))
    Tc_fit = np.concatenate((Tc[:-1], Tc3))
    h_fit = np.concatenate((h[:-1], h3))
    popt, _ = curve_fit(crit_h, Tc_fit, h_fit, maxfev=100000, p0=(A, B, C, D))
    Tc_plot = np.linspace(TcXY, np.max(Tc_fit), 100)
    h_est = crit_h(Tc_plot, *popt)
    #h_est2 = crit_h(Tc_fit, 1.2, 1.215, 100)
    #h_est3 = crit_h(Tc_fit, 2, 1.215, 100)
    print(h_est)
    print("correct popt?")
    print(popt)



    fig, ax = plt.subplots(1, 1)

    ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points")
    ax.plot(Tc[:-1], h[:-1] , marker='s', **blue_point_kwargs, label="$j_\parallel = 3.1$")
    ax.plot(Tc3, h3, marker='s', **get_point_kwargs_color(colors[5]), label="$j_\parallel =  110 \cdot 10^3$")
    ax.plot(Tc_plot, h_est, label="$h \propto e^{-AT^2 e^{B / T}}$", color="C1")
    ax.plot(TcXY, 0, marker="^", **get_point_kwargs_color("C1"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")

    #ax.plot(Tc[:-1], h_est2[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    #ax.plot(Tc[:-1], h_est3[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h~/~T$")
    ax.set_xlabel("$T~/~J_\parallel$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "labelhorizontalalignment": "right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-h(T).png", dpi=500, format="png")
    plt.show()

    Tc_desired = 0.232

    # We want to show the original critical temperature as horizontal lign
    h_to_use = find_x_star(h, Tc, Tc_desired)

    fig, ax = plt.subplots(1, 1)

    ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points")
    ax.plot(h[:-1], Tc[:-1], marker='s', **blue_point_kwargs, label="$j_\parallel = 3.1$")



    #ax.plot(h2, Tc2, marker='s', **get_point_kwargs_color(colors[4]), label="$J_{exp} 60$")
    ax.plot(h3, Tc3, marker='s',
            **get_point_kwargs_color(colors[5]), label="$j_\parallel =  110 \cdot 10^3$")
    #ax.plot(h_est, Tc_plot, color="C1")
    #mark_point(ax, h_to_use, Tc_desired, label=f"$h/ J_\parallel = {h_to_use:.3f}$")
    #ax.set_xscale("log")
    xlims = ax.get_xlim()
    print(xlims)
    ax.plot(0, TcXY, marker="^", **get_point_kwargs_color("C1"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    ax.set_xlim(xlims)

    ax.set_xlabel("$h~/~T$")
    ax.set_ylabel("$T~/~J_\parallel$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "legendlocation": "lower right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-T(h).png", dpi=500, format="png")
    plt.show()

if __name__ == "__main__":
    main()