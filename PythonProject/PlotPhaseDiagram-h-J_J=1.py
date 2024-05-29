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
    J_para = 1
    h =  np.array([0.01, 0.05, 0.1, 0.25, 0.5])#, 10]
    Tc = np.array([0.978, 1.0211, 1.0494, 1.1010, 1.1255]) #1.2741]


    h /= Tc

    Tc /= J_para

    A, B, C, D = 1, 1, 100, 0

    # popt, _ = curve_fit(crit_h, Tc_fit[:-1], h[:-1], maxfev=10000)
    TcXY = 1 / (2 * lambertw(1))
    Tc_fit = Tc
    h_fit = h
    popt, _ = curve_fit(crit_h, Tc_fit, h_fit, maxfev=100000, p0=(A, B, C, D))
    Tc_plot = np.linspace(TcXY, np.max(Tc_fit), 100)
    h_est = crit_h(Tc_plot, *popt)
    print(h_est)
    print("correct popt?")
    print(popt)

    h_given = 0.416
    J_given = 3

    # h_given = 0.283
    # J_given = 3.11

    h_given = 1
    h_given = np.linspace(0.5, 0.25, 1000)
    J_given = 3.11
    T_given = np.linspace(0.4 * J_given, 0.2 * J_given, 1000)
    h_T_given_plot = h_given / T_given      # those are my x values
    T_given_J_given = T_given / J_given     # those are the y values
    # h path?
    h_T_quench_path = 1 - T_given_J_given

    # Angle
    h_est_angle = crit_h(T_given_J_given, *popt)

    # angle_int, _ = find_first_intersection(T_given_J_given[::-1], T_given_J_given[::-1], h_T_given_plot[::-1], h_est_angle[::-1])
    # print("Path and PB intersect at T / J = ", angle_int)
    # angle = angle_between_curves(T_given_J_given[::-1], h_T_given_plot[::-1], h_est_angle[::-1], angle_int)
    # print("angle = ", angle)
    # angle_int2, _ = find_first_intersection(T_given_J_given[::-1], T_given_J_given[::-1], h_T_quench_path[::-1], h_est_angle[::-1])
    # angle2 = angle_between_curves(T_given_J_given[::-1], h_T_quench_path[::-1], h_est_angle[::-1], angle_int2)


    fig, ax = plt.subplots(1, 1, figsize=(6, 12))

    #ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points",)
    ax.plot(Tc[:-1], h[:-1] , marker='s', **blue_point_kwargs, markeredgewidth=1, markersize=7, label="$j_\parallel = 3.11$")
    ax.plot(Tc_plot, h_est, label="$h / T \propto e^{-AT^2 e^{B / T}}$", color=colors[1])
    ax.plot(TcXY, 0, marker="^", **get_point_kwargs_color(colors[0], markeredgewidth=1), markersize=7,)
    #ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7,)
    ax.plot([], [], marker="^", **get_point_kwargs_color("black"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    #ax.plot(Tc100 / 10, h100 / Tc100, marker="s", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7, label="$J_\parallel / J_\perp = 100$")
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    #ax.plot(T_given_J_given, h_T_given_plot, color="C1", alpha=0.5, label=f"$\phi = {angle}$")
    #ax.plot(T_given_J_given, h_T_quench_path, color="C2", alpha=0.5, label=f"$\phi_2 = {angle2}$")
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    # ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # point from J/J = 100
    # ax.plot(0.17, 0.245, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    #ax.plot(Tc[:-1], h_est2[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    #ax.plot(Tc[:-1], h_est3[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h~/~T$")
    ax.set_xlabel("$T~/~J_\parallel$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "labelhorizontalalignment": "center",
        "labelverticalalignment": "bottom",
        #"labelhorizontalalignment": "right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-h(T).png", dpi=500, format="png")
    plt.savefig("phase_transition-h(T)-svg.svg", format="svg")

    plt.show()

    # plot over J / T
    #angle_int, _ = find_first_intersection(1 / T_given_J_given[::-1], 1 / T_given_J_given[::-1], h_T_given_plot[::-1], h_est_angle[::-1])
    #angle = angle_between_curves(1 / T_given_J_given[::-1], h_T_given_plot[::-1], h_est_angle[::-1], angle_int)
#
    #angle_int2, _ = find_first_intersection(1 / T_given_J_given[::-1],
    #                                        1 / T_given_J_given[::-1],
    #                                        h_T_quench_path[::-1],
    #                                        h_est_angle[::-1])
    #angle2 = angle_between_curves(1 / T_given_J_given[::-1], h_T_quench_path[::-1],
    #                              h_est_angle[::-1], angle_int2)

    fig, ax = plt.subplots(1, 1, figsize=(6, 12))

    #ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points",)
    ax.plot(1 / Tc[:-1], h[:-1] , marker='s', **blue_point_kwargs, markeredgewidth=1, markersize=7, label="$j_\parallel = 3.11$")
    ax.plot(1 / Tc_plot, h_est, label="$h / T \propto e^{-AT^2 e^{B / T}}$", color=colors[1])
    ax.plot(1 / TcXY, 0, marker="^", **get_point_kwargs_color(colors[0], markeredgewidth=1), markersize=7,)
    #ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7,)
    ax.plot([], [], marker="^", **get_point_kwargs_color("black"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    #ax.plot(Tc100 / 10, h100 / Tc100, marker="s", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7, label="$J_\parallel / J_\perp = 100$")
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    # ax.plot(1 / T_given_J_given, h_T_given_plot, color="C1", alpha=0.5, label=f"$\phi = {angle}$")
    # ax.plot(1 / T_given_J_given, h_T_quench_path, color="C2", alpha=0.5, label=f"$\phi = {angle2}$")

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    # ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # point from J/J = 100
    # ax.plot(0.17, 0.245, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    #ax.plot(Tc[:-1], h_est2[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    #ax.plot(Tc[:-1], h_est3[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h~/~T$")
    ax.set_xlabel("$J_\parallel / T$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "labelhorizontalalignment": "center",
        "labelverticalalignment": "bottom",
        #"labelhorizontalalignment": "right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-h(T).png", dpi=500, format="png")
    plt.savefig("phase_transition-h(T)-svg.svg", format="svg")

    plt.show()

    # plot versus J_perp / T
    angle_int, _ = find_first_intersection(1 / T_given_J_given[::-1] / 31, 1 / T_given_J_given[::-1] / 31, h_T_given_plot[::-1], h_est_angle[::-1])
    angle = angle_between_curves(1 / T_given_J_given[::-1] / 31, h_T_given_plot[::-1], h_est_angle[::-1], angle_int)
    angle_int2, _ = find_first_intersection(1 / T_given_J_given[::-1] / 31,
                                            1 / T_given_J_given[::-1] / 31,
                                            h_T_quench_path[::-1],
                                            h_est_angle[::-1])
    angle2 = angle_between_curves(1 / T_given_J_given[::-1] / 31, h_T_quench_path[::-1],
                                  h_est_angle[::-1], angle_int2)

    fig, ax = plt.subplots(1, 1, figsize=(6, 12))

    # ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points",)
    ax.plot(1 / Tc[:-1] / 31, h[:-1], marker='s', **blue_point_kwargs,
            markeredgewidth=1, markersize=7, label="$j_\parallel = 3.11$")
    ax.plot(1 / Tc_plot / 31, h_est, label="$h / T \propto e^{-AT^2 e^{B / T}}$",
            color=colors[1])
    ax.plot(1 / TcXY / 31, 0, marker="^",
            **get_point_kwargs_color(colors[0], markeredgewidth=1),
            markersize=7, )
    # ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7,)
    ax.plot([], [], marker="^", **get_point_kwargs_color("black"),
            label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # ax.plot(Tc100 / 10, h100 / Tc100, marker="s", **get_point_kwargs_color("C1", markeredgewidth=1), markersize=7, label="$J_\parallel / J_\perp = 100$")
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.plot(1 / T_given_J_given / 31, h_T_given_plot, color="C1", alpha=0.5,
            label=f"$\phi = {angle}$")
    ax.plot(1 / T_given_J_given / 31, h_T_quench_path, color="C2", alpha=0.5, label=f"$\phi = {angle2}$")

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    # ax.plot(TcXY100, 0, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # point from J/J = 100
    # ax.plot(0.17, 0.245, marker="^", **get_point_kwargs_color("C4"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # ax.plot(Tc[:-1], h_est2[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    # ax.plot(Tc[:-1], h_est3[:-1], label="$h \propto e^{-AT^2 e^{B / T}}$")
    ax.set_ylabel("$h~/~T$")
    ax.set_xlabel("$J_\perp / T$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "labelhorizontalalignment": "center",
        "labelverticalalignment": "bottom",
        # "labelhorizontalalignment": "right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-h(T).png", dpi=500, format="png")
    plt.savefig("phase_transition-h(T)-svg.svg", format="svg")

    plt.show()

    exit()

    Tc_desired = 0.232


    # I want to have h / T in the range that I can cover which is from 0 to two
    # so T goes from infinity to 0.2

    # so from those i can get the intersection
    h_T_given_point, T_J_given_point = find_first_intersection(h_T_given_plot, h, T_given_J_given, Tc)

    # We want to show the original critical temperature as horizontal lign
    h_to_use = find_x_star(h, Tc, Tc_desired)

    fig, ax = plt.subplots(1, 1)

    ax.plot([], [], marker='s', **get_point_kwargs_color("black"), label="critical points")
    ax.plot(h[:-1], Tc[:-1], marker='s', **blue_point_kwargs, label="$j_\parallel = 3.1$")



    #ax.plot(h2, Tc2, marker='s', **get_point_kwargs_color(colors[4]), label="$J_{exp} 60$")
    ax.plot(h3, Tc3, marker='s',
            **get_point_kwargs_color(colors[5]), label="$j_\parallel =  110 \cdot 10^3$")
    #ax.plot(h_est, Tc_plot, color="C1")
    #ax.set_xscale("log")
    xlims = ax.get_xlim()
    print(xlims)
    ax.plot(0, TcXY, marker="^", **get_point_kwargs_color(colors[0]))

    ax.plot(0, TcXY100, marker="^", **get_point_kwargs_color("C1"))
    ax.plot([], [], marker="^", **get_point_kwargs_color("black"), label="$T_c^{\mathrm{XY}}~ / ~J_\parallel $")
    # point from J/J = 100
    ax.plot(0.245, 0.17, marker="s", **get_point_kwargs_color("C1"))
    ax.plot(h100 / Tc100, Tc100 / 10, marker="s", **get_point_kwargs_color("C1"), label="$J_\parallel / J_\perp = 100$")

    ax.set_xlim(xlims)
    ylims = ax.get_ylim()
    ax.plot(h_T_given_plot, T_given_J_given,  color="C3", alpha=0.5)
    ax.set_ylim(ylims)

    #mark_point(ax, h_T_given_point, T_J_given_point, label=f"$T / J_\parallel = {T_J_given_point:.3f}$")
    # here we plot this given plot thing

    ax.set_xlabel("$h~/~T$")
    ax.set_ylabel("$T~/~J_\parallel$")

    config = {
        "labelrotation": 90,
        "increasefontsize": 0.3,
        "labelhorizontalalignment": "center",
        "labelverticalalignment": "bottom",
        "legendlocation": "lower right"
    }

    configure_ax(fig, ax, config)
    plt.savefig("phase_transition-T(h).png", dpi=500, format="png")
    plt.show()

if __name__ == "__main__":
    main()