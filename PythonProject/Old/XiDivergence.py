from matplotlib import pyplot as plt
import numpy as np
from FunctionsAndClasses import *

def xi_div(eps, xi0, nu):
    return xi0 / (np.abs(eps) ** nu)

def xi_flattened(eps, xi0, nu, omega):
    return eps ** (omega) * ( 1 + xi0 * np.abs(eps) ** nu)

def simple_lorentz(eps, xi0, eps0, k):
    return xi0 / (1 + k * (eps - eps0) ** 2)

def main():
    res = 50
    U_xi = 1.95
    #f_minus = 1
    #f_plus = U_xi * f_minus
    f_minus = 3
    f_plus = 1
    nu = 1
    omega = 2
    eps0 = -0.4

    eps_below = np.linspace(-1, 0, res, endpoint=False)
    eps_above = np.sort(np.linspace(1, 0, int(res * 1.5), endpoint=False))
    eps_flat = np.linspace(-1, 1.5, 2 * res)
    xi_below = xi_div(eps_below, f_minus, nu)
    xi_above = xi_div(eps_above, f_plus, nu)
    # xi_flat_below = xi_flattened(eps_below, f_minus, nu, omega)
    # xi_flat_above = xi_flattened(eps_above, f_plus,  nu, omega)
    xi_flat = simple_lorentz(eps_flat, 20 * f_minus, eps0, 10)

    eps = np.concatenate((eps_below, eps_above))
    xi = np.concatenate((xi_below, xi_above))

    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
    ax.plot(eps, xi, label=r"$\xi = f^\pm / \varepsilon^{\nu}$    ")
    # ax.plot([], [], label=r"$\xi_L$", color=colors[4])
    # ax.plot(eps_flat, xi_flat, label=r"$\xi = f^\pm / \varepsilon^{\nu}$")
    ax.set_ylim((0.01, 25))
    ax.set_xlim((-0.8, 0.8))
    ax.set_xlabel(r"$\frac{T - T_c}{T_c}$", fontsize=14)
    ax.set_ylabel(r"$\xi$", fontsize=14)
    ax.set_title(r"critical divergence of the correlation length $\xi$")

    config = {
        "nr_x_major_ticks": 10,
        "nr_y_major_ticks": 4,
        "legendfontsize": 14
    }
    configure_ax(fig, ax, config)
    plt.savefig("xi-divergence.pdf", format="pdf")
    # plt.savefig("xi-divergence", format="png", dpi=1000)
    # save_plot("./", "xi-divergence.pdf", "png", dpi=1000)
    plt.show()

    res = 500
    f_minus = 0.008
    U_xi = 1
    f_plus = U_xi * f_minus
    nu = 1
    z = 2
    mu = nu * z
    eps_below = np.linspace(-0.9, 0, res, endpoint=False)
    eps_above = np.sort(np.linspace(1., 0, res , endpoint=False))
    tau_below = xi_div(eps_below, f_minus, mu)
    tau_above = xi_div(eps_above, f_plus, mu)

    t = np.concatenate((eps_below, eps_above))
    tau = np.concatenate((tau_below, tau_above))
    fig, ax = plt.subplots(1, 1, figsize=(1 + 8, 1 + 2))
    ax.plot(t, tau, label=r"$\tau$")
    ax.plot(t, np.abs(t), linestyle="--", label=r"$\tau(t) = t$")
    ax.set_xlabel(r"$t$", fontsize=14)
    ax.set_ylabel(r"$\tau$", fontsize=14)
    ax.set_ylim(0, 0.4)
    ax.set_xlim(-0.8, 0.8)
    inter = get_intersection_index(t, tau)
    # mark_point(ax, t[inter], tau[inter])
    # mark_point(ax, -t[inter], tau[inter], c="C0")
    # find out gradient at t[inter]
    dtau_dt = np.gradient(tau, t)
    slope = -dtau_dt[inter]

    t_linear =  np.linspace(-t[inter], t[inter], res)
    xi_linear = tau[inter] + slope * (t_linear + t[inter])

    i_inter, j_inter = get_intersection_index(xi_linear[int(res * 0.1):], tau, t_linear[int(res * 0.1):], t)
    i_inter += int(res * 0.1)
    print(i_inter , j_inter)
    print(xi_linear[i_inter], tau[j_inter])

    #ax.plot(t_linear[:i_inter], xi_linear[:i_inter], c="black")



    ax.set_title(r"  ")

    config = {
        "nr_y_major_ticks" : 2,
        "nr_x_major_ticks" : 10,
        "legendfontsize" : 14
    }

    configure_ax(fig, ax, config)
    save_plot("./", "xi-divergence-symmetric.pdf", "pdf")
    plt.show()


if __name__ == "__main__":
    main()