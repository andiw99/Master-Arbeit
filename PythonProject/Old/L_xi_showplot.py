from matplotlib import pyplot as plt
import numpy as np
from FunctionsAndClasses import *

def L_xi_func(eps, A, B, L, nu):
    return A * eps ** 2 + B * eps * L ** (1/nu)

def main():
    res = 150
    A = 2
    B = 1
    nu = 2
    omega = 2
    eps0 = -0.4
    Ls = [10, 100, 1000]
    L_xis = []
    eps = np.linspace(-1, 1.5, res)
    for L in Ls:
        L_xis.append(L_xi_func(eps, A, B, L, nu))


    fig, ax = plt.subplots(1, 1)
    for L_xi in L_xis:
        ax.plot(eps, L_xi)

    ax.set_xlabel(r"$\frac{T - T_c}{T_c}$", fontsize=14)
    ax.set_title(r"Finite Size Scaling of $L / \xi$")
    configure_ax(fig, ax)
    save_plot("./", "L_xi-FFS.pdf", "pdf")
    plt.show()

if __name__ == "__main__":
    main()