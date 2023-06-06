import numpy as np
from matplotlib import pyplot as plt


def classical_energy(N, T):
    return N * T


def w_k(k, J, N):
    return np.sqrt(4 * J * (1 + np.cos(np.pi * k / (N + 1))))


def qm_energy(J, N, T):
    k = np.arange(1, N + 1)
    w = w_k(k, J, N)

    E = np.sum(1/2 * w * (np.exp(w / T) + 1) / (np.exp(w / T) - 1))

    return E


def zero_point_energy(J, N):
    k = np.arange(1, N + 1)
    w = w_k(k, J, N)
    E = np.sum(1/2 * w)
    return E


def main():
    N = 10000
    T = 100
    J = 10

    T_vals = np.linspace(0.5, 100, 1000)

    E_class = np.zeros(len(T_vals))
    E_qm = np.zeros(len(T_vals))
    J_T = np.zeros(len(T_vals))

    for i, T_val in  enumerate(T_vals):
        E_class[i] = classical_energy(N, T_val)
        E_qm[i] = qm_energy(J, N, T_val)
        J_T[i] = J / T_val

        if((E_qm[i] - E_class[i]) / E_qm[i] > 0.02):
            J_T_krit = J_T[i]
            krit_ind = i

    fig, ax = plt.subplots(1, 1)
    ax.plot(J_T, E_class, label="classical energy")
    ax.plot(J_T, E_qm, label = "qm energy")
    ax.plot(J/T_vals[[0, -1]], [zero_point_energy(J, N), zero_point_energy(J, N)], ls="dashed", c="grey", label="zero point energy")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("E")
    ax.set_xlabel("J/T")
    ax.set_title("comparison between classical and qm oscillator chain")
    ax.axvline(J_T_krit, ls=":", c="grey", label=f"2% error at J/T = {J_T_krit:.2f}")
    ax.legend()
    plt.show()

    print(E_class[krit_ind])
    print(E_qm[krit_ind])

if __name__ == "__main__":
    main()