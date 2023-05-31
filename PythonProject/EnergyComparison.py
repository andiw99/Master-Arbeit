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




def main():
    N = 10000
    T = 30
    J = 10

    print(classical_energy(N, T))
    print(qm_energy(J, N, T))

if __name__ == "__main__":
    main()