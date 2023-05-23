import numpy as np


def construct_ev(N, i):
    j = np.arange(1, N + 1)
    a = np.sin(i * j * np.pi / (N + 1))
    # corresponding inverse EV
    b = 2 / (N + 1) * np.sin(i * j * np.pi / (N + 1))
    return a, b




def main():
    N = 10000
    a, b = construct_ev(N, 1)
    a2, b2 = construct_ev(N, 2)

    print(np.sum(a * b))
    print(np.sum(a2 * b))

if __name__ == "__main__":
    main()