import matplotlib.pyplot as plt
import numpy as np

from FunctionsAndClasses import *

def determine_phase_diagram(BJ_min, BJ_max, res):
    BJ_x = np.linspace(BJ_min, BJ_max, res)         # determine BJ_x
    # determine the BJ_y in dependence of BJ_x
    # we should have two branches? always?
    # how do we save the critical points? Just two arrays?
    BJ_y_pos = []      # this is for the case that BJ_y is smaller than BJ_x
    BJ_y_neg = []      # this is for the case that BJ_y is larger than BJ_x

    for BJ in BJ_x:
        BJ_small_y = BJ_small(BJ)
        BJ_large_y = BJ_large(BJ)

        if np.abs(BJ_large_y) < np.abs(BJ):
            BJ_y_pos.append(BJ_small_y)
            BJ_y_neg.append(-BJ_small_y)
        else:
            BJ_y_pos.append(BJ_large_y)
            BJ_y_neg.append(-BJ_large_y)

        trans_eq_small = BJ_assert(BJ, BJ_small_y)
        if trans_eq_small > 1e-5:
            print(f"Something wrong for BJ_x = {BJ}, BJ_y = {BJ_small_y}")
            print(f"transcendet equation equals = {trans_eq_small}")
        trans_eq_large = BJ_assert(BJ, BJ_large_y)
        if trans_eq_large > 1e-5:
            print(f"Something wrong for BJ_x = {BJ}, BJ_y = {BJ_large_y}")
            print(f"transcendet equation equals = {trans_eq_large}")
    return BJ_x, BJ_y_pos, BJ_y_neg

def main():
    BJ_min = -10.0
    BJ_max = 10.0
    res = 500

    BJ_x, BJ_y_pos, BJ_y_neg = determine_phase_diagram(BJ_min, BJ_max, res)

    print(BJ_x)
    print(BJ_y_pos)
    print(BJ_y_neg)
    fig, ax = plt.subplots(1, 1, figsize=(8 * np.max(BJ_x) / np.max(BJ_y_pos), 8))
    ax.plot(BJ_x, BJ_y_pos, color="C1", label="Phase Boundary")
    ax.plot(BJ_x, BJ_y_neg, color="C1")
    # the path that you quench is the path with BJ_para / B_J_perp = 31  BJ_para = 31 BJ_perp
    BJ_perp = np.linspace(-9 / 31, 0, res)
    BJ_para = 31 * BJ_perp
    ax.plot(BJ_perp, BJ_para, res, color="C0", label="Quench path")
    ax.set_ylim(-10, 10)
    ax.set_ylabel(r"$\beta J_\parallel$")
    ax.set_xlabel(r"$\beta J_\perp$")
    ax.set_title("Phase Diagram of the XY-Model")
    configure_ax(fig, ax)
    plt.show()

if __name__ == "__main__":
    main()