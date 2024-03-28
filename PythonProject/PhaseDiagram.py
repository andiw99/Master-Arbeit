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

def determine_phase_diagram_ising(BJ_min, BJ_max, res):
    BJ_x = np.linspace(BJ_min, BJ_max, res)         # determine BJ_x
    # determine the BJ_y in dependence of BJ_x
    # we should have two branches? always?
    # how do we save the critical points? Just two arrays?
    BJ_y_pos = []      # this is for the case that BJ_y is smaller than BJ_x
    BJ_y_neg = []      # this is for the case that BJ_y is larger than BJ_x

    for BJ in BJ_x:
        BJ_y = T_c_est_Ising_eff(BJ)
        BJ_y_pos.append(BJ_y)
        BJ_y_neg.append(-BJ_y)

    return BJ_x, BJ_y_pos, BJ_y_neg

def main():
    BJ_min = -6.0
    BJ_max = 6.0
    res = 1000

    BJ_x, BJ_y_pos, BJ_y_neg = determine_phase_diagram(BJ_min, BJ_max, res)

    print(BJ_x)
    print(BJ_y_pos)
    print(BJ_y_neg)
    fig, ax = plt.subplots(1, 1, figsize=(8 , 8))
    ax.plot(BJ_x, BJ_y_pos, color="C1", label="Phase Boundary")
    ax.plot(BJ_x, BJ_y_neg, color="C1")
    # the path that you quench is the path with BJ_para / B_J_perp = 31  BJ_para = 31 BJ_perp
    BJ_perp = 10 * np.linspace(-9 / 31, 0, res)
    BJ_para = 31 * BJ_perp
    ax.plot(BJ_perp, BJ_para, color="C0", label="Quench path")

    tangent = -12.34 * BJ_perp - 7.03
    ax.plot(BJ_perp, tangent, color=colors[5], label="tangent")

    ax.set_ylim(BJ_min, BJ_max)
    ax.set_xlim(BJ_min, BJ_max)
    ax.set_ylabel(r"$\beta J_\parallel$")
    ax.set_xlabel(r"$\beta J_\perp$")
    ax.set_title("Phase Diagram of the XY-Model")
    configure_ax(fig, ax)
    plt.show()

    res = 501
    BJ_min = -2.0
    BJ_max = 2.0
    BJ_x, BJ_y_pos, BJ_y_neg = determine_phase_diagram_ising(BJ_min, BJ_max, res)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.plot(BJ_x, BJ_y_pos, color="C1", label="Phase Boundary")
    ax.plot(BJ_x, BJ_y_neg, color="C1")
    # the path that you quench is the path with BJ_para / B_J_perp = 31  BJ_para = 31 BJ_perp
    BJ_perp = 1 * np.linspace(-9 / 31, 0, res)
    BJ_para = 31 * BJ_perp
    ax.plot(BJ_perp, BJ_para, color="C0", label="Quench path")

    # derivative line
    tangent = -10.25 * BJ_perp - 2.01
    ax.plot(BJ_perp, tangent, color=colors[5], label="tangent")
    ax.set_ylim(BJ_min, BJ_max)
    ax.set_xlim(BJ_min, BJ_max)
    ax.set_ylabel(r"$\beta J_\parallel$")
    ax.set_xlabel(r"$\beta J_\perp$")
    ax.set_title("Phase Diagram of the Ising-Model")
    configure_ax(fig, ax)
    plt.show()

    BJ_min = -6.0
    BJ_max = 6.0
    res = 4501
    BJ_x_Ising, BJ_y_pos_Ising, BJ_y_neg_Ising = determine_phase_diagram(BJ_min, BJ_max, res)
    BJ_x, BJ_y_pos, BJ_y_neg = determine_phase_diagram_ising(BJ_min, BJ_max, res)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.plot(BJ_x, BJ_y_pos, color=colors[4], label="Ising")
    ax.plot(BJ_x_Ising, BJ_y_pos_Ising, color=colors[0], label="XY")
    ax.plot(BJ_x, BJ_y_neg, color=colors[4])
    ax.plot(BJ_x_Ising, BJ_y_neg_Ising, color=colors[0])

    # the path that you quench is the path with BJ_para / B_J_perp = 31  BJ_para = 31 BJ_perp
    BJ_perp = 1 * np.linspace(-9 / 20, 0, res)
    BJ_para = 15 * BJ_perp
    ax.plot(BJ_perp, BJ_para, linewidth=2, color="C1", label=r"$J_\parallel / J_\perp = 15$")

    # derivative line
    #tangent = -10.25 * BJ_perp - 2.01
    #ax.plot(BJ_perp, tangent, color=colors[5], label="tangent")
    ax.set_ylim(BJ_min, BJ_max)
    ax.set_xlim(BJ_min, BJ_max)
    ax.set_ylabel(r"$\beta J_\parallel$")
    ax.set_xlabel(r"$\beta J_\perp$")
    ax.set_title("Phase Diagram of the Ising-Model")

    config = {
        "increasefontsize": 0.75,
        "labelrotation": 90,
        "labelverticalalignment": "bottom"
    }

    configure_ax(fig, ax, config)
    fig.savefig("./Phase-diagram-Ising-Model.svg", format="svg")
    plt.show()

if __name__ == "__main__":
    main()