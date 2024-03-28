import numpy as np

import matplotlib.pyplot as plt
from FunctionsAndClassesLight import *
from glob import glob


def main():

    simpath = "../../Generated content/Final/z-measurement-small/h=5/5/z"

    Tc = 2.85
    fold = 40
    folder_avg_function=process_folder_avg_balanced

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_5_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=5/1/z"

    Tc = 1.975000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 7, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_02_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=0.2-2/1/z"

    Tc = 1.970000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(80, 208, 3, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_eta_05_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=0.5/1/z"

    Tc = 1.970000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()

def plot_1_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/1/z"

    Tc = 1.975000
    fold = 120
    folder_avg_function=process_folder_avg_balanced

    sizes = np.linspace(72, 120, 3, dtype=np.int64)
    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 2, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()


def plot_eta_1_h_05():
    simpath = "../../Generated content/Final/z-measurement-small/0.5/z"

    Tc = 1.731000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 2, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    plt.show()



if __name__ == "__main__":
    main()