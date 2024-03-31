import numpy as np

import matplotlib.pyplot as plt
from FunctionsAndClassesLight import *
from glob import glob


def main():
    plot_eta_05_h_1()
    #plot_eta_5_h_1()
    #plot_eta_1_h_05()
    #plot_1_h_1()
    #plot_eta_1_h_5()
    exit()
    #plot_eta_1_h_5()
    simpath = "../../Generated content/Final/z-measurement-small/0.5/z"

    Tc = 1.725
    fold =  200
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()
def plot_eta_1_h_10():
    #simpath = "../../Generated content/Final/z-measurement-small/0.5/z"
    simpath = "/bigdata/StrongFieldQED/Weitzel/z-measurement-small/h=10/10/z"

    Tc = 3.3250
    fold =  200
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()
def plot_eta_1_h_5():
    #simpath = "../../Generated content/Final/z-measurement-small/0.5/z"
    simpath = "/bigdata/StrongFieldQED/Weitzel/z-measurement-small/h=5-3/5/z"

    Tc = 2.875
    fold =  200
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()
def plot_eta_1_h_10():
    #simpath = "../../Generated content/Final/z-measurement-small/0.5/z"
    simpath = "/bigdata/StrongFieldQED/Weitzel/Content without sync/z-measurement-small/h=1010/z"

    Tc = 3.325
    fold =  200
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=Tc, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()
def plot_eta_5_h_1():
    simpath = "/bigdata/StrongFieldQED/Weitzel/Content without sync/z-measurement-small/eta=5/1/z"

    Tc = 1.975000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 7, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")
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
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()

def plot_eta_05_h_1():
    simpath = "../../Generated content/Final/z-measurement-small/eta=0.5/1/z"

    Tc = 1.970000
    fold = 80
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()

def plot_1_h_1():
    simpath = "/bigdata/StrongFieldQED/Weitzel/Content without sync/z-measurement-small/1-2/z"

    Tc = 1.975000
    fold = 120
    folder_avg_function=read_folder_avg

    sizes = np.linspace(72, 120, 3, dtype=np.int64)
    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()


def plot_eta_1_h_05():
    simpath = "../../Generated content/Final/z-measurement-small/0.5/z"
    simpath = "../../../Generated content without sync/Final/z-measurement-small/0.5/z"
    Tc = 1.731000
    fold = 40
    folder_avg_function=read_folder_avg

    sizes = np.linspace(48, 144, 4, dtype=np.int64)

    size_cum_dic, size_times_dic = get_results_time_resolved(sizes, simpath, Tc=None, file_ending="mag", value_name="U_L",
                                  process_folder_avg_func=folder_avg_function)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold)
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.png", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-plotfile.svg", format="svg")


    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold // 4, xlim=0.15)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25.svg", format="svg")

    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.8 / 6.4 * 10))
    fit_and_plot(fig, ax, size_cum_dic, size_times_dic, fold_nr=fold, xlim=0.25)
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.png-plotfile", format="png")
    fig.savefig(simpath + f"/cum-over-time-scan-0.25-more-fold.svg", format="svg")

    plt.show()



if __name__ == "__main__":
    main()