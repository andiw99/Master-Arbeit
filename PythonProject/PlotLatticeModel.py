import csv
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
# import matplotlib; matplotlib.use("TkAgg")

def read_csv(filepath):
    df = pd.read_csv(filepath, header=None)
    return df

def plot_csv(df, ax, fig=None, row=1):
    ax.plot(df[0], df[row])
    ax.set_xlabel("t/s")
    ax.set_ylabel("x")
    ax.set_ylim(-1.2, 1.2)


def read_multiple_csv(filepaths):
    trajectories = []
    for filepath in filepaths:
        df = read_csv(filepath)
        trajectories.append(df)
    return trajectories

def plot_multiple(trajectories, nr_subplots=2):
    fig, axes  = plt.subplots(2, 2)
    for tra in trajectories:
        i = 1
        for axs in axes:
            for ax in axs:
                plot_csv(tra, ax, row=i)
                i += 1
    plt.tight_layout()
    plt.show()

def main():
    filepath = "../LearningProject/Brownian Motion.csv"
    root = "../LearningProject/Brownian Motion Stochastic Lattice/"
    filepaths = []
    for filename in os.listdir(root):
        f = os.path.join(root, filename)
        if os.path.isfile(f):
            filepaths.append(f)

    trajectories = read_multiple_csv(filepaths)
    plot_multiple(trajectories)


if __name__ == "__main__":
    main()
