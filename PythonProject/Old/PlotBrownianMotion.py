import csv
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
# import matplotlib; matplotlib.use("TkAgg")

def read_csv(filepath):
    df = pd.read_csv(filepath, header=None)
    return df

def plot_csv(df, ax, fig=None):
    ax.plot(df[0], df[1])
    ax.set_xlabel("t/s")
    ax.set_ylabel("x")


def read_multiple_csv(filepaths):
    trajectories = []
    for filepath in filepaths:
        df = read_csv(filepath)
        trajectories.append(df)
    return trajectories

def plot_multiple(trajectories):
    fig, ax  = plt.subplots(1, 1)
    for tra in trajectories:
        plot_csv(tra, ax)

    plt.show()

def main():
    filepath = "../LearningProject/Brownian Motion.csv"
    root = "../LearningProject/Brownian Motion in Potential Data/"
    filepaths = []
    for filename in os.listdir(root):
        f = os.path.join(root, filename)
        if os.path.isfile(f):
            filepaths.append(f)

    trajectories = read_multiple_csv(filepaths)
    plot_multiple(trajectories)


if __name__ == "__main__":
    main()
