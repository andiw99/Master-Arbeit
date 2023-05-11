import numpy as np
from matplotlib import pyplot as plt

def double_well(x, alpha, beta):
    return 1/2 * alpha*  (x ** 4 - beta * x **2)


def kos():
    icks_frequency = 1

    # Plot points
    fig, ax = plt.subplots(figsize=(10, 10))


    # Set identical scales for both axes

    # Set bottom and left spines as x and y axes of coordinate system
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)



    # Create custom major ticks to determine position of tick labels

    # Create minor ticks placed at each integer to enable drawing of minor grid
    # lines: note that this has no effect in this example with ticks_frequency=1

    # Draw major and minor grid lines
    ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

    # Draw arrows
    arrow_fmt = dict(markersize=4, color='black', clip_on=False)
    ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
    ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

    return fig, ax


def plot_potential(alpha, beta):
    xmin = np.sqrt(beta / 2)
    Vmin = - 1/8 * alpha * beta ** 2

    x = np.linspace(-1.7 * xmin, 1.7 * xmin, 500)

    V = double_well(x, alpha, beta)

    # Enter x and y coordinates of points and colors


    fig, ax = kos()

    ax.plot(x, V, c="C1", lw=2)
    ax.scatter(xmin, Vmin)
    ax.plot([xmin, xmin], [0, Vmin], c="C0", ls='--', lw=2, alpha=0.75)
    ax.plot([0, xmin], [Vmin, Vmin], c="C0", ls='--', lw=2, alpha=0.75)
    # Create 'x' and 'y' labels placed at the end of the axes
    ax.set_xlabel('q', size=14, labelpad=-35, x=1.00)
    ax.set_ylabel('V(q)', size=14, labelpad=-50, y=0.98, x=0.5, rotation=0)
    plt.tight_layout()


def plot_potential_interaction(alpha, beta, J):
    xmin = np.sqrt(beta / 2)
    Vmin = - 1/8 * alpha * beta ** 2

    x = np.linspace(-1.7 * xmin, 1.7 * xmin, 500)

    V = double_well(x, alpha, beta) + 2 * J * (x - xmin) ** 2

    # Enter x and y coordinates of points and colors


    fig, ax = kos()

    ax.plot(x, V, c="C1", lw=2)
    ax.scatter(xmin, Vmin)
    ax.plot([xmin, xmin], [0, Vmin], c="C0", ls='--', lw=2, alpha=0.75)
    ax.plot([0, xmin], [Vmin, Vmin], c="C0", ls='--', lw=2, alpha=0.75)
    # Create 'x' and 'y' labels placed at the end of the axes
    ax.set_xlabel('q', size=14, labelpad=-35, x=1.00)
    ax.set_ylabel('V(q)', size=14, labelpad=-50, y=0.98, x=0.5, rotation=0)
    plt.tight_layout()


def plot_transformation(a):

    eps = 0.01

    x = np.linspace(-1 + eps, 1 + eps, 500)

    f = a * np.arctanh(x)
    f2 = 3 * a * np.arctanh(x)

    fig, ax = kos()

    ax.plot(x, f, c="C1", label=f"a = {a:.2f}")
    limits = ax.get_ylim()
    ax.plot(x, f2, c="C0", label=f"a = {3*a:.2f}")
    ax.set_ylim(limits)
    ax.set_xlabel(r"$\frac{2 \vartheta}{\pi}$", size=14, labelpad=-5, x=1.00)
    ax.set_ylabel(r"$q(\vartheta)$", size=14, labelpad=-50, y=0.98, x=0.5, rotation=0)
    ax.legend(fontsize=14)



def main():
    alpha = 10
    beta = 1
    a = 1
    J = alpha * beta / 17
    #plot_potential(alpha, beta)
    #plot_transformation(a)
    plot_potential_interaction(alpha, beta, J)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()