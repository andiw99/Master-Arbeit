from FunctionsAndClasses import *
from itertools import product
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def signal(x):
    return np.exp(-x**2/2)


def calc_corr_func(x):
    """
    calculates lookup arrays for the correlation function for the possible
    distances in a lattice with size n x n, values stored in x
    :param x: values of the order parameter
    :return: array[distance], array[value_corr_func]
    """

    n = x.shape[0]
    dists = distances_in_2D_lattice(n)

    corr_func = np.zeros(len(dists))        # initialize array with len of
                                            # possible distances

    # For every dist we have to calculate the correlation function
    # TODO when calculating the possible disances we already know which
    # combinations of delta i and delta j result in the distance, we could
    # remember this?
    # TODO i think this is the most inefficient calculation ever but lets see
    # if the performance is okay
    for (ind, dist) in enumerate(dists):
        # find all lattice combinations and add them up
        corr_value = 0
        nr_calcs = 0
        # look at all starting lattice sites
        for (i, j) in product(range(n), range(n)):
            # Now we want to iter through all coming lattice sites:
            for (k, h) in product(range(n), range(n)):
                # check if they are really later sites:
                # I Think this calculation is even wrong at least for PBC
                # but lets
                # see if it shows the right scaling
                if (k * n + h >= i * n + j):
                    # check if we got the correct distance
                    # i think we can just use == since we compare integers
                    actual_dist = squared_dist((i, j), (k, h))
                    if actual_dist == dist:
                        corr_value += x[i][j] * x[k][h]
                        nr_calcs += 1
        corr_value /= nr_calcs          # averaging
        corr_func[ind] = corr_value     # saving

    return dists, corr_func


def calc_corr_func_v2(x):
    """
    calculates lookup arrays for the correlation function for the possible
    distances in a lattice with size n x n, values stored in x
    :param x: values of the order parameter
    :return: array[distance], array[value_corr_func]
    """

    n = x.shape[0]
    dist_tuples = distance_tuples(n)
    dists = distances_in_2D_lattice(n)

    corr_func = np.zeros(len(dists))        # initialize array with len of
                                            # possible distances
    # we iter through all tuples, adding up the ones with the same distance
    ind = 0
    corr_ind = 0
    while ind < len(dist_tuples):
        dist = dist_tuples[ind][0]
        corr_value = 0
        nr_calcs = 0
        # Solange die distanz des nächsten tuples sich nicht verändert gehen
        # wir durch die tuple
        while (ind < len(dist_tuples)) and (dist_tuples[ind][0] == dist):
            # extract the differences we are looking at
            di = dist_tuples[ind][1]
            dj = dist_tuples[ind][2]
            # now we go through every lattice site and calc the correlation
            # with all neighbors that have di, dj
            for (i, j) in product(range(n), range(n)):
                # think of the modulo operation?
                # TODO problem we value distances with di = dj stronger
                # TODO und es passt auch sonst gar nicht weil di = 19
                # wieder nachbar ist
                corr_value += x[i][j] * x[(i + di) % n][(j + dj) % n]
                corr_value += x[i][j] * x[(i - di) % n][(j + dj) % n]
                corr_value += x[i][j] * x[(i + di) % n][(j - dj) % n]
                corr_value += x[i][j] * x[(i - di) % n][(j - dj) % n]
                nr_calcs += 4
            ind += 1
        # We have updated ind by 1 and now we have another dist meaning we break
        # now we save the calculated value
        corr_func[corr_ind] = corr_value / nr_calcs
        corr_ind += 1
        # print(ind, dist_tuples[ind][0], corr_ind, dists[corr_ind])

    return dists, corr_func


def calc_corr_manhatten(x, cutoff=0):
    """
    calculates the correlation function using the manhatten metrik
    :param x: same as above
    :return: same as above
    """
    # TODO neue Idee für corr func berechnung:
    # 1.) Gehe jeden Gitterplatz i,j durch
    # 2.) Berechne Distanz zu jedem anderen Gitterplatz k,h
    # 3.) Addiere corr((i,j), (k,h)) zu corr_func[squared_dist((i, j),(k, h))]
    # 4.) Average am Ende
    # Eventuell musst du dir dafür die additionen auch merken
    # die possible distances reichen jetzt von 0 bis (n-1) + (n-1)
    n = x.shape[0]
    # possible dists in one direction
    # is halfed because of pbc
    dists = np.arange(0, np.minimum(int(n/2), cutoff//2 + 1))
    # with manhatten we can have double that so we need twice as much space
    corr_func = np.zeros(2 * len(dists) - 1)
    nr_calcs = np.zeros(2 * len(dists) - 1)
    for (i, j) in product(range(n), range(n)):
        for dist_i in dists:
            for dist_j in dists:
                # with periodic boundary conditions we will be fine if we only
                # look at the +dist terms?
                # the -dist terms are included by x[i-dist][j-dist]
                corr_func[dist_i + dist_j] += \
                    x[i][j] * x[(i + dist_i) % n][(j + dist_j) % n]
                # keep track of number of calculations
                nr_calcs[dist_i + dist_j] += 1
    # average
    corr_func /= nr_calcs
    # TODO normalize to 1?
    corr_func /= np.max(corr_func)
    return np.arange(2 * len(dists) - 1), corr_func


def transform_corr_func(corr_func):
    """
    Transforms the correlation function in a way so that we only take te values
    which are greater than zero and after that just set every value to be zero
    :param dists:
    :param corr_func:
    :return:
    """
    # index at which the correlation function is closest to zero
    # leider zu einfach
    # ind = np.argmin(np.abs(corr_func))
    # corr_func[ind:] = 0
    # Wenn es keinen Vorzeichenwechsel gibt ist ind der inex vom kleinsten eintrag
    ind = np.argmin(np.abs(corr_func))
    for i in range(len(corr_func)-1):
        if(corr_func[i] > 0) and corr_func[i+1] < 0:
            ind = i
            break

    corr_func[ind:] = 0
    return corr_func



def calc_corr(x, cutoff, dist_func):
    # 1.) Gehe jeden Gitterplatz i,j durch
    # 2.) Berechne Distanz zu jedem anderen Gitterplatz k,h
    # 3.) Addiere corr((i,j), (k,h)) zu corr_func[squared_dist((i, j),(k, h))]
    # problem: squared dist ist nicht äquidistant aber wir können:
    # corr_func[dists.ind_of_arg(squared_dist((i, j), (k, h)] machen
    # 4.) Average am Ende
    # Eventuell musst du dir dafür die additionen auch merken
    # die possible distances reichen jetzt von 0 bis (n-1) + (n-1)
    n = x.shape[0]
    squared_dists = distances_in_2D_lattice(n)
    dists = list(squared_dists)     # so that we can use the index option
    dists = dists[:dists.index(cutoff ** 2) + 1]
    # with manhatten we can have double that so we need twice as much space
    corr_func = np.zeros(len(dists))
    nr_calcs = np.zeros(len(dists))
    for (i, j) in product(range(n), range(n)):
        for (k, h) in product(range(n), range(n)):
            # Now calculating without PBC
            curr_dist = dist_func((i, j), (k, h), n)
            if curr_dist <= (cutoff ** 2):
                dist_ind = dists.index(curr_dist)
                corr_func[dist_ind] += \
                    x[i][j] * x[k][h]
                nr_calcs[dist_ind] += 1
    # average
    corr_func /= nr_calcs
    # TODO normalize to 1?
    corr_func /= np.max(corr_func)
    return dists, corr_func


def distance_tuples(n):
    """
    calculates the possible distances in a 2D lattice
    :param n: n x n lattice
    :return: set of possible distances
    """

    dists = set()

    for i in range(n):
        for j in range(n):
            dist = i ** 2 + j ** 2
            dists.add((dist, i, j))
    return sorted(dists)


def squared_dist(coords1, coords2, n):
    squared_dist = (coords1[0] - coords2[0]) ** 2 +\
                   (coords1[1] - coords2[1]) ** 2
    return squared_dist


def squared_dist_pbc(coords1, coords2, n):
    #  TODO there has to be a better way to calc this?
    di_abs = abs(coords1[0] - coords2[0])
    di = min(di_abs, n - di_abs)
    dj_abs = abs(coords1[1] - coords2[1])
    dj = min(dj_abs, n - dj_abs)
    return di ** 2 + dj ** 2


def distances_in_2D_lattice(n):
    """
    calculates the possible distances in a 2D lattice
    :param n: n x n lattice
    :return: set of possible distances
    """

    dists = set()

    for i in range(n):
        for j in range(i + 1):
            dist = (i) ** 2 + j ** 2
            dists.add(dist)
    return sorted(dists)


def fourier_transform(corr_func, dists):
    struct_fact = np.abs(np.fft.fft(corr_func))
    q_values = np.fft.fftfreq(len(dists), d=dists[1] - dists[0])
    # shift the struct fact to zero?
    # struct_fact -= np.min(struct_fact)
    return q_values, struct_fact


def gaussian(x, a, x0, sigma):
    # Define a Gaussian function
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def fwhm(x, y):
    # Fit a Gaussian distribution to the data and extract the FWHM
    p0 = [np.max(y), x[np.argmax(y)], np.std(y)]
    popt, _ = curve_fit(gaussian, x, y, p0=p0)
    return 2 * np.sqrt(2 * np.log(2)) * np.abs(popt[2])


def exp_decay(r, xi):
    return np.exp(-r/xi)



def calc_corr_length(filepath):
    df = read_csv(filepath)
    n2 = int((df.shape[1] - 2) / 2)
    n = int(np.sqrt(n2))
    op_values = df.iloc[-1].iloc[1:n2 + 1]
    op_values = np.array(op_values, dtype=float).reshape((n, n))
    # squared_dists, corr_func = calc_corr_func(op_values)
    # plt.plot(np.sqrt(squared_dists), corr_func)
    cutoff = n // 2
    squared_dists, corr_func = calc_corr_manhatten(op_values, cutoff)
    corr_func = transform_corr_func(corr_func)
    dists = np.sqrt(squared_dists)
    dists **= 2
    q_values, struct_fact = fourier_transform(corr_func, dists)
    struct_fact = np.abs(struct_fact[np.argsort(q_values)])
    q_values = sorted(q_values)
    curve_width = fwhm(q_values, struct_fact)
    xi = 1 / curve_width
    n = op_values.shape[0]
    print("FWHM: ", curve_width)
    print("Correlation Length", xi)
    print("Number of Domains ", n ** 2 / xi ** 2)
    return (dists, corr_func), (q_values, struct_fact), xi




def main():
    root = "../../Generated content/Cooling Bath/Comparison/eta=5.00/T=70.00/dt=0.0050/n=75" \
           "/alpha=5.00/beta=10.00/J=50.00/"
    name = "tau=200.00/0.csv"
    filepath = root + name
    # read in the file and extract last row
    filepaths = new_files_in_dir(root, root, plot_all=True)
    # leeren container für xi und tau
    data = {
        "xi" : [],
        "xi_inter": [],
        "tau": [],
    }
    for filepath in filepaths:
        # we need to read in the parameters
        para_path = os.path.splitext(filepath)[0] + ".txt"
        parameters = read_parameters_txt(para_path)
        (dists, corr_func), (q_values, struct_fact), xi = calc_corr_length(filepath)
        # interpolieren für zusätzliche stützstellen:

        f = interp1d(dists, corr_func, kind="quadratic")
        additional_dists = np.linspace(0, np.max(dists), 100)
        additional_corr_func = f(additional_dists)
        add_q, add_struct = fourier_transform(additional_corr_func, additional_dists)
        add_struct = add_struct[np.argsort(add_q)]
        add_q = sorted(add_q)
        inter_xi = 1/fwhm(add_q, add_struct)
        print("correlation length with interpolated values", inter_xi)
        # den daten hinzufügen
        data["xi"].append(xi)
        data["tau"].append(parameters["tau"])
        data["xi_inter"].append(inter_xi)
        # plotten, plot speichern?
        fig, axes = plt.subplots(2, 1, figsize=(6, 7))
        axes[0].set_ylabel("C(r)")
        axes[0].set_xlabel("r")
        axes[0].plot(dists, corr_func, label="Corr. Function")
        axes[0].plot(additional_dists, additional_corr_func,
                     label="Interpolation")
        axes[0].plot(dists, exp_decay(dists, xi), label="Exponential Fit")
        axes[0].plot(dists, exp_decay(dists, inter_xi), label="Exp. Fit with Inter")
        axes[0].legend()
        axes[1].set_ylabel("S(q)")
        axes[1].set_xlabel("q")
        axes[1].plot(q_values, np.abs(struct_fact), label="Fourier Transform")
        axes[1].plot(add_q, add_struct, label="Interpolated")
        axes[0].set_title("Correlation Function with Manhattan Metrik")
        axes[1].legend()
        plt.show()

    sortargs = np.argsort(data["tau"])
    print(sortargs)
    data["tau"] = sorted(data["tau"])
    data["xi"] = np.array(data["xi"])[sortargs]
    data["xi_inter"] = np.array(data["xi_inter"])[sortargs]

    # write to csv:
    df = pd.DataFrame(data)
    df.to_csv(root + "Corr Length.txt", index=False)

    # plot it

    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\xi$")
    ax.plot(data["tau"], data["xi"], label="corr length")
    ax.plot(data["tau"], data["xi_inter"], label="with interpolation")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()
