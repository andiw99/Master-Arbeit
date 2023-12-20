from FunctionsAndClasses import *

def main():
    simulation_path = "../../Generated content/Test/FtAvgTest"

    t_xix = {}
    t_xiy = {}

    for size in os.listdir(simulation_path):
        sizepath = os.path.join(simulation_path, size)
        for temp in os.listdir(sizepath):
            t_xix[temp] = {}
            t_xiy[temp] = {}
            temppath = os.path.join(sizepath, temp)
            ft_k, ft_l = average_ft(temppath)
            for t in ft_k:
                p_k = get_frequencies_fftw_order(len(ft_k[t]))
                popt_x, perr_x = fit_lorentz(p_k, ft_k[t])
                xix = np.abs(popt_x[0])
                t_xix[temp][t] = xix
            for t in ft_l:
                p_l = get_frequencies_fftw_order(len(ft_l[t]))
                popt_y, perr_y = fit_lorentz(p_l, ft_l[t])
                xiy = np.abs(popt_y[0])
                t_xiy[temp][t] = xiy

    print(t_xix)

    fig, ax = plt.subplots(1, 1)
    for temp in t_xix:
        t = list(t_xix[temp].keys())
        xix = list(t_xix[temp].values())
        print(t)
        ax.plot(t, xix, label=rf"$\xi_x$  T = {float(temp):.2f}")
    configure_ax(fig, ax)
    plt.show()

    fig, ax = plt.subplots(1, 1)
    for temp in t_xiy:
        t = list(t_xix[temp].keys())
        xiy = list(t_xiy[temp].values())
        print(t)
        ax.plot(t, xiy, label=rf"$\xi_y$  T = {float(temp):.2f}")
    configure_ax(fig, ax)
    plt.show()



if __name__ == "__main__":
    main()