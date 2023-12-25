from FunctionsAndClasses import *

def main():
    simulation_path = "../../Generated content/Silicon/Quench/Left-Right"

    t_xix = {}
    t_xiy = {}

    for size in os.listdir(simulation_path):
        if size != "plots":
            sizepath = os.path.join(simulation_path, size)
            if os.path.isdir(sizepath):
                for temp in os.listdir(sizepath):
                    temppath = os.path.join(sizepath, temp)
                    if os.path.isdir(temppath):
                        t_xix[temp] = {}
                        t_xiy[temp] = {}
                        ft_k, ft_l = average_ft(temppath)
                        if temp == "128.000000":
                            print(ft_l)
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

    print(t_xiy)
    print("\n\n")
    print(t_xiy.keys())

    fig, ax = plt.subplots(1, 1)
    for temp in t_xix:
        t = list(t_xix[temp].keys())
        xix = list(t_xix[temp].values())
        ax.plot(t, xix, linestyle="", marker=".", ms=2, label=rf"$\xi_x$  T = {float(temp):.2f}")
    configure_ax(fig, ax)
    plt.show()

    fig, ax = plt.subplots(1, 1)
    for temp in t_xiy:
        t = list(t_xix[temp].keys())
        xiy = list(t_xiy[temp].values())
        ax.plot(t, xiy, linestyle="", marker=".", label=rf"$\xi_y$  T = {float(temp):.2f}")
    configure_ax(fig, ax)
    plt.show()



if __name__ == "__main__":
    main()