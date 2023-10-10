from FunctionsAndClasses import *
from scipy.optimize import curve_fit


def fit_gamma_m(t, m, ppf=5):
    dt = ppf * (t[1] - t[0])    # ppf: points per fit
    gamma_arr = []
    t_arr = []

    nr_fits = len(t) // ppf

    # doing the fitting for every ppf points
    for i in range(nr_fits):
        lnt_fit = np.log(t[i * ppf:(i+1) * ppf])
        lnm_fit = np.log(m[i * ppf:(i+1) * ppf])

        popt, pcov = curve_fit(linear_fit, lnt_fit, lnm_fit)
        gamma_arr.append(-popt[0])
        t_arr.append(t[i * ppf + 2])

    return gamma_arr, t_arr


def calc_gamm_fluctuations(f, t):
    # We will do some kind of averaging so that the derivative does not
    # oscillate as strongly
    ln_f = np.log(f)
    ln_t = np.log(t)
    #gamma = second_order_num_diff(ln_t, ln_f)
    gamma = np.gradient(ln_f, ln_t)
    # TODO Just use gradient for now since you do not have equidistant
    # points so the central difference implementation becomes a bit harder
    return gamma


def avg(x, y, ppp):
    nr_points = len(x) // ppp
    y_avg = np.zeros(nr_points)
    x_avg = np.zeros(nr_points)

    for i in range(nr_points):
        for j in range(ppp):
            y_avg[i] += y[i * ppp + j]
            x_avg[i] += x[i * ppp + j]
        y_avg[i] /= ppp
        x_avg[i] /= ppp

    return x_avg, y_avg

def smooth(x, y, ppp):
    # ppp should be uneven and at least 3
    nr_points = len(x)
    y_smooth = np.zeros(nr_points)
    # smoothing works only if we got (ppp - 1) / 2 points before and after the
    # concerned point
    print(y)
    for i in range(nr_points):
        # we can also try to calculate the possible smoothing range
        pp = np.minimum(ppp, (np.minimum(2 * i, 2 * (nr_points - i - 1)) + 1))
        a = (pp - 1) // 2
        for j in range(pp):
            y_smooth[i] += y[i - a + j]
        y_smooth[i] /= pp
    return y_smooth

def main():
    root = "../../Generated content/NER Long/Selected"
    root_dirs = list_directory_names(root)
    file_extension = ".ner"
    print(root_dirs)
    ppf = 10
    ppp = 5
    ppp_smooth = 1
    times_smoothing = 0

    m_dic = {}
    t_dic = {}
    f_mm_dic = {}
    f_me_dic = {}
    gamma_m_dic = {}
    t_gamma_m_dic = {}

    for temp_dir in root_dirs:
        temp_path = os.path.join(root, temp_dir)

        # now we need all ner files in this folder
        files = os.listdir(temp_path)
        print(files)
        runs = 0
        for file in files:
            if os.path.splitext(file)[1] == file_extension:
                runs += 1
                file_path = os.path.join(temp_path, file)
                df = read_struct_func(file_path)
                t = np.array(df["t"])
                m = np.array(df["m"])
                f_mm = np.array(df["f_mm"])
                f_me = np.array(df["f_me"])

                if temp_dir in m_dic:
                    m_dic[temp_dir] += m
                    f_mm_dic[temp_dir] += f_mm
                    f_me_dic[temp_dir] += f_me
                else:
                    m_dic[temp_dir] = m
                    f_mm_dic[temp_dir] = f_mm
                    f_me_dic[temp_dir] = f_me
                    t_dic[temp_dir] = t
        print("nr of files read: ", runs)
        m_dic[temp_dir] /= runs     # averaging
        f_mm_dic[temp_dir] /= runs
        f_me_dic[temp_dir] /= runs

    # plotting

    fig, ax = plt.subplots(1, 1)

    ax.set_yscale("log")
    ax.set_xscale("log")
    for i, key in enumerate(m_dic.keys()):
        t_avg, m_avg = avg(t_dic[key][1:], m_dic[key][1:], 1)
        for j in range(times_smoothing):
            m_avg = smooth(t_avg, m_avg, ppp_smooth)
        ax.plot(t_avg, m_avg, ls="", marker="x", color=colors[i], label=f"T = {float(key):.3f}")
    configure_ax(fig, ax)
    plt.show()
    # fitting gamma_m
    for i, key in enumerate(m_dic.keys()):
        t_avg, m_avg = avg(t_dic[key][1:], m_dic[key][1:], 1)
        for j in range(times_smoothing):
            m_avg = smooth(t_avg, m_avg, ppp_smooth)
        gamma_arr, t_arr = fit_gamma_m(t_avg, m_avg, ppf=ppf)
        gamma_m_dic[key] = gamma_arr
        t_gamma_m_dic[key] = t_arr

    # plotting gamma_m
    fig, ax = plt.subplots(1, 1)

    for i, key in enumerate(m_dic.keys()):
        ax.plot(1 / np.array(t_gamma_m_dic[key][1:]), gamma_m_dic[key][1:], ls="", marker="x", color=colors[(2 * i) % len(colors) ], label=f"T = {float(key):.3f}")
    configure_ax(fig, ax)
    plt.show()

    # plotting f_mm and f_me
    fig, ax = plt.subplots(1, 1)

    for i, key in enumerate(m_dic.keys()):
        t_avg, f_mm_avg = avg(t_dic[key][1:], f_mm_dic[key][1:], ppp)
        for j in range(times_smoothing):
            f_mm_avg = smooth(t_avg, f_mm_avg, ppp_smooth)
        ax.plot(t_avg, f_mm_avg, ls="", marker="x",
                color=colors[(i) % len(colors)], label=f"T = {float(key):.3f}")
    ax.set_ylabel(r"$f_{mm}(t)$")
    ax.set_xlabel(r"t")
    ax.set_yscale("log")
    ax.set_xscale("log")
    configure_ax(fig, ax)

    fig, ax = plt.subplots(1, 1)

    for i, key in enumerate(m_dic.keys()):
        ax.plot(np.array(t_dic[key][1:]), f_me_dic[key][1:], ls="", marker="x",
                color=colors[(2 * i) % len(colors)], label=f"T = {float(key):.3f}")
    ax.set_ylabel(r"$f_{me}(t)$")
    ax.set_xlabel(r"t")
    ax.set_yscale("log")
    ax.set_xscale("log")
    configure_ax(fig, ax)
    plt.show()

    # calcing and plotting gamma_mm

    fig, ax = plt.subplots(1, 1)
    for i, key in enumerate(m_dic.keys()):
        t_avg, f_mm_avg = avg(t_dic[key][1:], f_mm_dic[key][1:], ppp)
        for j in range(times_smoothing):
            f_mm_avg = smooth(t_avg, f_mm_avg, ppp_smooth)
        gamma_mm= calc_gamm_fluctuations(f_mm_avg, t_avg)
        ax.plot(1 / t_avg, gamma_mm, ls="", marker="x", color=colors[2 * i], label=f"T = {float(key):.3f}")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(r"$\lambda_{mm}$")
    ax.set_xlabel(r"$\frac{1}{t}$")
    configure_ax(fig, ax)
    plt.show()

    fig, ax = plt.subplots(1, 1)
    for i, key in enumerate(m_dic.keys()):
        t_avg, f_me_avg = avg(t_dic[key][1:], f_me_dic[key][1:], ppp)
        for j in range(times_smoothing):
            f_mm_avg = smooth(t_avg, f_mm_avg, ppp_smooth)
        gamma_me = calc_gamm_fluctuations(f_me_avg, t_avg)
        ax.plot(1 / t_avg, gamma_me, ls="", marker="x", color=colors[2 * i], label=f"T = {float(key):.3f}")
    ax.set_yscale("log")
    ax.set_xscale("log")
    configure_ax(fig, ax)
    plt.show()
    exit()




if __name__ == "__main__":
    main()