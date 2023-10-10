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
    ln_f = np.log(f)
    ln_t = np.log(t)
    gamma = second_order_num_diff(ln_t, ln_f)

    return gamma



def main():
    root = "../../Generated content/NER Long/Selected"
    root_dirs = list_directory_names(root)
    file_extension = ".ner"
    print(root_dirs)
    ppf = 20

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
                    print("adding")
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
        ax.plot(t_dic[key][1:], m_dic[key][1:], ls="", marker="x", color=colors[i], label=f"T = {float(key):.3f}")
    configure_ax(fig, ax)
    plt.show()
    # fitting gamma_m
    for i, key in enumerate(m_dic.keys()):
        gamma_arr, t_arr = fit_gamma_m(t_dic[key], m_dic[key], ppf=ppf)
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
        ax.plot(np.array(t_dic[key][1:]), f_mm_dic[key][1:], ls="", marker="x",
                color=colors[(2 * i) % len(colors)], label=f"T = {float(key):.3f}")
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
        gamma_mm = calc_gamm_fluctuations(f_mm_dic[key][1:], t_dic[key][1:])
        ax.plot(1 / np.array(t_dic[key][1:]), gamma_mm, ls="", marker="x", color=colors[i], label=f"T = {float(key):.3f}")
    ax.set_yscale("log")
    ax.set_xscale("log")
    configure_ax(fig, ax)
    plt.show()

    fig, ax = plt.subplots(1, 1)
    for i, key in enumerate(m_dic.keys()):
        gamma_me = calc_gamm_fluctuations(f_me_dic[key][1:], t_dic[key][1:])
        ax.plot(1 / np.array(t_dic[key][1:]), gamma_me, ls="", marker="x", color=colors[i], label=f"T = {float(key):.3f}")
    ax.set_yscale("log")
    ax.set_xscale("log")
    configure_ax(fig, ax)
    plt.show()
    exit()




if __name__ == "__main__":
    main()