import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def exp(x, xi):
    return np.exp(-x/xi)

def gaussian(x, a, x0, sigma):
    # Define a Gaussian function
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def fwhm(x, y):
    # Fit a Gaussian distribution to the data and extract the FWHM
    p0 = [np.max(y), x[np.argmax(y)], np.std(y)]
    popt, _ = curve_fit(gaussian, x, y, p0=p0)
    return 2 * np.sqrt(2 * np.log(2)) * np.abs(popt[2])

def main():
    ns = [11, 100, 1000]
    end_x = 20
    for n in ns:
        x = np.linspace(0, end_x, n)
        xi = 1

        y = exp(x, xi)

        y_fft = np.fft.fft(y)
        k = np.fft.fftfreq(len(x), d=x[1] - x[0])

        y_fft = np.abs(y_fft[np.argsort(k)])
        k = sorted(k)

        fig, axes = plt.subplots(2, 1, figsize=(6, 7))
        axes[0].set_ylabel("C(r)")
        axes[0].set_xlabel("r")
        axes[0].plot(x, y, label="Function")
        axes[1].plot(k, y_fft, label="FFT Function")
        axes[1].set_ylabel("S(q)")
        axes[1].set_xlabel("q")
        axes[0].set_title("Correlation Function with Manhattan Metrik")

        width = fwhm(k, y_fft)
        xi_estimate = 1/width
        axes[0].plot(x, exp(x, xi_estimate), label="fit")
        axes[0].legend()

        plt.show()

        print("FWHM: ", width)

if __name__ == "__main__":
    main()