from FunctionsAndClasses import *
import functools
import matplotlib; matplotlib.use("TkAgg")


def start_sinus(event, ax, phi_t):
    """Plotte Sinus-Kurve."""
    x = np.linspace(0.0, 2.0*np.pi, 100)
    sinus = ax.plot(x, np.sin(x)) # Anfangsplot
    print("Hi")
    for phi in phi_t:
        sinus_t = np.sin(x-phi) # Neue Daten
        sinus[0].set_ydata(sinus_t) # Plotdaten aktualisieren,
        event.canvas.flush_events() # und dynamisch
        event.canvas.draw() # darstellen.

def main():
    """Hauptprogramm."""
    phi_t = np.linspace(0.0, 6.0, 100) # Phasenverschiebung
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Maustaste zum Starten klicken")
    klick_funktion = functools.partial(start_sinus, ax=ax,
    phi_t=phi_t)
    fig.canvas.mpl_connect("button_press_event", klick_funktion)
    plt.show()

    print("hä?")
    x = np.linspace(0, 10, 10, endpoint=True)

    y = x * x

    dy_dx = second_order_num_diff(x, y)
    dy_dx_np = np.gradient(y, x)
    dy_dx_analytic = 2 * x

    print(dy_dx)
    print(dy_dx_np)
    print(dy_dx_analytic)

if __name__ == "__main__":
    main()