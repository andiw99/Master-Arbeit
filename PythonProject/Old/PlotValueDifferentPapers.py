import matplotlib.pyplot as plt
import numpy as np
from FunctionsAndClasses import *

def main():

    # Data for experiments
    experiments = {
        "MC\n[Ref 1]": {"z": 2.167, "error": 0},
        "RG\n[Ref 2]": {"z": 2.14, "error": 0.02},
        "HT\n[Ref 3]": {"z": 2.183, "error": 0.005},
        "RG\n[Ref 4]": {"z": 2.15, "error": 0}
    }

    # Extract x-axis labels, corresponding values, and errors
    x_labels = list(experiments.keys())
    z_values = [experiment["z"] for experiment in experiments.values()]
    # z_values_error = [experiment["z"] for experiment in experiments.values() if experiment["error"] > 0]
    # z_values_no_error = [experiment["z"] for experiment in experiments.values() if experiment["error"] == 0]
    errors = [experiment["error"] for experiment in experiments.values()]

    # Calculate mean value
    mean_z = np.mean(z_values)

    # Plotting with error bars and mean line
    plt.figure(figsize=(4, 4))
    plt.errorbar(x_labels, z_values, marker="s", yerr=errors, capsize=4, **blue_point_kwargs)
    plt.axhline(mean_z, linestyle='--', label=r'$\overline{z} = $' + f"{mean_z:.2f}", color=colors[1])
    #plt.xlabel('Research Papers')
    plt.ylabel('z', rotation=0, ha="right")
    #plt.title('Results for Scalar Variable z for Different Experiments with Errors')
    plt.xticks(rotation=45, ha='center')  # Rotate x-axis labels for better readability
    plt.legend()
    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    plt.savefig("z-values.png", format="png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()