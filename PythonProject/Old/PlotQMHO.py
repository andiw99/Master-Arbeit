import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite
from scipy.special import factorial

# Parameters
n_max = 5  # Number of eigenfunctions to plot
x = np.linspace(-5, 5, 400)  # Range of x values

# Define the quantum harmonic oscillator eigenfunctions
def harmonic_eigenfunction(x, n):
    psi = np.exp(-x**2 / 2) * hermite(n)(x) / np.sqrt(2**n * factorial(n) * np.sqrt(np.pi))
    return psi

# Calculate eigenenergies
def eigenenergy(n):
    return n + 0.5

# Plot the eigenfunctions
plt.figure(figsize=(6, 6))
for n in range(n_max):
    eigenfunction = harmonic_eigenfunction(x, n)
    energy = eigenenergy(n)
    scaled_eigenfunction = eigenfunction * 0.6 + energy  # Scale and shift
    plt.plot(x, scaled_eigenfunction, label=f'n={n}')

# Add labels and legend
plt.title('Eigenfunctions of the Quantum Harmonic Oscillator')
plt.xlabel('Position (x)')
plt.ylabel('Wavefunction (ψ)')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig('quantum_harmonic_eigenfunctions_shifted.pdf')

# Show the plot
plt.show()
