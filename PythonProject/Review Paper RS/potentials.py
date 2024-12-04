import numpy as np
import matplotlib.pyplot as plt

# Define the double-well potential function
def double_well_potential(x, a=1, b=1):
    return a * x**4 - b * x**2

# Generate x values
x = np.linspace(-1.3, 1.3, 500)

# Compute the corresponding y values
a = 1
b = -10
y = double_well_potential(x, a, b)

# Create the plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x, y, color='#ff0000ff', linewidth=2)

# Remove axes
ax.axis('off')

# Set equal aspect ratio for proper visualization
#ax.set_aspect('equal', adjustable='box')

# Save the plot as an SVG file
plt.savefig('double_well_potential.svg', format='svg', bbox_inches='tight')

# Show the plot (optional)
plt.show()


import numpy as np
import matplotlib.pyplot as plt

# Define the Landau-Ginzburg potential function
def landau_ginzburg_potential(phi, alpha=-1, beta=1, gamma=0.5):
    return -alpha * phi**2 + beta * phi**4 + gamma * phi**3

# Generate phi values
#phi = np.linspace(-2.5, 1.7, 500)
# phi = np.linspace(-2.5, 2.5, 500)
phi = np.linspace(-1.2, 1.2, 500)

# Compute the corresponding potential values
alpha = 1
beta = 1
gamma = 0
V = landau_ginzburg_potential(phi, alpha, beta, gamma)

# Create the plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(phi, V, color='#ff0000ff', linewidth=2)

# Remove axes
ax.axis('off')

# Save the plot as an SVG file
plt.savefig('landau_ginzburg_potential.svg', format='svg', bbox_inches='tight')

# Show the plot
plt.show()

