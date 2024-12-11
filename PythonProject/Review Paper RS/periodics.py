import numpy as np
import matplotlib.pyplot as plt

# Generate x values for 4 periods of a cosine
T = 2 * np.pi  # Period of the cosine
x = np.linspace(0, 1 * T, 1000)  # Cover 4 periods

# Compute the cosine values
y = np.cos(x)

# Create the plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x, y, color='#ff0000ff', linewidth=2)

# Remove axes
ax.axis('off')

plt.tight_layout()

# Save the plot as an SVG file
plt.savefig('cosine.svg', format='svg', bbox_inches='tight', transparent=True)

# Show the plot
plt.show()


# Generate x values
T = 2 * np.pi  # Period of the first cosine
# 12.6 * T for long periodics
# 4.1 for short
x = np.linspace(0.125, 4.1* T, 1000)  # Cover 4 periods of the first cosine

# Compute the superposition of two cosines
A1, A2 = 1, 0.3  # Amplitudes of the two cosines
A3 = 0.08  # for this linear decreasation of thingy
y = A1 * np.cos(x) + A2 * np.sin(x / 2) - A3 * x  # Cosine with periods T and 2T

# Create the plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(x, y, color='#1e4375ff', linewidth=4)

# Remove axes
ax.axis('off')

plt.tight_layout()
# Save the plot as an SVG file
plt.savefig('superposition_cosines.svg', format='svg', bbox_inches='tight', transparent=True)

# Show the plot
plt.show()
