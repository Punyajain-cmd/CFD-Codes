import numpy as np
import matplotlib.pyplot as plt

# Given
M1 = 2.0            # upstream Mach number
theta_deg = 10      # wedge angle (degrees)
gamma = 1.4         # specific heat ratio

theta = np.radians(theta_deg)

# Calculate shock angle using approximate solution (small theta)
beta = np.arcsin(1/M1)

# Post-shock properties
Mn1 = M1 * np.sin(beta)
Mn2 = np.sqrt((1 + (gamma-1)/2 * Mn1**2) / (gamma*Mn1**2 - (gamma-1)/2))
M2 = Mn2 / np.sin(beta - theta)

p2_p1 = 1 + 2*gamma/(gamma+1) * (Mn1**2 -1)

# Print Results
print(f"Shock angle (beta): {np.degrees(beta):.2f} degrees")
print(f"Post-shock Mach number M2: {M2:.2f}")
print(f"Pressure ratio p2/p1: {p2_p1:.2f}")
