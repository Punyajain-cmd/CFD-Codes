import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx = 200
dx = 1.0 / nx
gamma = 1.4
nt = 200
dt = 0.0005

# Initial Conditions
rho = np.ones(nx)
rho[:nx//2] = 1.0
rho[nx//2:] = 0.125

u = np.zeros(nx)
p = np.ones(nx)
p[:nx//2] = 1.0
p[nx//2:] = 0.1

# Conserved Variables
E = p/(gamma-1) + 0.5*rho*u**2

for n in range(nt):
    F1 = rho*u
    F2 = rho*u**2 + p
    F3 = u*(E + p)
    
    # Forward Euler Time Integration
    rho[1:-1] -= dt/dx * (F1[2:] - F1[1:-1])
    u[1:-1] -= dt/dx * ((F2[2:] - F2[1:-1])/rho[1:-1])
    E[1:-1] -= dt/dx * (F3[2:] - F3[1:-1])
    p = (gamma-1)*(E - 0.5*rho*u**2)

# Plot
x = np.linspace(0, 1, nx)
plt.figure(figsize=(12,8))
plt.subplot(311)
plt.plot(x, rho)
plt.title('Density')
plt.grid()

plt.subplot(312)
plt.plot(x, u)
plt.title('Velocity')
plt.grid()

plt.subplot(313)
plt.plot(x, p)
plt.title('Pressure')
plt.grid()

plt.tight_layout()
plt.show()
