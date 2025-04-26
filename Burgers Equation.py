import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx = 101           # number of grid points
dx = 2.0 / (nx-1)  # grid spacing
nt = 50            # number of time steps
dt = 0.01          # time step size

# Initialize u
u = np.ones(nx)
u[int(0.5 / dx):int(1 / dx + 1)] = 2  # initial condition: u=2 between 0.5 and 1

# Time stepping loop
for n in range(nt):
    un = u.copy()
    u[1:] = un[1:] - un[1:] * dt / dx * (un[1:] - un[0:-1])

# Plot the results
x = np.linspace(0, 2, nx)
plt.plot(x, u)
plt.xlabel('x')
plt.ylabel('u')
plt.title('1D Inviscid Burgers Equation')
plt.grid(True)
plt.show()
