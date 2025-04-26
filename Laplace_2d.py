import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx, ny = 50, 50
lx, ly = 2.0, 2.0
dx, dy = lx/(nx-1), ly/(ny-1)
max_iter = 1000

# Initialize potential field
phi = np.zeros((ny, nx))

# Boundary conditions
phi[:, 0] = 0     # Left
phi[:, -1] = 0    # Right
phi[0, :] = 100   # Top
phi[-1, :] = 0    # Bottom

# Iterative solver (Gauss-Seidel)
for it in range(max_iter):
    phi_old = phi.copy()
    phi[1:-1,1:-1] = 0.25 * (phi_old[1:-1,2:] + phi_old[1:-1,0:-2] +
                              phi_old[2:,1:-1] + phi_old[0:-2,1:-1])

# Plot
x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

plt.contourf(X, Y, phi, 50, cmap='viridis')
plt.colorbar(label='Potential φ')
plt.title('2D Laplace Equation Solution')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
