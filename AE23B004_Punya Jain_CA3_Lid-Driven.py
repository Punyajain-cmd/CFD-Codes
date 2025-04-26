import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

nx = 51        # no.of nodes in grid
L = 1.0        # Domain Size
dx = L/(nx-1)
dy = L/(nx-1)  # length of sides of one grid
b =  5000       # No. of Iteration 
dt = 0.001     # Time interval(delta t) in s
Re = 100       # Reynold's No.
rho = 1.0      # Density
a = 1.0        # lid's velocity

N_PRESSURE_POISSON_ITERATIONS = 50

# Grid setup for plotting 
dx = L / (nx - 1)
x = np.linspace(0.0, L, nx)
y = np.linspace(0.0, L, nx)
X, Y = np.meshgrid(x, y)

# Initialize velocity and pressure fields
u = np.zeros((nx, nx))  # Horizontal velocity
v = np.zeros((nx, nx))  # Vertical velocity
P = np.zeros((nx, nx))  # Pressure field

# Helper Functions
# Compute the central difference in the x-direction.
def CDSx(f):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (f[1:-1, 2:] - f[1:-1, :-2]) / (2 * dx)
    return diff

# Compute the central difference in the y-direction.
def CDSy(f):
    d = np.zeros_like(f)
    d[1:-1, 1:-1] = (f[2:, 1:-1] - f[:-2, 1:-1]) / (2 * dx)
    return d

# Compute the Laplacian of a field.
def lap(f):
    f1 = np.zeros_like(f)
    f1[1:-1, 1:-1] = ( f[1:-1, :-2] + f[:-2, 1:-1] + f[1:-1, 2:] + f[2:, 1:-1] - 4 * f[1:-1, 1:-1] ) / (dx**2)
    return f1

# # Stability Check
# max_dt = 0.5 * dx**2 / Re
# if dt > SAFETY_FACTOR * max_dt:
#     raise ValueError("Time step is too large for stability. Reduce dt.")

# Main loop for calculation
for _ in tqdm(range(b), desc="Solving Navier-Stokes"):
    # Compute derivatives
    du_dx = CDSx(u)
    du_dy = CDSy(u)
    dv_dx = CDSx(v)
    dv_dy = CDSy(v)
    lap_u = lap(u)
    lap_v = lap(v)

    # provisional velocity fields
    u_p = u + dt * (-u * du_dx - v * du_dy + 1/Re * lap_u)
    v_p = v + dt * (-u * dv_dx - v * dv_dy + 1/Re * lap_v)

    # Apply boundary conditions for tentative velocities
    u_p[0, :] = 0.0
    u_p[:, 0] = u_p[:, -1] = 0.0
    u_p[-1, :] = a
    v_p[0, :] = v_p[-1, :] = 0.0
    v_p[:, 0] = v_p[:, -1] = 0.0

    # Compute divergence of Provisional velocities
    divergence = rho / dt * ( CDSx(u_p) + CDSy(v_p) )

    # Solve for pressure using iterative Poisson equation
    for _ in range(N_PRESSURE_POISSON_ITERATIONS):
        P1 = np.zeros_like(P)
        P1[1:-1, 1:-1] = 0.25 * ( P[1:-1, :-2] + P[:-2, 1:-1] + P[1:-1, 2:] + P[2:, 1:-1] - dx**2 * divergence[1:-1, 1:-1] )
        # Pressure boundary conditions
        P1[:, -1] = P1[:, -2]
        P1[0, :] = P1[1, :]
        P1[:, 0] = P1[:, 1]
        P1[-1, :] = 0.0
        P = P1

    # Correct velocities using pressure gradient
    dp_dx = CDSx(P)
    dp_dy = CDSy(P)
    u = u_p - dt / rho * dp_dx
    v = v_p - dt / rho * dp_dy

    # Apply boundary conditions for corrected velocities
    u[0, :] = 0.0
    u[:, 0] = u[:, -1] = 0.0
    u[-1, :] = a
    v[0, :] = v[-1, :] = 0.0
    v[:, 0] = v[:, -1] = 0.0

# Post-Processing for finding vorticiy w
w = CDSx(v) - CDSy(u)

# Plot1 Streamlines
plt.figure()
plt.streamplot(X, Y, u, v, color="black")
plt.title("Streamlines")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# Plot2 Vorticity Contours
plt.figure()
plt.contourf(X, Y, w, cmap="coolwarm")
plt.colorbar(label="Vorticity")
plt.title("Vorticity Contours")
plt.xlabel("X")
plt.ylabel("Y") 
plt.show()

# Plot3 Pressure Distribution
plt.figure()
plt.contourf(X, Y, P, cmap="viridis")
plt.colorbar(label="Pressure")
plt.title("Pressure Distribution")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# Plot4 U-Velocity Profile
x_a = nx // 2
plt.figure()
plt.plot(u[:, x_a], y, label="U-Velocity")
plt.title("U-Velocity vs Y at X=0.5")
plt.xlabel("U-Velocity")
plt.ylabel("Y")
plt.legend()
plt.grid()
plt.show()

# Plot4 V-Velocity Profile
y_a = nx // 2
plt.figure()
plt.plot(x, v[y_a, :], label="V-Velocity")
plt.title("V-Velocity vs X at Y=0.5")
plt.xlabel("X")
plt.ylabel("V-Velocity")
plt.legend()
plt.grid()
plt.show()
