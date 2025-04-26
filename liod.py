import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

nx = 51        # no.of nodes in grid
L = 1.0        # Domain Size
dx = L/(nx-1)
dy = L/(nx-1)  # length of sides of one grid
b = 1000       # No. of Iteration 
dt = 0.001     # Time interval(delta t) in s
Re = 100       # Reynold's No.
rho = 1.0      # Density
a = 1.0        # lid's velocity

N_PRESSURE_POISSON_ITERATIONS = 50

x = np.linspace(0.0, L, nx)
y = np.linspace(0.0, L, nx)
X, Y = np.meshgrid(x, y) 

u = np.zeros_like(X)
v = np.zeros_like(X)
P = np.zeros_like(X)
u_p = np.zeros_like(X)
v_p = np.zeros_like(X)
u1 = np.zeros_like(X)
v1 = np.zeros_like(X)
P1 = np.zeros_like(X)
# defining discretisation function
# central difference scheme in x
def CDSx(f):
    f1 = np.zeros_like(f)
    f1[1:-1, 1:-1] = (f[1:-1, 2:] - f[1:-1, :-2]) / (2 * dx)
    return f1
# central difference scheme in y
def CDSy(f):
    f1 = np.zeros_like(f)
    f1[1:-1, 1:-1] = (f[2:, 1:-1] - f[:-2, 1:-1]) / (2 * dy)
    return f1

# for laplace equation 
def lap(f):
    f1 = np.zeros_like(f)
    f1[1:-1, 1:-1] = ((f[2:, 1:-1] - 2*f[1:-1, 1:-1] + f[:-2, 1:-1]) / dx**2 +
                      (f[1:-1, 2:] - 2*f[1:-1, 1:-1] + f[1:-1, :-2]) / dy**2)
    return f1

# main loop for ploting
for _ in tqdm(range(b)):
    du_dx = CDSx(u)
    dv_dx = CDSx(v)
    du_dy = CDSy(u)
    dv_dy = CDSy(v)
    d2u = lap(u)
    d2v = lap(v)

    # provisional velocity

    u_p = u - dt*(u * du_dx + v * du_dy) + dt * (1/Re) * d2u
    v_p = v - dt*(u * dv_dx + v * dv_dy) + dt * (1/Re) * d2v

    # Boundary condition for u_p and v_p
    u_p[0,:] = a
    u_p[-1,:] = 0.0
    u_p[:,-1] = 0.0
    u_p[:,0] = 0.0
    v_p[0,:] = 0
    v_p[-1,:] = 0.0
    v_p[:,-1] = 0.0
    v_p[:,0] = 0.0

    du_p_dx = CDSx(u_p)
    dv_p_dy = CDSy(v_p)
    cont = (du_p_dx + dv_p_dy) / dt

    for j in range(N_PRESSURE_POISSON_ITERATIONS):
        # pressure term 
        P1[1:-1, 1:-1] = (1/2*(dx**2 + dy**2)) * ((dy**2) * (P[1:-1, :-2] + P[1:-1, 2:]) + (dx**2) * (P[:-2, 1:-1] + P[2:, 1:-1]) - (dy * dx)**2 * cont[1:-1,1:-1])

        # Boundary condition for P1
        P1[:, -1] = P1[:, -2]
        P1[0, :] = P1[1, :]
        P1[:, 0] = P1[:, 1]
        P1[-1, :] = 0.0

        P = P1

    dP_dx = CDSx(P1)
    dP_dy = CDSy(P1)

    # true velocity by adding correction term as in terms of u_p, v_p and P1

    u1 = u_p - dt * dP_dx
    v1 = v_p - dt * dP_dy

    # Boundary condition for u_p and v_p
    u1[0,:] = a
    u1[-1,:] = 0.0
    u1[:,-1] = 0.0
    u1[:,0] = 0.0
    v1[0,:] = 0
    v1[-1,:] = 0.0
    v1[:,-1] = 0.0
    v1[:,0] = 0.0

    # as loop repeat we want to update old u,v,P with values of new u,v,P
    u = np.clip(u1, -1e6, 1e6)
    v = np.clip(v1, -1e6, 1e6)
    P = np.clip(P1, -1e6, 1e6)

# ploting 1: streamlines 
x = np.linspace(0.0, L, nx)
y = np.linspace(0.0, L, nx)
X, Y = np.meshgrid(x, y)

plt.figure()
plt.streamplot(X, Y, u1, v1, color="black")
plt.title("Streamlines")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()