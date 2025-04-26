import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx, ny = 64, 64
lx, ly = 2*np.pi, 2*np.pi
dx, dy = lx/nx, ly/ny
nu = 0.01
dt = 0.001
nt = 500

# Grid
x = np.linspace(0, lx, nx)
y = np.linspace(0, ly, ny)
X, Y = np.meshgrid(x, y)

# Initial conditions
u = np.sin(X) * np.cos(Y)
v = -np.cos(X) * np.sin(Y)

for n in range(nt):
    un = u.copy()
    vn = v.copy()
    u = un + nu*dt*(np.roll(un,1,axis=0) -2*un + np.roll(un,-1,axis=0))/dx**2 \
            + nu*dt*(np.roll(un,1,axis=1) -2*un + np.roll(un,-1,axis=1))/dy**2
    v = vn + nu*dt*(np.roll(vn,1,axis=0) -2*vn + np.roll(vn,-1,axis=0))/dx**2 \
            + nu*dt*(np.roll(vn,1,axis=1) -2*vn + np.roll(vn,-1,axis=1))/dy**2

# Plot
plt.streamplot(x, y, u, v, density=2)
plt.title('Taylor-Green Vortex at t = {:.2f}'.format(nt*dt))
plt.xlabel('x')
plt.ylabel('y')
plt.show()
