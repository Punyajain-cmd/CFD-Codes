import numpy as np
import matplotlib.pyplot as plt
import tqdm

# defining grid length
dx = 1/(51-1); dy = 1/(51-1); dt = 0.01; rho = 1; Re = 100; k = (dx/dy)**2
# define u, v, P, u', v', P', u1, v1, P1
u = np.zeros((51,51))
v = np.zeros((51,51))
P = np.zeros((51,51))

u_p = np.zeros((51,51))
v_p = np.zeros((51,51))

u1 = np.zeros((51,51))
v1 = np.zeros((51,51))
P1 = np.zeros((51,51))

# Boundary condition
u[-1,1:-1] = u_p[-1,1:-1] = u1[-1,1:-1] = 1
P[:,:] = 0

print(u)
# solving 
m = np.zeros((51,51))
n = np.zeros((51,51))
f = np.zeros((51,51))
f1 = np.zeros((51,51))
m_p = np.zeros((51,51))
n_p = np.zeros((51,51))
m1 = np.zeros((51,51))
n1 = np.zeros((51,51))
def pro_vel(m, n, m_p, n_p):
    def sq(a):
        return a*a 
    def d(a, b, c):
        return (a - b)/c 
    def d2(a, b, c, dp):
        return (a - 2*b + c)/dp
    m_p[1:-1,1:-1] = m_p[1:-1,1:-1] - dt*m[1:-1,1:-1]*(d((m[1:-1,2:]), (m[1:-1,:-2]), 2*dx) + n[1:-1,1:-1]*d(m[2:,1:-1], m[:-2,1:-1], 2*dy)) + dt*(d2(m[1:-1,2:], m[1:-1,1:-1], m[1:-1,:-2], (dx**2)))/Re + dt*(d2(m[2:,1:-1], m[1:-1,1:-1], m[:-2,1:-1], (dy**2)))/Re
    n_p[1:-1,1:-1] = n_p[1:-1,1:-1] - dt*m[1:-1,1:-1]*(d(n[1:-1,2:], n[1:-1,:-2], 2*dx) + n[1:-1,1:-1]*d((n[2:,1:-1]), (n[:-2,1:-1]), 2*dy)) + dt*(d2(n[1:-1,2:], n[1:-1,1:-1], n[1:-1,:-2], (dx**2)))/Re + dt*(d2(n[2:,1:-1], n[1:-1,1:-1], n[:-2,1:-1], (dy**2)))/Re
    
    return m_p, n_p

def pressure(m, n, f, f1, a):
    def d(a, b, c):
        return (a - b)/c 
    f1[1:-1,1:-1] = f[1:-1,:-2] + f[1:-1,2:] + k*(f[:-2,1:-1] + f[2:,1:-1]) - (dx**2)*(d( m[1:-1,2:],m[1:-1,:-2], 2*dx)+ d(n[2:,1:-1], n[:-2,1:-1], 2*dy))/dt
    f1[1:-1,1:-1] = f1[1:-1,1:-1]/(2*(1+k))
    
    f1[:,0] = f1[:,1]
    f1[50,:] = f1[49,:]
    f1[0,:] = f1[1,:]
    f1[:,50] = 0
    return f1

def correction(m, n, f, m1, n1):
    def d(a, b, c):
        return (a - b)/c 
    m1[1:-1,1:-1] = m[1:-1,1:-1] - dt*(d(f[1:-1,2:], f[1:-1,:-2], 2*dx))
    n1[1:-1,1:-1] = n[1:-1,1:-1] - dt*(d(f[2:,1:-1], f[:-2,1:-1], 2*dy))
    return m1, n1

b=0
for _ in tqdm.tqdm(range(200)):
    (u_p, v_p) = pro_vel(u, v, u_p, v_p)
    for _ in range(50):
        P1 = pressure(u_p, v_p, P, P1)
        P = np.copy(P1)
    (u1, v1) = correction(u_p, v_p, P1, u1, v1)
    
    residual_u = np.max(np.abs(u1 - u))
    residual_v = np.max(np.abs(v1 - v))
    residual_P = np.max(np.abs(P1 - P))
    b += 1
    # if residual_u < 1e-6 and residual_v < 1e-6 and residual_P < 1e-6:
    #     print(f"Residuals: u={residual_u}, v={residual_v}, P={residual_P}")
    #     break
    # elif (b==1000):
    #     break
    u = np.copy(u1)
    v = np.copy(v1)
    P = np.copy(P1)

print(b)
# print(u1)
# print(v1)
# print(P1)

# Define grid points for plotting
x = np.linspace(0, 1, 51)  # x-coordinates (corresponding to dx)
y = np.linspace(0, 1, 51)   # y-coordinates (corresponding to dy)

# Create meshgrid for streamline plotting
X, Y = np.meshgrid(x, y)

u1p = np.zeros((51,51))
v1p = np.zeros((51,51))
u1p[:,:] = u1[:,:]
v1p[:,:] = v1[:,:]
# for i in range(51):
#     u1p[i,:] = u1[50-i,:]
#     v1p[i,:] = v1[50-i,:]
# Calculate velocity magnitude
velocity_magnitude = np.sqrt(u1p[:,:]**2 + v1p[:,:]**2)

# Streamline plot velocity_magnitude.T, cmap='viridis'
plt.figure(figsize=(6, 6))
plt.streamplot(X, Y,u1p, v1p, density=2, color='black')
plt.title('Streamline Plot')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid()
plt.show()

# u velocity 
plt.plot(u1[:,25], np.linspace(0, 1, 51), color = 'blue')
plt.show()

# v velocity 
plt.plot(np.linspace(0, 1, 51), u1[25,:], color = 'blue')
plt.show()

# plt.figure(figsize=(5,10))
# pressure = plt.pcolormesh(np.arange(0,51,1), np.arange(51,0,-1), u1p, shading='auto', cmap='jet')
# plt.colorbar(pressure, label='Temperature')
# plt.xlabel('X-axis')
# plt.ylabel('Y-axis')
# plt.title('Count of temperature of 2D plot')
# plt.show()