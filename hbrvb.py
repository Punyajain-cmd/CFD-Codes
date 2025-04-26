import numpy as np
import matplotlib.pyplot as plt

# defining grid length
dx = 10/(51-1); dy = 1/(51-1); dt = 0.01; rho = 1; Re = 100; k = (dx/dy)**2
# define u, v, P, u', v', P', u1, v1, P1
u = np.zeros((51,52))
v = np.zeros((51,52))
P = np.zeros((51,52))

u_p = np.zeros((51,52))
v_p = np.zeros((51,52))

u1 = np.zeros((51,52))
v1 = np.zeros((51,52))
P1 = np.zeros((51,52))

# Boundary condition
u[:,0] = 1
u[0, 1:] = u[50, 1:] = 0
v[:, :] = 0
u1[:,0] = 1
u_p[:,0] = 1
v[0,:] = v[50,:] = 0
P1[:,51] = 1 # atm
P[:,:] = 1 # atm
P[:,0] = 0


# solving 
m = np.zeros((51,52))
n = np.zeros((51,52))
f = np.zeros((51,52))
f1 = np.zeros((51,52))
m_p = np.zeros((51,52))
n_p = np.zeros((51,52))
m1 = np.zeros((51,52))
n1 = np.zeros((51,52))
def pro_vel(m, n, m_p, n_p):
    def sq(a):
        return a*a 
    def d(a, b, c):
        return (a - b)/c 
    def d2(a, b, c, dp):
        return (a - 2*b + c)/dp
    m_p[1:49,1:50] = m_p[1:49,1:50] - dt*(d(sq(m[1:49,2:51]), sq(m[1:49,0:49]), 2*dx) + d(m[2:50,1:50]*n[2:50,1:50], m[0:48,1:50]*n[0:48,1:50], 2*dy)) + dt*(d2(m[1:49,2:51], m[1:49,1:50], m[1:49,0:49], 2*(dx**2)))/Re + dt*(d2(m[2:50,1:50], m[1:49,1:50], m[0:48,1:50], 2*(dy**2)))/Re
    # n_p[1:49,1:50] = n_p[1:49,1:50] - dt*(d(m[1:49,2:51]*n[1:49,2:51], m[1:49,0:49]*n[1:49,0:49], 2*dx) + d(sq(n[2:50,1:50]), sq(n[0:48,1:50]), 2*dy)) + dt*(d2(n[1:49,2:51], n[1:49,1:50], n[1:49,0:49], 2*(dx**2)))/Re + dt*(d2(n[2:50,1:50], n[1:49,1:50], n[0:48,1:50], 2*(dy**2)))/Re
    m_p[1:49,50] = m_p[1:49,49]
    m_p[1:49,51] = m_p[1:49,50]

    return m_p, n_p

def pressure(m, n, f, f1):
    def d(a, b, c):
        return (a - b)/c 
    f1[1:49,1:50] = f[1:49,0:49] + f[1:49,2:51] + k*(f[0:48,1:50] + f[2:50,1:50]) - (dx**2)*(d( m[1:49,2:51],m[1:49,0:49], 2*dx)+ d(n[2:50,1:50], n[0:48,1:50], 2*dy))/dt
    f1[1:49,1:50] = f1[1:49,1:50]/(2*(1+k))
    
    f1[:,0] = f1[:,1]
    f1[50,:] = f1[49,:]
    f1[0,:] = f1[1,:]
    f1[:,51] = f1[:,50]
    return f1

def correction(m, n, f, m1, n1):
    def d(a, b, c):
        return (a - b)/c 
    m1[1:50,1:51] = m[1:50,1:51] - dt*(d(f[1:50,1:51], f[1:50,0:50], 2*dx))
    # for j in range(1,51):
    #     n1[j,:] = n[j,:] - dt*(d(f[j,:], f[j-1,:], dy))
    m1[:,51] = m1[:,50]
    n1[:,0] = n1[:,1]
    return m1, n1

a=0
b=0
c=0
while a==0:
    (u_p, v_p) = pro_vel(u, v, u_p, v_p)
    P1 = pressure(u_p, v_p, P, P1)
    (u1, v1) = correction(u_p, v_p, P1, u1, v1)
    
    for i in range(52):
        for j in range(51):
            a +=  abs(u1[j,i] - u[j,i])
    # for i in range(1,51):
    #     for j in range(51):
    #         c +=  abs(u1[j,i+1] - u1[j,i])/dx
    # a = np.linalg.norm(u1 - u) / np.linalg.norm(u1)
    b += 1
    if (a<1e-7):
        break
    else :
        a = 0
        u = np.copy(u1)
        v = np.copy(v1)
        P = np.copy(P1)

print(b)
# print(u1)
# print(v1)
# print(P1)

for i in range(1,50):
    c=0
    for j in range(51):
        c +=  abs(u1[j,i+1] - u1[j,i])/(dx)
    # print(c)
    if (c<1):
        print(i)
        break

        

plt.plot(u1[:,40], np.arange(50,-1,-1), color = 'blue')
plt.show()

# Define grid points for plotting
x = np.linspace(0, 10, 51)  # x-coordinates (corresponding to dx)
y = np.linspace(0, 1, 51)   # y-coordinates (corresponding to dy)

# Create meshgrid for streamline plotting
X, Y = np.meshgrid(x, y)

u1p = np.zeros((51,51))
v1p = np.zeros((51,51))
u1p[:,:] = u1[:,1:]
v1p[:,:] = v1[:,1:]
# Calculate velocity magnitude
velocity_magnitude = np.sqrt(u1p[:,1:]**2 + v1p[:,1:]**2)

# # Streamline plot
# plt.figure(figsize=(10, 4))
# plt.streamplot(X, Y,u1p, v1p, density=1.5, color=velocity_magnitude.T, cmap='viridis')
# plt.colorbar(label='Velocity Magnitude')
# plt.title('Streamline Plot')
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.grid()
# plt.show()

plt.figure(figsize=(5,10))
pressure = plt.pcolormesh(np.arange(0,51,1), np.arange(51,0,-1), u1p, shading='auto', cmap='jet')
plt.colorbar(pressure, label='Temperature')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Count of temperature of 2D plot')
plt.show()
    
   

