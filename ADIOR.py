import numpy as np
import matplotlib.pyplot as plt

# Initialize the temperature arrays
t = np.zeros((41, 21))
t1 = np.empty((41, 21))
x = 1/(21 - 1)    # delta x = l/(n_x-1)
y = 2/(41 - 1)    # delta x = l/(n_x-1)
beta = x/y
k = beta*beta
g = -2*(1 + k)
omega=1.3     # relaxation factor
omega1 = 1 - omega

# Set initial values for the temperature array
for j in range(41):
    for i in range(21):
        t[j,i] = 30

for i in range(21):
    t1[0,i] = 100  # Top boundary condition
    t1[40,i] = 30  # Bottom boundary condition
for j in range(40):
    t1[j+1,0] = 30  # Left boundary condition
    t1[j+1,20] = 30  # Right boundary condition

# creating two array which are having same condition as t1
t1_2 = np.copy(t1)
t6 = np.copy(t1)

# Construct the B matrix for solving the linear system
B = np.zeros((19, 19))
for i in range(19):
    for j in range(19):
        if j == i:
            B[i,j] = g
        elif abs(j-i) == 1:
            B[i,j] = omega

# Construct another matrix B1 for the second part of the problem
B1 = np.zeros((39, 39))
for i in range(39):
    for j in range(39):
        if j == i:
            B1[i,j] = g
        elif abs(j - i) == 1:
            B1[i,j] = omega*k

# Initialize variables for the iterative process
x = np.zeros(19)
x1 = np.zeros(39)
c = np.zeros(19)
c1 = np.zeros(39)
a = float('inf')  # Large initial value to enter the while loop
b = 0  # Iteration counter

# Iterative loop for solving the equations
while a > 0.01:  # Check convergence criterion
    a = 0  # Reset convergence parameter at the start of each iteration
    
    for j in range(39):
        for i in range(19):
            if i == 0 or i==18:
                c[i] = (-1)*(omega1)*((-g)*(t[j+1,i+1])) - omega*k*(t[j+2,i+1]+t1_2[j,i+1]) - (30 * omega)
            else:
                c[i] = (-1)*(omega1)*((-g)*(t[j+1,i+1])) - omega*k*(t[j+2,i+1]+t1_2[j,i+1])
        
        # Solve the system of equations using a stable solver
        x = np.linalg.solve(B, np.transpose(c))

        # Update the values in t1_2 based on the solution x
        for i in range(19):
            t1_2[j+1,i+1] = x[i]
    
    for i in range(19):
        for j in range(39):
            if j == 0 :
                c1[j] = (-1)*(omega1)*((-g)*(t1_2[j+1,i+1])) - omega*(t1_2[j+1,i+2]+t1[j+1,i]) - 100*k*omega
            elif j == 38:
                c1[j] = (-1)*(omega1)*((-g)*(t1_2[j+1,i+1])) - omega*(t1_2[j+1,i+2]+t1[j+1,i]) - 30*k*omega
            else:
                c1[j] = (-1)*(omega1)*((-g)*(t1_2[j+1,i+1])) - omega*(t1_2[j+1,i+2]+t1[j+1,i])
        
        # Solve the second system of equations
        x1 = np.linalg.solve(B1, np.transpose(c1))

        # Update the values in t1 based on the solution x1
        for j in range(39):
            t1[j+1,i+1] = x1[j]

    # Check convergence by calculating the difference between t1 and t
    for i in range(41):
        for j in range(21):
            a += abs(t1[i,j] - t[i,j])

    b += 1  # Increment the iteration count

    # Break the loop if convergence criterion is met
    if a < 0.01:
        break
    else:
        t = np.copy(t1)
        t1 = np.copy(t6)
        t1_2 = np.copy(t6)  # Update t for the next iteration

print('Number of iterations for ADI method with relaxation parameter(omega):', omega, 'is:', b)

q=0
for i in range(40):
   for j in range (21):
        q += -50*(t1[i+1,j]-t1[i,j])/y

print(q) 


# Plot the 2D heatmap of the temperature distribution
plt.figure(figsize=(5, 10))
heatmap = plt.pcolormesh(np.arange(0, 2.1, 0.1), np.arange(0, 4.1, 0.1), t1, shading='auto', cmap='jet')
plt.colorbar(heatmap, label='Temperature')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Count of temperature of 2D plot')
plt.show()

