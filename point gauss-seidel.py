import numpy as np
import matplotlib.pyplot as plt

# Initialize the temperature arrays
t = np.zeros((41,21))
t1 = np.zeros((41,21))
x = 1/(21-1)       # delta x = l/(n_x-1)
y = 2/(41-1)    # delta x = l/(n_x-1)
beta = x/y
k = beta*beta
g = -2*(1+k)

# Set initial values for the temperature array
for j in range(41):
    for i in range(21):
        t[j,i] = 30

for i in range(21):
    t1[0,i] = 100  # Top boundary condition
    t1[40,i] = 30  # Bottom boundary condition
for j in range(40):
    t1[j + 1, 0] = 30  # Left boundary condition
    t1[j + 1, 20] = 30  # Right boundary condition

# Initialize iteration counters
a = 0  # Variable to track the convergence criterion
b = 0  # Variable to count the number of iterations

# Gauss-Seidel iteration loop for solving the 2D heat diffusion problem
while a>(-1):
    # Loop over all internal grid points to update the temperature using Gauss-Seidel method
    for i in range(1,40):  # Loop over rows
        for j in range(1,20):  # Loop over column
             # Update the temperature at each internal grid point based on neighboring values
            t1[i,j] = (t[i+1,j] + t1[i-1,j] + k*(t[i,j+1]) + k*(t1[i,j-1]))/(-g)
    
    # Check convergence by computing the sum of absolute differences between the old and new temperatures
    for i in range(41):
        for j in range(21):
            a = a + abs(t1[i,j] - t[i,j])
    b = b + 1  # Increment the iteration count

     # Break the loop if the convergence criterion is met
    if (a<0.01):
        break
    else :
        # Reset the convergence criterion and update the old temperature values with the new ones
        a=0
        for i in range(41):
            for j in range(21):
                t[i,j] = t1[i,j]

print('no. of iterration taken by Point gauss-seidel method is: ',b)

# heat flux per unit second function 
q=0
for i in range(40):
   for j in range (21):
        q += -50*(t1[i+1,j]-t1[i,j])/y

print(q)

# Plot the temperature distribution using a 2D color plot (heatmap)
plt.figure(figsize=(5,10))
heatmap = plt.pcolormesh(np.arange(0,21,1), np.arange(0,41,1), t1, shading='auto', cmap='jet')
plt.colorbar(heatmap, label='Temperature')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Count of temperature of 2D plot by Point gauss-seidel method')
plt.show()


    