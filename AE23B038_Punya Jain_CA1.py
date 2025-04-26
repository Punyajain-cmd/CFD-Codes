import numpy as np
import matplotlib.pyplot as plt

T=np.zeros(11)
T1=np.zeros(11)
T1[0]=1.00
T1[10]=0.00
del_x=0.1; del_t=0.001; alpha=1.0; gamma= alpha*del_t/(del_x*del_x)
def f(n):
    for J in range (0, int(n/del_t)):
        for i in range (1,10):
            T1[i] = gamma*(T[i-1]) + (1-2*gamma)*(T[i]) + gamma*(T[i+1])
    
        for k in range(0,10):
            T[k]=T1[k]
    return T

plt.plot(np.arange(0,11), f(0), label="for t=0s", color='green')
plt.plot(np.arange(0,11), f(0.1), label="for t=0.1s", color='blue' ) 
plt.plot(np.arange(0,11), f(0.5), label="for t=0.5s", color='orange')
plt.plot(np.arange(0,11), f(1.0), label="for t=1.0s", color='red')
plt.plot(np.arange(0,11), f(2.0), label="for t=2.0s", color='purple')
plt.plot(np.arange(0,11), f(10.0), label="for t=10.0s", color='cyan')
plt.legend()
plt.xlabel('grid point (*1o^-1m)')
plt.ylabel('temperature')
plt.show()