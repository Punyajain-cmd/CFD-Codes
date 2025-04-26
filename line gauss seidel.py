import numpy as np
import matplotlib.pyplot as plt

t = np.zeros((41,21))
t1 = np.empty((41,21))
x = 1/(21-1)
y = 2/(41-1)
beta = x/y
k = beta*beta
g = (-2)*(1+k)

for j in range(41):
    for i in range(21):
        t[j,i]= 30

for i in range(21):
    t1[0,i]= 100
    t1[40,i]= 30
for j in range(40):
    t1[j+1,0]= 30
    t1[j+1,20]= 30

t6 = np.copy(t1)

B = np.zeros((19,19))
for i in range(19):
    for j in range(19):
        if j==i:
            B[i,j] = g
        elif abs(j-i) == 1:
            B[i,j] = 1

# det = np.linalg.det(B)
# print(B)

x = np.zeros(19)
c = np.zeros(19)
a=0
b=0
while a>(-1):
    for j in range(39):
        for i in range(19):
            if i==0 or i==18:
                c[i] = (-k)*(t[j+2,i+1] + t1[j,i+1]) - 30
            else :
                c[i] = (-k)*(t[j+2,i+1] + t1[j,i+1])
        x = np.dot(np.linalg.inv(B),(np.transpose(c)))

        for i in range(19):
            t1[j+1,i+1] = x[i]

    for i in range(41):
        for j in range(21):
            a = a + abs(t1[i,j] - t[i,j])
    b = b + 1
    if (a<0.01):
        break
    else :
        a=0
        t = np.copy(t1)
        t1 = np.copy(t6)

print('no. of iterration taken by line gauss seidel method is: ',b)

# heat flux per unit second function 
q=0
for i in range(40):
   for j in range (21):
        q += -50*(t1[i+1,j]-t1[i,j])/y

print(q)

plt.figure(figsize=(5,10))
heatmap = plt.pcolormesh(np.arange(0,21,1), np.arange(0,41,1), t1, shading='auto', cmap='jet')
plt.colorbar(heatmap, label='Temperature')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Count of temperature of 2D plot by line gauss seidel method')
plt.show()
