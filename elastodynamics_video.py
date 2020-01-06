import numpy as np
from matplotlib import pyplot as plt

n = 100
length = 1
k = 100
rho = 1
dt = 0.0001
distribution = 0.5

timesteps = int(4e4)
n_plot = int(5e1)

p = np.zeros([n+1, 2])
v = np.zeros([n-1, 2])
p[:, 0] = np.linspace(0, length, n+1)
p[:int(n/2), 1] = distribution*np.arange(0, int(n/2))/int(n/2)
p[int(n/2):, 1] = distribution*np.arange(int(n/2), -1, -1)/int(n/2)

# plt.plot(p[:, 0], p[:, 1], label='t=0s')

l0 = length/n
mass = rho*l0
ya = np.zeros(timesteps*5)
yb = np.zeros(timesteps*5)
t = np.arange(timesteps*5)*dt

for i in range(timesteps*5):
    d1 = p[2:, :] - p[1:-1, :]
    r1 = np.linalg.norm(d1, axis=1)
    # r1[r1<length/n] = 1
    r1 = r1.reshape([-1, 1])
    f1 = k*d1*(r1 - l0)/r1

    d2 = p[:-2, :] - p[1:-1, :]
    r2 = np.linalg.norm(d2, axis=1)
    # r2[r2<length/n] = 1
    r2 = r2.reshape([-1, 1])
    f2 = k*d2*(r2 - l0)/r2

    a = (f1 + f2)/mass
    v += a*dt
    p[1:-1, :] += v*dt

    ya[i] = p[50, 1]
    yb[i] = p[25, 1]
    if i % n_plot == 0:
        plt.plot(p[:, 0], p[:, 1])
        plt.ylim([-0.45, 0.6])
        plt.xlim([-0.05, 1.05])
        plt.xlabel('x/m')
        plt.ylabel('y/m')
        plt.savefig(f't1/fig_{i}.png')
        plt.close()
        print(i)
