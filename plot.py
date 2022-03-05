#!/usr/bin/env python
#%%
import numpy as np
import matplotlib.pyplot as plt

pos = np.fromfile(f'results/pos.dat', dtype=np.int32)
pos = pos.reshape([-1, 2])
for i in range(1, 4): 
    a = np.fromfile(f'results/{i}a.dat', dtype=np.int32)
    b = np.fromfile(f'results/{i}b.dat', dtype=np.int32)

    pa = pos[a]
    pb = pos[b]
    pa = np.append(pa, [pa[0]], axis=0)
    pb = np.append(pb, [pb[0]], axis=0)
    xa, ya = zip(*pa)
    xb, yb = zip(*pb)
    plt.plot(xa, ya)
    plt.plot(xb, yb)
    plt.scatter(xa, ya)
    plt.scatter(xb, yb)
    plt.show()
