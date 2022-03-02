#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

a = np.fromfile('results/a.dat', dtype=np.int32)
b = np.fromfile('results/b.dat', dtype=np.int32)
pos = np.fromfile('results/pos.dat', dtype=np.int32)
pos = pos.reshape([-1, 2])

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
