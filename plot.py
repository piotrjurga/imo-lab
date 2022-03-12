#!/usr/bin/env python
#%%
import numpy as np
import matplotlib.pyplot as plt

instances = ["kroA100", "kroB100"]
methods = ["greedy_simple", "greedy_loop", "regret_loop"]
names_map = {"greedy_simple": "Metoda najbliższego sąsiada", "greedy_loop": "Metoda rozbudowy cyklu", "regret_loop": "2-żal"}

for instance in instances:
    pos = np.fromfile(f'results/{instance}/pos.dat', dtype=np.int32)
    pos = pos.reshape([-1, 2])
    for method in methods: 
        a = np.fromfile(f'results/{instance}/{method}-a.dat', dtype=np.int32)
        b = np.fromfile(f'results/{instance}/{method}-b.dat', dtype=np.int32)

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
        ax = plt.gca()
        ax.axis('off')
        plt.savefig(f'results/instance-{instance}-method-{method}', dpi=300, pad_inches=0, bbox_inches='tight')
        plt.show()
        plt.clf()
