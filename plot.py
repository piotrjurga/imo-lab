#!/usr/bin/env python
#%%
import numpy as np
import matplotlib.pyplot as plt
import itertools

if 1:
    # dumb debug code
    pos = np.fromfile(f'results/pos.dat', dtype=np.int32)
    pos = pos.reshape([-1, 2])
    a = np.fromfile('results/a.dat', dtype=np.int32)
    b = np.fromfile('results/b.dat', dtype=np.int32)
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
    exit()

instances = ["kroA100", "kroB100"]
methods = ["greedy_simple", "greedy_loop", "regret_loop"]
initial_solutions = ["random_loop", "greedy_loop"]
methods = ["random_walk", "greedy-neighbour_search_node", "greedy-neighbour_search_edge", "best-neighbour_search_node", "best-neighbour_search_edge"]

solutions = list(itertools.product([initial_solutions[0]], methods)) + list(itertools.product([initial_solutions[1]], methods))

for instance in instances:
    pos = np.fromfile(f'results/{instance}/pos.dat', dtype=np.int32)
    pos = pos.reshape([-1, 2])
    for solution in solutions:
        solution = "-".join(list(solution))
        a = np.fromfile(f'results/{instance}/{solution}-a.dat', dtype=np.int32)
        b = np.fromfile(f'results/{instance}/{solution}-b.dat', dtype=np.int32)

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
        plt.savefig(f'results/{instance}-{solution}', dpi=300, pad_inches=0, bbox_inches='tight')        
        plt.show()
        plt.clf()
