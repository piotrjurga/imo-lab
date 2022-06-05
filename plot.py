#!/usr/bin/env python
# %%
import numpy as np
import matplotlib.pyplot as plt
import itertools

# %%
if 0:
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

#instances = ["kroA200", "kroB200"]
instances = ["kroA", "kroB"]
#methods = ["HEA-LS", "HEA+LS"]
methods = ["HEA+LS"]

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
        plt.savefig(f'results/{instance}-{method}', dpi=300, pad_inches=0, bbox_inches='tight')        
        plt.show()
        plt.clf()

# %%
scores = np.fromfile(f'results/scores.dat', dtype=np.int32)
edge_average = np.fromfile(f'results/edge_average.dat', dtype=np.float64)
vertex_average = np.fromfile(f'results/vertex_average.dat', dtype=np.float64)
edge_best = np.fromfile(f'results/edge_best.dat', dtype=np.float64)
vertex_best = np.fromfile(f'results/vertex_best.dat', dtype=np.float64)

print(vertex_average.mean(), vertex_average.min(), vertex_average.max())

edge_average /= 200
vertex_average /= 200
edge_best /= 200
vertex_best /= 200

minscore = scores.argmin()
best_mask = np.ones(scores.shape, bool)
best_mask[minscore] = False

fig, axs = plt.subplots(2, 2, figsize=(16,12))
axs = axs.ravel()
axs[0].scatter(scores, edge_average)
axs[1].scatter(scores, vertex_average)
axs[2].scatter(scores[best_mask], edge_best[best_mask])
axs[3].scatter(scores[best_mask], vertex_best[best_mask])

axs[0].set_title(f"Common edge average (r={np.corrcoef(scores, edge_average)[0][1]:.2f})")
axs[1].set_title(f"Common vertex average (r={np.corrcoef(scores, vertex_average)[0][1]:.2f})")
axs[2].set_title(f"Common edge best (r={np.corrcoef(scores[best_mask], edge_best[best_mask])[0][1]:.2f})")
axs[3].set_title(f"Common vertex best (r={np.corrcoef(scores[best_mask], vertex_best[best_mask])[0][1]:.2f})")
plt.setp(axs, ylim=(0, 1))
# plt.tight_layout()
fig.suptitle("kroB200", fontsize=20)
fig.savefig("kroB200-sim", facecolor='white')
fig.show()

# %%
