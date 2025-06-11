import numpy as np
import matplotlib.pyplot as plt
import os

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

path = "Simulation/Bifurkation"
files = os.listdir(path)
files.sort(key=lambda x: float(x.split("_")[1][:-4]))
As = np.array([float(f.split("_")[1][:-4]) for f in files])
OFFSET = 0

fig, ax = plt.subplots(1, 1, figsize=(9.6, 4.8))

xs = []
ys = []

for i, file in enumerate(files):
    data = np.loadtxt(path + "/" + file, delimiter="\t")
    y = np.round(data, 6)
    ps = np.unique(y)
    xs.extend([As[i] for _ in range(len(ps))])
    ys.extend(ps)

xs = np.array(xs) * 2
ys = np.array(ys)

y = ys[np.mod(xs * 1e3, 1) < 1e-8]
x = xs[np.mod(xs * 1e3, 1) < 1e-8]

hist, _, _ = np.histogram2d(x, y, bins=(800, 800))
hist = hist.T
hist = np.log(hist + OFFSET)
# edges = np.abs(sobel(hist, axis=0)) + 1
print(np.min(hist), np.min(hist[hist > 0]), np.max(hist))

ax.imshow(
    hist,
    vmin=-0.5,
    vmax=3,
    origin="lower",
    cmap="inferno",
    extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
)
ax.set_aspect("auto")

########################################
axins = inset_axes(
    ax,
    "30%",
    "65%",
    bbox_to_anchor=[0.01, 0.01, 1, 1],
    bbox_transform=ax.transAxes,
    loc="lower left",
)
mark_inset(ax, axins, loc1=4, loc2=2, fc="none", ec="#ddd", ls="--", lw=1.5, alpha=0.5)

y = ys[xs > 5.5]
x = xs[xs > 5.5]
y = y[x < 7.0]
x = x[x < 7.0]
y = y[np.mod(x * 2e3, 1) < 1e-8]
x = x[np.mod(x * 2e3, 1) < 1e-8]

hist, _, _ = np.histogram2d(x, y, bins=(500, 500))
hist = hist.T
hist = np.log(hist + OFFSET)
# edges = np.abs(sobel(hist, axis=0)) + 1
print(np.min(hist), np.min(hist[hist > 0]), np.max(hist))

axins.imshow(
    hist,
    vmin=-0.5,
    vmax=3.5,
    origin="lower",
    cmap="inferno",
    extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
)
axins.set_aspect("auto")

########################################
axinsins = inset_axes(
    axins,
    "12%",
    "35%",
    bbox_to_anchor=[0.33, 0.12, 1, 1],
    bbox_transform=ax.transAxes,
    loc="lower left",
)
mark_inset(
    axins, axinsins, loc1=4, loc2=2, fc="none", ec="#ddd", ls="--", lw=1.5, alpha=0.5
)

y = ys[xs > 6.15]
x = xs[xs > 6.15]
y = y[x < 6.25]
x = x[x < 6.25]
x = x[y < -30]
y = y[y < -30]
y = y[np.mod(x * 1e4, 1) < 1e-8]
x = x[np.mod(x * 1e4, 1) < 1e-8]

hist, _, _ = np.histogram2d(x, y, bins=(200, 300))
hist = hist.T
hist = np.log(hist + OFFSET)
# edges = np.abs(sobel(hist, axis=0)) + 1
print(np.min(hist), np.min(hist[hist > 0]), np.max(hist))

axinsins.imshow(
    hist,
    vmin=-0.5,
    vmax=3.0,
    origin="lower",
    cmap="inferno",
    extent=[np.min(x), np.max(x), np.min(y), np.max(y)],
)
axinsins.set_aspect("auto")


ax.grid(False)
axins.grid(False)
axins.spines[:].set_color("#ddd")  # Matches grid
axins.spines[:].set_linewidth(1.5)
axins.set_xticks([])
axins.set_yticks([])
axinsins.grid(False)
axinsins.spines[:].set_linewidth(1.5)
axinsins.spines[:].set_color("#ddd")  # Matches grid
axinsins.set_xticks([])
axinsins.set_yticks([])
ax.set_xlabel("$A$ in $U_s$")
ax.set_ylabel("$i$ in $I_0$")
# plt.savefig("Figures/bifurkation_sim.png", dpi=300, bbox_inches="tight")
plt.show()
