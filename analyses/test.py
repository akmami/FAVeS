import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

# -------- CONFIG --------
INPUT_FILE = "../out/fa-blend/fuzzy-seeds-output.txt"
NUM_GRID_POINTS = 50  # smoothness of the surface
# ------------------------

# Regex patterns (anchored to avoid collisions)
kw_pattern = re.compile(r"^k:\s*(\d+),\s*w:\s*(\d+)")
uniq_pattern = re.compile(r"^Uniqueness %:\s*([0-9.]+)")
runiq_pattern = re.compile(r"^Relaxed-Uniqueness %:\s*([0-9.]+)")


data = []

with open(INPUT_FILE, "r") as f:
    lines = f.readlines()

k = w = uniq = runiq = None

for line in lines:
    kw_match = kw_pattern.search(line)
    if kw_match:
        k = int(kw_match.group(1))
        w = int(kw_match.group(2))

    uniq_match = uniq_pattern.search(line)
    if uniq_match:
        uniq = float(uniq_match.group(1))

    runiq_match = runiq_pattern.search(line)
    if runiq_match:
        runiq = float(runiq_match.group(1))

    # Once we have a full record, store it
    if k is not None and w is not None and uniq is not None and runiq is not None:
        data.append((k, w, uniq, runiq))
        k = w = uniq = runiq = None

# Convert to numpy arrays
data = np.array(data)
K = data[:, 0]
W = data[:, 1]
U = data[:, 2]
RU = data[:, 3]


# Create a fine grid for smooth interpolation
k_lin = np.linspace(K.min(), K.max(), NUM_GRID_POINTS)
w_lin = np.linspace(W.min(), W.max(), NUM_GRID_POINTS)
K_grid, W_grid = np.meshgrid(k_lin, w_lin)

# Interpolate Uniqueness
U_grid = griddata((K, W), U, (K_grid, W_grid), method='cubic')
RU_grid = griddata((K, W), RU, (K_grid, W_grid), method='cubic')

# -------- PLOT 1: Uniqueness --------
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")
surf1 = ax.plot_surface(K_grid, W_grid, U_grid, cmap="viridis", edgecolor="k", alpha=0.9)
ax.set_title("Uniqueness (%) vs k and w")
ax.set_xlabel("k")
ax.set_ylabel("w")
ax.set_zlabel("Uniqueness (%)")
fig.colorbar(surf1, shrink=0.6, aspect=10)
plt.tight_layout()
plt.show()

# -------- PLOT 2: Relaxed Uniqueness --------
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")
surf2 = ax.plot_surface(K_grid, W_grid, RU_grid, cmap="plasma", edgecolor="k", alpha=0.9)
ax.set_title("Relaxed Uniqueness (%) vs k and w")
ax.set_xlabel("k")
ax.set_ylabel("w")
ax.set_zlabel("Relaxed Uniqueness (%)")
fig.colorbar(surf2, shrink=0.6, aspect=10)
plt.tight_layout()
plt.show()