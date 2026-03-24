import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


# -----------------------------
# Load JSON
# -----------------------------
with open("../out/fa-minimizer/min-seeds-output.txt", "r") as f:
    full_data = json.load(f)

target_k = 21
target_w = 11

# Filter only k=21, w=11
data = [d for d in full_data
        if d["k"] == target_k and d["w"] == target_w]

# Extract unique parameters
b_values = sorted(set(d["b"] for d in data))
n_values = sorted(set(d["n"] for d in data))

# -----------------------------
# Styling
# -----------------------------

# Color per b
base_colors = ["#5A6F8E", "#6BAF92", "#D8C99B",
               "#C97C5D", "#8E5A99", "#4C956C"]

color_map = {
    b: base_colors[i % len(base_colors)]
    for i, b in enumerate(b_values)
}

# Line style per n
linestyle_map = {
    3: "-",
    5: "--",
    7: ":"
}

r = np.arange(6)  # r = 0..5

# --- Color legend (b values) ---
color_handles = [
    mpatches.Patch(color=color_map[b], label=f"b={b}")
    for b in b_values
]

# --- Line style legend (n values) ---
style_handles = [
    mlines.Line2D(
        [],
        [],
        color="black",
        linestyle=linestyle_map[n],
        linewidth=2,
        label=f"n={n}"
    )
    for n in n_values
]

# Combine
handles = color_handles + style_handles


# ==========================================================
# 1️⃣ UNIQUENESS PLOT
# ==========================================================

plt.figure(figsize=(10, 6))

plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13
})

for b in b_values:
    for n in n_values:

        entry = next((d for d in data
                      if d["b"] == b and d["n"] == n), None)

        if entry is None:
            continue

        uniqueness = np.array(
            [x["unique_seed_ratio"] for x in entry["details"]]
        )

        print(b, n, r, uniqueness)

        plt.plot(
            r,
            uniqueness,
            color=color_map[b],
            linestyle=linestyle_map.get(n, "-"),
            linewidth=2,
            label=f"b={b}, n={n}"
        )

plt.xlabel("Relaxation Radius (r)")
plt.ylabel("Unique Seed Ratio")
plt.ylim(0.4, 1.0)
plt.xticks(r)
plt.title(f"Uniqueness vs Relaxation (k={target_k}, w={target_w})")
plt.legend(
    handles=handles,
    loc="lower right",
    frameon=True,
    ncol=2,
    fontsize=13
)
plt.tight_layout()
plt.savefig("../out/uniqueness_k21_w11.png", dpi=300)
plt.close()

# ==========================================================
# 2️⃣ GAP PLOT
# ==========================================================

plt.figure(figsize=(10, 6))

for b in b_values:
    for n in n_values:

        entry = next((d for d in data
                      if d["b"] == b and d["n"] == n), None)

        if entry is None:
            continue

        gap = np.array(
            [x["gap_ratio"] for x in entry["details"]]
        )
        print(b, n, r, gap)

        plt.plot(
            r,
            gap,
            color=color_map[b],
            linestyle=linestyle_map.get(n, "-"),
            linewidth=2,
            label=f"b={b}, n={n}"
        )

plt.xlabel("Relaxation Radius (r)")
plt.ylabel("Gap Ratio")
plt.ylim(0.03, 0.08)
plt.xticks(r)
plt.title(f"Gap Ratio vs Relaxation (k={target_k}, w={target_w})")
plt.legend(
    handles=handles,
    loc="upper right",
    ncol=2,
    frameon=True,
    fontsize=13
)
plt.tight_layout()
plt.savefig("../out/figures/minimizers_gap_k21_w11.png", dpi=300)
plt.close()