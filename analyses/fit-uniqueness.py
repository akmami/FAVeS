#!/usr/bin/env python3
"""
Fit and plot the null -> repeat-aware uniqueness model

For every genome it:
  * overlays the measured p(k) against the null    p_rand(k) = exp(-N Q^k),
  * fits the repeat-aware model                    p(k) = exp(-(1-a) N Q^k)*(1 - a (k/k0)^-b),
  * marks the crossover k_c = ln N / ln(1/Q),
  * plots the repeat spectrum  S_rep(k) = Pr(L >= k)  from the match-length histogram.
"""
import sys, json, math, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

PALETTE = ["#059669", "#2563eb", "#d97706", "#dc2626", "#7c3aed", "#0891b2", "#be185d"]


def load(path):
    d = json.load(open(path))
    ks = np.array([u["k"] for u in d["uniqueness"]], float)
    ps = np.array([u["p_unique"] for u in d["uniqueness"]], float)
    short = str(d.get("short", "")).strip()
    name = short if short else d.get("genome", path).split("/")[-1].split(".")[0]
    hist = d.get("match_length_hist", [])
    L = np.array([e["L"] for e in hist], float)
    Lc = np.array([e["count"] for e in hist], float)
    return dict(name=name, N=float(d["valid_len"]), Q=float(d["Q"]),
                kc=float(d["k_c"]), ks=ks, ps=ps, L=L, Lc=Lc,
                overflow=float(d.get("match_length_overflow", 0)))


def model(k, N, Q, a, b, k0):
    k = np.asarray(k, float)
    Srep = np.where(k < k0, 1.0, np.power(k / k0, -b))
    return np.exp(-(1.0 - a) * N * np.power(Q, k)) * (1.0 - a * Srep)


def fit(g, k0):
    if not HAVE_SCIPY:
        return None
    m = g["ks"] >= 1
    ks, ps = g["ks"][m], g["ps"][m]
    f = lambda k, a, b: model(k, g["N"], g["Q"], a, b, k0)
    try:
        popt, _ = curve_fit(f, ks, ps, p0=[0.5, 1.0],
                            bounds=([0.0, 0.05], [0.999, 20.0]), maxfev=40000)
        return float(popt[0]), float(popt[1])
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json", nargs="+")
    ap.add_argument("--k0", type=int, default=12)
    ap.add_argument("--out", default="uniqueness")
    ap.add_argument("--linear-x", action="store_true",
                    help="use a linear x-axis for the curves figure (default: log)")
    ap.add_argument("--xmin", type=float, default=5.0,
                    help="left x-limit for the log curves figure (default: 5)")
    args = ap.parse_args()

    genomes = [load(p) for p in args.json]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.dpi": 130})

    # Figure 1: uniqueness null vs repeat-aware fit
    fig, ax = plt.subplots(figsize=(9, 6))
    kmax = max(g["ks"].max() for g in genomes)
    kgrid = np.arange(1, kmax + 1)
    rows = []
    for i, g in enumerate(genomes):
        c = PALETTE[i % len(PALETTE)]
        null = np.exp(-g["N"] * np.power(g["Q"], kgrid))
        ax.plot(kgrid, null, ":", color=c, lw=1.3, alpha=0.7)
        ax.axvline(g["kc"], color=c, lw=0.7, alpha=0.3)
        step = max(1, len(g["ks"]) // 40)
        ax.scatter(g["ks"][::step], g["ps"][::step], color=c, s=14, zorder=5)
        ab = fit(g, args.k0)
        if ab:
            a, b = ab
            ax.plot(kgrid, model(kgrid, g["N"], g["Q"], a, b, args.k0),
                    "-", color=c, lw=2.0, label=f"{g['name']} (α={a:.2f}, β={b:.2f})")
            plateau = 1.0 - a
        else:
            ax.plot([], [], "-", color=c, lw=2.0, label=g["name"])
            a = b = plateau = float("nan")
        rows.append((g["name"], g["N"], g["Q"], g["kc"], a, b, plateau))
    # title="k-mer uniqueness: null (dotted) vs repeat-aware fit (solid)\npoints = exact p(k)=CDF(L); colored verticals = crossover $k_c$",
    ax.set(xlabel="k-mer length (k)", ylabel="Fraction of unique k-mers")
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_ylim(-0.02, 1.02)
    if not args.linear_x:
        ax.set_xscale("log")
        ax.set_xlim(left=max(1.0, args.xmin), right=kmax)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{args.out}_curves.pdf", bbox_inches="tight")

    # Figure 2: repeat spectrum S_rep(k) = Pr(L >= k)
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    for i, g in enumerate(genomes):
        if g["L"].size == 0:
            continue
        c = PALETTE[i % len(PALETTE)]
        tot = g["Lc"].sum() + g["overflow"]
        # survival: fraction of positions with match length >= L
        surv = 1.0 - (np.cumsum(g["Lc"]) - g["Lc"]) / tot
        ax2.plot(g["L"], surv, "-", color=c, lw=1.8, label=g["name"])
    ax2.set(title="Repeat spectrum: survival of maximal match length  S(k) = Pr(L >= k)",
            xlabel="match length L (bp)", ylabel="Pr(L >= k)")
    ax2.set_yscale("log"); ax2.set_xscale("log")
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(f"{args.out}_repeat_spectrum.pdf", bbox_inches="tight")

    # summary table
    print(f"{'genome':<12}{'N':>13}{'Q':>7}{'k_c':>7}{'alpha':>8}{'beta':>8}{'plateau':>9}")
    for name, N, Q, kc, a, b, pl in rows:
        af = f"{a:.3f}" if a == a else "  -  "
        bf = f"{b:.3f}" if b == b else "  -  "
        pf = f"{pl:.3f}" if pl == pl else "  -  "
        print(f"{name:<12}{N:>13.0f}{Q:>7.3f}{kc:>7.2f}{af:>8}{bf:>8}{pf:>9}")
    print(f"\nSaved: {args.out}_curves.pdf")
    print(f"Saved: {args.out}_repeat_spectrum.pdf")
    if not HAVE_SCIPY:
        print("[note] scipy missing: fitting skipped. pip install scipy --break-system-packages")


if __name__ == "__main__":
    main()
