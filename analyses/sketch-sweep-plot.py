#!/usr/bin/env python3
"""
Plot the sketching-invariance sweep (analysis 6).

Reads out/sweep/<genome>.<method>.json (each an array of sketch-seeds runs over
a k sweep) and plots per-seed uniqueness vs k, one line per sketching method,
one panel per genome. A companion figure shows seed density (thinning).

Usage:
  python3 sweep-plot.py                     # all genomes found -> small multiples
  python3 sweep-plot.py <genome> [short]    # single genome -> its own figure
Options:
  --sweep-dir DIR   (default ../out/sweep)
  --out-dir DIR     (default ../out)
"""
import sys, json, argparse, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.lines import Line2D

# method -> (label, color)   colors consistent with sketch-plot.py
METHOD_STYLE = {
    "kmer":       ("k-mer",       "#333333"),
    "minimizer":  ("Minimizers",  "#6BAF92"),
    "syncmer":    ("Syncmer",     "#D8C99B"),
    "blend":      ("Blend",       "#5A6F8E"),
    "fracmh":     ("FracMinHash", "#4477AA"),
    "strobemer":  ("Strobemer",   "#8E5A99"),
    "lcp":        ("LCP",         "#4C956C"),
}
METHOD_ORDER = list(METHOD_STYLE.keys())
# small genomes first (extend/edit freely); unknown names sort after, alphabetically
GENOME_ORDER = ["ecoliK12", "mtub", "sacCer3", "pf3D7", "dm6", "human"]


def _int(x):
    return int(str(x).replace(",", "")) if x is not None else 0


def load_sweep(sweep_dir):
    """-> {genome: {method: {'k':[...], 'uniq':[...], 'dens':[...], 'short':str}}}"""
    data = {}
    for p in sorted(Path(sweep_dir).glob("*.*.json")):
        stem = p.name[:-5]                 # strip .json
        if "." not in stem:
            continue
        genome, method = stem.rsplit(".", 1)
        try:
            runs = json.load(open(p))
        except Exception as e:
            print(f"  skip {p.name}: {e}")
            continue
        runs = [r for r in runs if isinstance(r, dict) and "k" in r]
        runs.sort(key=lambda r: r["k"])
        if not runs:
            continue
        ks   = [r["k"] for r in runs]
        uniq = [float(r.get("unique_seeds_ratio", 0.0)) for r in runs]
        dens = [(_int(r.get("total_seeds")) / r["valid_regions"]) if r.get("valid_regions") else float("nan")
                for r in runs]
        short = str(runs[0].get("shortname", "")).strip() or genome
        data.setdefault(genome, {})[method] = dict(k=ks, uniq=uniq, dens=dens, short=short)
    return data


def order_genomes(genomes):
    known = [g for g in GENOME_ORDER if g in genomes]
    rest  = sorted(g for g in genomes if g not in GENOME_ORDER)
    return known + rest


def order_methods(methods):
    return [m for m in METHOD_ORDER if m in methods] + \
           sorted(m for m in methods if m not in METHOD_ORDER)


def savefig(fig, prefix):
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")


def legend_for(methods):
    """proxy handles/labels for a set of methods, in canonical order."""
    handles, labels = [], []
    for m in order_methods(methods):
        label, color = METHOD_STYLE.get(m, (m, None))
        handles.append(Line2D([0], [0], color=color, marker="o", lw=2, ms=4))
        labels.append(label)
    return handles, labels


def draw_panel(ax, methods_data, metric, title):
    for m in order_methods(methods_data):
        d = methods_data[m]
        label, color = METHOD_STYLE.get(m, (m, None))
        ax.plot(d["k"], d[metric], "o-", color=color, lw=2, ms=4, label=label)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.set_xlabel("k")


def figure(data, genomes, metric, ylabel, title, percent, logy, out_prefix):
    n = len(genomes)
    ncol = min(3, n) if n > 1 else 1
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.8 * nrow),
                             squeeze=False, sharex=True)
    for i, g in enumerate(genomes):
        ax = axes[i // ncol][i % ncol]
        draw_panel(ax, data[g], metric, next(iter(data[g].values()))["short"])
        ax.set_ylabel(ylabel)
        if percent:
            ax.yaxis.set_major_formatter(PercentFormatter(1.0))
        if logy:
            ax.set_yscale("log")
    for j in range(n, nrow * ncol):            # hide empty panels
        axes[j // ncol][j % ncol].axis("off")
    # one shared legend, built from every method present across panels
    allm = set().union(*(set(data[g]) for g in genomes))
    handles, labels = legend_for(allm)
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    savefig(fig, out_prefix)
    plt.close(fig)
    print(f"wrote {out_prefix}.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome", nargs="?", help="single genome key (default: all)")
    ap.add_argument("shortname", nargs="?", help="unused; label comes from JSON")
    ap.add_argument("--sweep-dir", default="../out/sweep")
    ap.add_argument("--out-dir", default="../out")
    args = ap.parse_args()

    data = load_sweep(args.sweep_dir)
    if not data:
        print(f"no sweep JSON in {args.sweep_dir} (run: make sweep)"); sys.exit(1)

    out = Path(args.out_dir)
    if args.genome:                            # single genome
        if args.genome not in data:
            print(f"no sweep data for '{args.genome}' in {args.sweep_dir}"); sys.exit(1)
        (out / "sweep").mkdir(parents=True, exist_ok=True)
        short = next(iter(data[args.genome].values()))["short"]
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
        draw_panel(ax[0], data[args.genome], "uniq", f"Uniqueness")
        ax[0].set_ylabel("unique-seed ratio"); ax[0].yaxis.set_major_formatter(PercentFormatter(1.0))
        draw_panel(ax[1], data[args.genome], "dens", f"Seed density")
        ax[1].set_ylabel("seeds / bp"); ax[1].set_yscale("log")
        h, l = legend_for(set(data[args.genome]))
        fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=9, frameon=False)
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        savefig(fig, str(out / f"{args.genome}_sweep"))
        print(f"wrote {out/(args.genome+'_sweep')}.pdf")
    else:                                      # all genomes, small multiples
        genomes = order_genomes(list(data))
        figure(data, genomes, "uniq", "unique-seed ratio",
               "Sketching invariance: per-seed uniqueness vs k",
               percent=True, logy=False, out_prefix=str(out / "sweep_uniqueness"))
        figure(data, genomes, "dens", "seeds / bp",
               "Seed density (thinning) vs k",
               percent=False, logy=True, out_prefix=str(out / "sweep_density"))


if __name__ == "__main__":
    main()