import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


# ==========================================================
# configuration
# ==========================================================
#     if it satisfies ALL of the parameters listed in the dict.
#         {"k": 21, "w": 11, "b": 50, "n": 5}   -> records with exactly these
#     A parameter value may also be a list, meaning "any of these":
#         {"k": 21, "w": 11, "b": [25, 50], "n": 5}   -> b==25 OR b==50

OUT_DIR = Path("../out")
OUT_DIR.mkdir(parents=True, exist_ok=True)

PARAM_KEYS = ["k", "w", "s", "b", "n", "w_min", "w_max", "l", "d"]
STYLE_KEYS = {"color", "style", "marker", "label"}

METHOD_CONFIGS = [
    {
        "name": "Blend seeds",
        "path": "../out/fa-blend/blend-seeds-output.txt",
        "targets": [
            {"k": 21, "w": 11, "b": 50, "n": 5, "color": "#5A6F8E", "style": "-"},
        ],
    },
    {
        "name": "Minimizers",
        "path": "../out/fa-minimizer/minimizer-seeds-output.txt",
        "targets": [
            {"k": 21, "w": 11, "color": "#6BAF92", "style": "-"},
        ],
    },
    {
        "name": "Syncmer",
        "path": "../out/fa-syncmer/syncmer-seeds-output.txt",
        "targets": [
            {"k": 21, "s": 10, "color": "#D8C99B", "style": "-"},
        ],
    },
    {
        "name": "Strobemer",
        "path": "../out/fa-strobemer/strobemer-seeds-output.txt",
        "targets": [
            {"n": 2, "k": 22, "w_min": 6, "w_max": 16, "color": "#8E5A99", "style": "-"},
        ],
    },
    {
        "name": "LCP",
        "path": "../out/fa-lcp/lcp-seeds-output.txt",
        "targets": [
            {"l": 4, "d": 1, "color": "#4C956C", "style": "-"},
        ],
    },
]

PLOTS = [
    {
        "metric": "unique_seed_ratio",
        "ylabel": "Unique Seed Ratio",
        "title": "Uniqueness vs Relaxation",
        "ylim": (0, 1.0),
        "filename": "all_methods_uniqueness.png",
        "legend_loc": "lower right",
    },
    {
        "metric": "gap_ratio",
        "ylabel": "Gap Ratio",
        "title": "Gap Ratio vs Relaxation",
        "ylim": (0, 0.6),
        "filename": "all_methods_gap.png",
        "legend_loc": "upper right",
    },
]


# ==========================================================
# Styling defaults (used only when a target omits color/style/marker)
# ==========================================================
BASE_COLORS = [
    "#5A6F8E", "#6BAF92", "#D8C99B",
    "#C97C5D", "#8E5A99", "#4C956C",
    "#4477AA", "#CC6677", "#228833",
]

LINESTYLES = ["-", "--", ":", "-."]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
})

# ==========================================================
# Loading / normalization
# ==========================================================

def load_json(path):
    path = Path(path)
    with path.open("r") as f:
        return json.load(f)

def safe_value(record, key, default=None):
    return record[key] if key in record else default

def as_value_set(target):
    if isinstance(target, (list, tuple, set)):
        return set(target)
    return {target}

def split_target(target, method_index, target_index):
    param_filters = {}
    for key, val in target.items():
        if key in STYLE_KEYS:
            continue
        if key not in PARAM_KEYS:
            raise ValueError(
                f"unknown key '{key}'. Expected a parameter "
                f"({', '.join(PARAM_KEYS)}) or a style key "
                f"({', '.join(sorted(STYLE_KEYS))})."
            )
        param_filters[key] = val

    style = {
        "color": target.get("color", BASE_COLORS[method_index % len(BASE_COLORS)]),
        "style": target.get("style", LINESTYLES[target_index % len(LINESTYLES)]),
        "marker": target.get("marker", MARKERS[method_index % len(MARKERS)]),
        "label": target.get("label"),  # None -> auto-generated later
    }
    return param_filters, style

def record_matches(rec, param_filters):
    for param in PARAM_KEYS:
        if param in rec and param not in param_filters:
            raise ValueError(
                f"record contains parameter '{param}'={rec[param]} but the "
                f"target did not pin it. Add '{param}' to this target dict."
            )

    for param, value in param_filters.items():
        if param not in rec:
            return False
        if rec[param] not in as_value_set(value):
            return False

    return True

def make_row(method_name, rec, style):
    row = {
        "method": method_name,
        "details": rec.get("details", []),
        "color": style["color"],
        "style": style["style"],
        "marker": style["marker"],
        "label": style["label"],  # may be None until auto-filled
    }
    for param in PARAM_KEYS:
        row[param] = safe_value(rec, param)
    return row

def load_all_methods(method_configs):
    all_rows = []

    for m_idx, config in enumerate(method_configs):
        method_name = config["name"]
        path = config["path"]
        targets = config.get("targets", [])

        # be forgiving if I wrote a single dict instead of a list.
        if isinstance(targets, dict):
            targets = [targets]

        records = load_json(path)

        for t_idx, target in enumerate(targets):
            try:
                param_filters, style = split_target(target, m_idx, t_idx)
                for rec in records:
                    if record_matches(rec, param_filters):
                        all_rows.append(make_row(method_name, rec, style))
            except ValueError as err:
                raise ValueError(
                    f"[{method_name}] target #{t_idx + 1}: {err}"
                ) from err

    # fill in any labels the targets didn't set explicitly.
    varying = compute_varying_params(all_rows)
    for row in all_rows:
        if row["label"] is None:
            row["label"] = make_auto_label(row, varying)

    return all_rows


# ==========================================================
# Labels / titles
# ==========================================================

def compute_varying_params(rows):
    varying = {}
    for method in sorted(set(r["method"] for r in rows)):
        method_rows = [r for r in rows if r["method"] == method]
        vary = []
        for param in PARAM_KEYS:
            values = set(r[param] for r in method_rows if r[param] is not None)
            if len(values) > 1:
                vary.append(param)
        varying[method] = vary
    return varying


def make_auto_label(row, varying_params):
    parts = [row["method"]]
    for param in varying_params.get(row["method"], []):
        if row[param] is not None:
            parts.append(f"{param}={row[param]}")
    return ", ".join(parts)

# ==========================================================
# plotting
# ==========================================================

def extract_metric_from_details(row, metric):
    details = row["details"]

    values = []
    for d in details:
        if metric in d:
            values.append(d[metric])

    return np.array(values, dtype=float)


def plot_metric(
    rows,
    metric,
    ylabel,
    title,
    filename,
    ylim=None,
    legend_loc="best",
):
    plt.figure(figsize=(10, 6))

    max_r = 0
    legend_handles = {}

    for row in rows:
        y = extract_metric_from_details(row, metric)

        if len(y) == 0:
            continue

        r = np.arange(len(y))
        max_r = max(max_r, len(y) - 1)

        plt.plot(
            r,
            y,
            color=row["color"],
            linestyle=row["style"],
            marker=row["marker"],
            linewidth=2,
            markersize=5,
        )

        label = row["label"]
        if label not in legend_handles:
            legend_handles[label] = mlines.Line2D(
                [], [],
                color=row["color"],
                linestyle=row["style"],
                marker=row["marker"],
                linewidth=2,
                markersize=6,
                label=label,
            )

        print(label, "r=", r, metric, "=", y)

    plt.xlabel("Relaxation Radius (r)")
    plt.ylabel(ylabel)

    if ylim is not None:
        plt.ylim(*ylim)

    plt.xticks(np.arange(max_r + 1))

    plt.title(title)

    plt.legend(
        handles=list(legend_handles.values()),
        loc=legend_loc,
        frameon=True,
        ncol=2,
        fontsize=12,
    )

    plt.tight_layout()

    out_path = OUT_DIR / filename
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Saved: {out_path}")


# ==========================================================
# main
# ==========================================================

def main():
    rows = load_all_methods(METHOD_CONFIGS)

    if not rows:
        raise ValueError(
            "No records matched. Check the `targets` list of each method "
            "in METHOD_CONFIGS."
        )

    for plot_cfg in PLOTS:
        plot_metric(
            rows=rows,
            metric=plot_cfg["metric"],
            ylabel=plot_cfg["ylabel"],
            title=plot_cfg["title"],
            ylim=plot_cfg["ylim"],
            filename=plot_cfg["filename"],
            legend_loc=plot_cfg["legend_loc"],
        )


if __name__ == "__main__":
    main()
