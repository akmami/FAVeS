import json
from collections import defaultdict


def format_millions(value):
    """
    Convert integer to millions with one decimal (like 538.5)
    """
    return value # f"{value/1e6:.1f}"


def generate_latex_table(json_file, target_b):
    # Load JSON
    with open(json_file, "r") as f:
        full_data = json.load(f)

    # Filter by b
    data = [d for d in full_data if d["b"] == target_b]

    if not data:
        print(f"No entries found for b={target_b}")
        return

    # Group by (k, w)
    grouped = defaultdict(list)
    for d in data:
        grouped[(d["k"], d["w"])].append(d)

    # Sort groups by k, w
    groups = sorted(grouped.items())

    # Collect sorted n values
    n_values = sorted(set(d["n"] for d in data))

    # Header
    col_count = len(groups) * len(n_values)

    print(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}l" +
          "c" * col_count +
          r"@{\extracolsep{\fill}}}")
    print(r"\toprule")

    # First header row
    header1 = " & "
    for (k, w), group_data in groups:
        header1 += rf"\multicolumn{{{len(n_values)}}}{{c}}{{k={k}, w={w}}} & "
    print(header1.rstrip("& ") + r" \\")
    
    # Clines
    col_start = 2
    clines = ""
    for _ in groups:
        col_end = col_start + len(n_values) - 1
        clines += rf"\cline{{{col_start}-{col_end}}}"
        col_start = col_end + 1
    print(clines)

    # Second header row
    header2 = " & "
    for _ in groups:
        for n in n_values:
            header2 += f"n={n} & "
    print(header2.rstrip("& ") + r" \\")
    print(r"\midrule")

    # Row generator helper
    def make_row(label, key, formatter=lambda x: x):
        row = label + " & "
        for (_, _), group_data in groups:
            group_data = sorted(group_data, key=lambda x: x["n"])
            for n in n_values:
                entry = next((d for d in group_data if d["n"] == n), None)
                if entry:
                    row += formatter(entry[key]) + " & "
                else:
                    row += "- & "
        print(row.rstrip("& ") + r" \\")

    # Table rows
    make_row("Total \\# Seeds (M$^{*}$)", "total_seeds", format_millions)
    make_row("Distinct Seeds (M$^{*}$)", "distinct_seeds", format_millions)
    make_row("Unique Seeds (M$^{*}$)", "unique_seeds", format_millions)
    make_row("Avg Length (bp)", "avg_len", lambda x: f"{x:.2f}")
    make_row("StdDev Length (bp)", "std_dev_len", lambda x: f"{x:.2f}")
    make_row("Min/Max Length (bp)", "min_max_len")

    print(r"\bottomrule")
    print(r"\end{tabular*}")

    print()
    print("% ===== Unique Seed Ratio Table (All b merged) =====\n")

    # Build grouping for full dataset
    grouped_by_b = defaultdict(dict)
    for d in full_data:
        key = (d["k"], d["w"], d["n"])
        grouped_by_b[d["b"]][key] = d["unique_seeds_ratio"]

    all_b_values = sorted(grouped_by_b.keys())
    all_groups = sorted(set((d["k"], d["w"]) for d in full_data))
    all_n_values = sorted(set(d["n"] for d in full_data))

    col_count = len(all_groups) * len(all_n_values)

    print(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}l" +
        "c" * col_count +
        r"@{\extracolsep{\fill}}}")
    print(r"\toprule")

    # First header row
    header1 = " & "
    for (k, w) in all_groups:
        header1 += rf"\multicolumn{{{len(all_n_values)}}}{{c}}{{k={k}, w={w}}} & "
    print(header1.rstrip("& ") + r" \\")

    # Clines
    col_start = 2
    clines = ""
    for _ in all_groups:
        col_end = col_start + len(all_n_values) - 1
        clines += rf"\cline{{{col_start}-{col_end}}}"
        col_start = col_end + 1
    print(clines)

    # Second header row
    header2 = " & "
    for _ in all_groups:
        for n in all_n_values:
            header2 += f"n={n} & "
    print(header2.rstrip("& ") + r" \\")
    print(r"\midrule")

    # Rows = b values
    for b in all_b_values:
        row = f"b={b} & "
        for (k, w) in all_groups:
            for n in all_n_values:
                key = (k, w, n)
                ratio = grouped_by_b[b].get(key, "-")
                if ratio != "-":
                    row += f"{ratio:.3f} & "
                else:
                    row += "- & "
        print(row.rstrip("& ") + r" \\")

    print(r"\bottomrule")
    print(r"\end{tabular*}")

if __name__ == "__main__":
    generate_latex_table("../out/fa-blend/fuzzy-seeds-output.txt", target_b=32)