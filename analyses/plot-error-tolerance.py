import json
from collections import defaultdict

# ==========================================
# Load JSON
# ==========================================
with open("../out/err-sim/simulation-output.txt", "r") as f:
    full_data = json.load(f)

# ==========================================
# Configuration
# ==========================================
target_tuples = [(15, 10), (21, 11)]
b_values = [32, 37, 42]
n_values = [3, 5, 7]
error_rates = [0.001, 0.005, 0.01, 0.02]

metric = "precision"

# ==========================================
# Index data for fast lookup
# ==========================================
data_map = defaultdict(dict)

for entry in full_data:
    k, w, b, n = entry["k"], entry["w"], entry["b"], entry["n"]

    if (k, w) not in target_tuples:
        continue

    for stat in entry["stats"]:
        err = stat["error_rate"]
        value = stat[metric]
        data_map[(k, w, b, n, err)] = value

# ==========================================
# Generate tables (one per error rate)
# ==========================================
for err in error_rates:

    print("\n" + "="*80)
    print(f"% Error rate = {err}")
    print("="*80)

    print("\t\t\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lcccccc@{\\extracolsep{\\fill}}}")
    print("\t\t\t\\toprule")

    # Header row
    header = "& "
    for (k, w) in target_tuples:
        header += f"\\multicolumn{{3}}{{c}}{{k={k}, w={w}}} & "
    header = header.rstrip(" &")
    header += " \\\\"
    print(f"\t\t\t{header}")

    print("\t\t\t\\cline{2-4}\\cline{5-7}")
    print("\t\t\t& n=3 & n=5 & n=7 & n=3 & n=5 & n=7 \\\\")
    print("\t\t\t\\midrule")

    # Table body
    for b in b_values:
        row = f"b={b} "

        for (k, w) in target_tuples:
            for n in n_values:
                value = data_map.get((k, w, b, n, err), None)

                if value is None:
                    row += "& -- "
                else:
                    row += f"& {value:.3f} "

        row += "\\\\"
        print(f"\t\t\t{row}")

    print("\t\t\t\\bottomrule")
    print("\t\t\\end{tabular*}")
