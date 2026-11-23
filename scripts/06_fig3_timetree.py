#!/usr/bin/env python3
"""Figure 3: Time-scaled phylogeny colored by Nextclade lineage.

The x-axis is anchored on the node dates stored in timetree.nexus
(root_year = date - depth, constant across nodes). It must NOT be anchored on
'deepest tip = latest sampling date': clock outliers carry apparent dates far
beyond the sampling period (e.g. PP234964 at ~2046), which would shift the
whole tree. Tips with an apparent date after LATEST_TIP_YEAR are omitted.
"""

from pathlib import Path
import re
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

# --- CONFIG ---
TREE_PATH = Path("05_timetree/full_stoch/timetree.nexus")
NEXTCLADE_PATH = Path("02_qc/out/nextclade.tsv")
OUT_DIR = Path("06_figures")

LATEST_TIP_YEAR = 2026.0  # latest possible sampling date; later apparent dates are clock outliers
LINEAGE_COLORS = {
    "1V_E.1": "#1f77b4",
    "1V_D.1.1": "#ff7f0e",
    "1V_D.1.4": "#2ca02c",
    "1V_A": "#d62728",
    "1V_F": "#9467bd",
    "1V_F.1": "#8c564b",
    "1V_F.2": "#e377c2",
    "1V_H": "#17becf",
}


def read_node_dates(path: Path) -> dict:
    """Map node/tip name -> decimal date from [&date=...] comments."""
    txt = path.read_text()
    pat = r"([A-Za-z0-9_.]+):([0-9.eE\-]+)\[&[^\]]*?date=([0-9.]+)\]"
    return {m.group(1): float(m.group(3)) for m in re.finditer(pat, txt)}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tree = Phylo.read(TREE_PATH, "nexus")
    node_dates = read_node_dates(TREE_PATH)

    depths = tree.depths()
    offsets = [node_dates[t.name] - depths[t] for t in tree.get_terminals() if t.name in node_dates]
    root_year = float(np.median(offsets))

    # Drop tips whose apparent date is impossible (clock outliers)
    dropped = [t for t in tree.get_terminals()
               if node_dates.get(t.name, 0) > LATEST_TIP_YEAR]
    dropped_names = [t.name for t in dropped]
    for t in dropped:
        tree.prune(t)

    tree.ladderize()
    depths = tree.depths()

    nc = pd.read_csv(NEXTCLADE_PATH, sep="\t", low_memory=False)
    clade_map = dict(zip(nc["seqName"], nc["clade"]))

    y_pos = {}
    for idx, term in enumerate(tree.get_terminals()):
        y_pos[term] = idx
    for clade in tree.get_nonterminals(order="postorder"):
        y_pos[clade] = sum(y_pos[c] for c in clade.clades) / len(clade.clades)

    fig, ax = plt.subplots(figsize=(10, 12))

    for parent in tree.get_nonterminals():
        px = root_year + depths[parent]
        ys = [y_pos[c] for c in parent.clades]
        ax.plot([px, px], [min(ys), max(ys)], color="#444444", lw=0.4)
        for child in parent.clades:
            cx = root_year + depths[child]
            ax.plot([px, cx], [y_pos[child]] * 2, color="#444444", lw=0.4)

    tip_x, tip_y, tip_c = [], [], []
    for term in tree.get_terminals():
        tip_x.append(root_year + depths[term])
        tip_y.append(y_pos[term])
        tip_c.append(LINEAGE_COLORS.get(clade_map.get(term.name, "Other"), "#bbbbbb"))
    ax.scatter(tip_x, tip_y, c=tip_c, s=6, zorder=3, edgecolors="none")

    ax.set_title("Time-Scaled Phylogeny of Brazilian DENV-1 (TreeTime)", fontweight="bold", loc="left")
    ax.set_xlabel(f"Year  (omitted: {len(dropped_names)} clock-outlier tip(s) with apparent date > {LATEST_TIP_YEAR})")
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=col, markersize=6, label=lin)
        for lin, col in LINEAGE_COLORS.items()
    ] + [Line2D([0], [0], marker="o", color="w", markerfacecolor="#bbbbbb", markersize=6, label="Other")]
    ax.legend(handles=handles, title="Nextclade Lineage", loc="upper left", frameon=False, fontsize=8)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"fig3_timetree.{ext}", dpi=300, bbox_inches="tight")
    print(f"root_year={root_year:.2f}; tips drawn={len(tip_x)}; omitted={dropped_names}")
    print(f"tip x range: {min(tip_x):.1f}-{max(tip_x):.1f}")
    print(f"Saved {OUT_DIR}/fig3_timetree.[png|pdf]")


if __name__ == "__main__":
    main()
