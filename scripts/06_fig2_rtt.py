#!/usr/bin/env python3
"""Figure 2: Root-to-tip regression and clock outliers.

TreeTime's rtt.csv lists tips AND internal nodes. Internal-node dates are
inferred from the regression itself, so they sit on the line by construction
and must not be plotted or fitted. Year-only tips appear as 'YYYY.0-YYYY.999'
and are placed at the midpoint of that interval.
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- CONFIG ---
RTT_PATH = Path("05_timetree/clock_quick/rtt.csv")
META_PATH = Path("01_metadata/metadata_clean.csv")
NEXTCLADE_PATH = Path("02_qc/out/nextclade.tsv")
OUT_DIR = Path("06_figures")

ACC_COL = "name"
STATE_COL = "uf"
EXCLUDE_FROM_FIT = ["PP234964", "PV682577"]


def parse_date(value):
    """Return (decimal year, is_range). Ranges 'a-b' become their midpoint."""
    s = str(value).strip()
    if s.count(".") >= 2 and "-" in s:
        a, b = s.split("-")
        return (float(a) + float(b)) / 2, True
    return float(s), False


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rtt = pd.read_csv(RTT_PATH, comment="#", skipinitialspace=True)
    rtt.columns = [c.strip().lower() for c in rtt.columns]
    rtt = rtt.rename(columns={"root-to-tip distance": "rtt_dist"})

    # Tips only: internal nodes are named NODE_xxxxxxx
    rtt = rtt[~rtt["name"].str.startswith("NODE_")].copy()
    parsed = rtt["date"].apply(parse_date)
    rtt["num_date"] = [p[0] for p in parsed]
    rtt["is_range"] = [p[1] for p in parsed]
    rtt["rtt_dist"] = pd.to_numeric(rtt["rtt_dist"], errors="coerce")
    rtt = rtt.dropna(subset=["num_date", "rtt_dist"])

    meta = pd.read_csv(META_PATH)[[ACC_COL, STATE_COL]]
    nc = pd.read_csv(NEXTCLADE_PATH, sep="\t", low_memory=False)[["seqName", "clade"]]
    df = rtt.merge(meta, left_on="name", right_on=ACC_COL, how="left").merge(
        nc, left_on="name", right_on="seqName", how="left"
    )

    # Fit on all tips except the two records judged mis-dated / off-lineage
    fit_df = df[~df["name"].isin(EXCLUDE_FROM_FIT)]
    slope, intercept = np.polyfit(fit_df["num_date"], fit_df["rtt_dist"], 1)
    r2 = np.corrcoef(fit_df["num_date"], fit_df["rtt_dist"])[0, 1] ** 2
    root_date = -intercept / slope

    is_pp = df["name"] == "PP234964"
    is_pv = df["name"] == "PV682577"
    is_es_1va = (df[STATE_COL] == "ES") & (df["clade"] == "1V_A") & (df["num_date"] >= 2023)
    is_special = is_pp | is_pv | is_es_1va
    is_anchor = (df["num_date"] < 2000) & ~is_special
    is_main = ~(is_special | is_anchor)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.scatter(df.loc[is_main, "num_date"], df.loc[is_main, "rtt_dist"], s=14, alpha=0.45,
               color="#4c72b0", label=f"DENV-1 tips (n = {int(is_main.sum())})")
    ax.scatter(df.loc[is_anchor, "num_date"], df.loc[is_anchor, "rtt_dist"], s=40, color="#55a868",
               edgecolor="black", zorder=4, label=f"Pre-2000 tips (n = {int(is_anchor.sum())}, year-only dates)")
    ax.scatter(df.loc[is_es_1va, "num_date"], df.loc[is_es_1va, "rtt_dist"], s=40, color="#dd8452",
               edgecolor="black", zorder=4, label="ES 2023 1V_A outliers")
    ax.scatter(df.loc[is_pv, "num_date"], df.loc[is_pv, "rtt_dist"], s=60, color="#c44e52", marker="D",
               edgecolor="black", zorder=5, label="PV682577 (1V_H, likely mis-dated)")
    ax.scatter(df.loc[is_pp, "num_date"], df.loc[is_pp, "rtt_dist"], s=70, color="#8172b3", marker="^",
               edgecolor="black", zorder=5, label="PP234964 (1III_A.2)")

    x_line = np.linspace(df["num_date"].min() - 1, df["num_date"].max() + 1, 200)
    ax.plot(x_line, slope * x_line + intercept, color="black", linestyle="--", linewidth=1.5)

    for acc, label, offset in [
        ("PP234964", "PP234964\n(1III_A.2)", (-60, 25)),
        ("PV682577", "PV682577\n(1V_H)", (-60, -35)),
    ]:
        row = df[df["name"] == acc]
        if not row.empty:
            ax.annotate(label, xy=(row["num_date"].values[0], row["rtt_dist"].values[0]),
                        xytext=offset, textcoords="offset points", fontsize=8,
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    stats_txt = (
        f"Rate: {slope:.3e} subs/site/yr\n"
        f"$R^2$ = {r2:.3f}  (n = {len(fit_df)} tips)\n"
        f"Root date ≈ {root_date:.1f}"
    )
    ax.text(0.03, 0.95, stats_txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc"))

    ax.set_title("Root-to-Tip Genetic Divergence vs. Collection Date", fontweight="bold", loc="left")
    ax.set_xlabel("Collection Year (year-only dates at mid-year)")
    ax.set_ylabel("Root-to-Tip Distance (substitutions/site)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.82), fontsize=8, frameon=True)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"fig2_rtt.{ext}", dpi=300, bbox_inches="tight")
    print(f"tips plotted: {len(df)}; fit n={len(fit_df)}; rate={slope:.3e}; r2={r2:.3f}; root={root_date:.1f}")
    print(f"Saved {OUT_DIR}/fig2_rtt.[png|pdf]")


if __name__ == "__main__":
    main()
