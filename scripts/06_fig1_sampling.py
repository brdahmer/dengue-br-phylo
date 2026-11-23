#!/usr/bin/env python3
"""Figure 1: Temporal lineage turnover and geographic sampling by state."""

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

# --- CONFIG ---
META_PATH = Path("01_metadata/metadata_clean.csv")
NEXTCLADE_PATH = Path("02_qc/out/nextclade.tsv")
OUT_DIR = Path("06_figures")

ACC_COL = "name"
YEAR_COL = "year"
DATE_COL = "collection_date"
STATE_COL = "uf"

TOP_LINEAGES = ["1V_E.1", "1V_D.1.1", "1V_D.1.4", "1V_A", "1V_F", "1V_F.1", "1V_F.2", "1V_H"]
ALL_BR_STATES = [
    "AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", "MT", "MS",
    "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", "RS", "RO", "RR", "SC",
    "SP", "SE", "TO",
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(META_PATH)
    nc = pd.read_csv(NEXTCLADE_PATH, sep="\t", low_memory=False)
    nc = nc[nc["qc.overallStatus"].isin(["good", "mediocre"])][["seqName", "clade"]]

    df = meta.merge(nc, left_on=ACC_COL, right_on="seqName", how="inner")

    if YEAR_COL not in df.columns:
        df[YEAR_COL] = pd.to_numeric(df[DATE_COL].astype(str).str[:4], errors="coerce")
    df = df.dropna(subset=[YEAR_COL])
    df[YEAR_COL] = df[YEAR_COL].astype(int)

    df["lineage_grp"] = df["clade"].where(df["clade"].isin(TOP_LINEAGES), "Other")

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1.3, 1]})

    # Panel A: Sequences per year stacked by lineage
    piv = df.pivot_table(index=YEAR_COL, columns="lineage_grp", values=ACC_COL, aggfunc="count", fill_value=0)
    cols = [c for c in TOP_LINEAGES if c in piv.columns] + (["Other"] if "Other" in piv.columns else [])
    piv = piv[cols]
    piv.plot(kind="bar", stacked=True, ax=ax1, width=0.8, colormap="tab10")
    ax1.set_title("A. DENV-1 Genomes per Year by Lineage", fontweight="bold", loc="left")
    ax1.set_xlabel("Collection Year")
    ax1.set_ylabel("Number of QC-Passing Genomes")
    ax1.legend(title="Nextclade Lineage", fontsize=8, frameon=False)
    ax1.tick_params(axis="x", rotation=45, labelsize=8)

    # Panel B: Sequences by state (including 0-count states like MG)
    state_counts = df[STATE_COL].value_counts().reindex(ALL_BR_STATES, fill_value=0).sort_values(ascending=True)
    colors = ["#d62728" if s == "MG" else ("#1f77b4" if s == "RS" else "#7f7f7f") for s in state_counts.index]
    state_counts.plot(kind="barh", ax=ax2, color=colors)
    ax2.set_title("B. Geographic Distribution by State (UF)", fontweight="bold", loc="left")
    ax2.set_xlabel("Number of QC-Passing Genomes")
    ax2.set_ylabel("")
    ax2.tick_params(axis="y", labelsize=8)

    # Annotate MG = 0
    mg_idx = list(state_counts.index).index("MG")
    ax2.text(5, mg_idx, "0 (absent)", va="center", fontsize=8, color="#d62728", fontweight="bold")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"fig1_sampling.{ext}", dpi=300, bbox_inches="tight")
    print(f"Saved {OUT_DIR}/fig1_sampling.[png|pdf]")


if __name__ == "__main__":
    main()