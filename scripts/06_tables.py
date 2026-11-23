#!/usr/bin/env python3
"""Generate Markdown tables for Lineage Summary (Table 2) and Clock Outliers (Table 3)."""

from pathlib import Path
import pandas as pd

# --- CONFIG ---
META_PATH = Path("01_metadata/metadata_clean.csv")
NEXTCLADE_PATH = Path("02_qc/out/nextclade.tsv")
DATES_TSV = Path("05_timetree/full_stoch/dates.tsv")  # TreeTime inferred dates & residuals
OUT_PATH = Path("06_figures/tables.md")

ACC_COL = "name"
YEAR_COL = "year"
DATE_COL = "collection_date"
STATE_COL = "uf"


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(META_PATH)
    nc = pd.read_csv(NEXTCLADE_PATH, sep="\t", low_memory=False)
    nc = nc[nc["qc.overallStatus"].isin(["good", "mediocre"])][["seqName", "clade"]]

    df = meta.merge(nc, left_on=ACC_COL, right_on="seqName", how="inner")
    if YEAR_COL not in df.columns:
        df[YEAR_COL] = pd.to_numeric(df[DATE_COL].astype(str).str[:4], errors="coerce")

    total = len(df)
    md = ["# Report Tables\n", "## Table 2: Lineage Summary\n"]
    md.append("| Lineage | Sequences (n) | Share (%) | Year Range | States (n) |")
    md.append("|---|---:|---:|---|---:|")

    for clade, grp in df.groupby("clade"):
        n = len(grp)
        pct = n / total * 100
        yr_min, yr_max = int(grp[YEAR_COL].min()), int(grp[YEAR_COL].max())
        yr_str = f"{yr_min}" if yr_min == yr_max else f"{yr_min}–{yr_max}"
        n_states = grp[STATE_COL].dropna().nunique()
        md.append(f"| `{clade}` | {n} | {pct:.2f}% | {yr_str} | {n_states} |")

    # Table 3: Clock Outliers (from TreeTime dates.tsv where date == '--')
    if DATES_TSV.exists():
        tt = pd.read_csv(DATES_TSV, sep="\t", comment="#")
        tt.columns = [c.strip() for c in tt.columns]
        # TreeTime marks clock-filtered tips with '--' in the numeric date column
        outliers = tt[tt.iloc[:, 1].astype(str).str.contains("--", na=False)]
        merged = outliers.merge(df, left_on=tt.columns[0], right_on=ACC_COL, how="left")

        md.extend(["\n## Table 3: Clock Outliers Flagged by TreeTime (`--clock-filter 4`)\n"])
        md.append("| Accession | Lineage | State | Stated Date |")
        md.append("|---|---|---|---|")
        for _, r in merged.iterrows():
            md.append(f"| `{r[tt.columns[0]]}` | `{r.get('clade', '')}` | {r.get(STATE_COL, '')} | {r.get(DATE_COL, '')} |")

    OUT_PATH.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()