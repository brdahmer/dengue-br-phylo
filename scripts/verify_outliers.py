#!/usr/bin/env python3
"""Cross-reference TreeTime outliers.tsv against curated and raw metadata."""

from pathlib import Path
import pandas as pd

# --- CONFIG ---
META_CLEAN = Path("01_metadata/metadata_clean.csv")
META_RAW = Path("01_metadata/all_denv_sequences_metadata.csv")  # set to None to skip
OUTLIERS = Path("05_timetree/full_stoch/outliers.tsv")
OUT_CSV = Path("05_timetree/full_stoch/outliers_verified.csv")


def main() -> None:
    meta = pd.read_csv(META_CLEAN)

    if META_RAW and META_RAW.exists():
        raw = pd.read_csv(META_RAW, low_memory=False)
        raw["name"] = raw["Accession"].astype(str).str.split(".").str[0]
        extra = [
            c for c in ["name", "Isolate_Strain", "GenBank_Title", "BioSample", "Release_Date", "Publications"]
            if c in raw.columns
        ]
        meta = meta.merge(raw[extra].drop_duplicates(subset=["name"]), on="name", how="left")

    outliers = (
        pd.read_csv(OUTLIERS, sep="\t", index_col=0)
        .reset_index()
        .rename(columns={"index": "name"})
    )
    outliers["name"] = outliers["name"].astype(str).str.strip()
    out_df = outliers.merge(meta, on="name", how="left")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_CSV, index=False)

    print(f"CLOCK OUTLIERS (n = {len(out_df)})\n")
    preview_cols = [
        c for c in ["name", "collection_date", "given_date", "apparent_date", "residual", "uf", "municipality", "organization"]
        if c in out_df.columns
    ]
    print(out_df[preview_cols].to_string(index=False))
    print(f"\nSaved -> {OUT_CSV}")


if __name__ == "__main__":
    main()