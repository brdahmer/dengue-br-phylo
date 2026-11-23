#!/usr/bin/env python3
"""Cross-reference qc_bad.txt against curated and raw metadata."""

from pathlib import Path
import pandas as pd

# --- CONFIG ---
META_CLEAN = Path("01_metadata/metadata_clean.csv")
META_RAW = Path("01_metadata/all_denv_sequences_metadata.csv")  # set to None to skip
QC_BAD = Path("02_qc/qc_bad.txt")
OUT_CSV = Path("02_qc/qc_bad_verified.csv")


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

    bad_ids = [line.strip() for line in QC_BAD.read_text().splitlines() if line.strip()]
    bad_df = pd.DataFrame({"name": bad_ids}).merge(meta, on="name", how="left")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    bad_df.to_csv(OUT_CSV, index=False)

    print(f"QC BAD SEQUENCES (n = {len(bad_df)})")
    for col in ["year", "uf", "bioproject", "organization"]:
        if col in bad_df.columns:
            print(f"\nBy {col}:")
            print(bad_df[col].fillna("NA").value_counts().head(10).to_string())

    print(f"\nSaved -> {OUT_CSV}")


if __name__ == "__main__":
    main()