from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import pandas as pd

TSV_PATH = "02_qc/out/nextclade.tsv"
META_PATH = "00_data/denv1_br_metadata.csv"  # optional, for inspecting outliers
OUT_PATH = Path("02_qc/summary.txt")

buf = StringIO()
with redirect_stdout(buf):
    df = pd.read_csv(TSV_PATH, sep="\t", low_memory=False)
    print(f"Loaded {len(df)} rows, {df.shape[1]} columns\n")

    # 1. Overall QC status
    print("=== QC OVERALL STATUS ===")
    print(df["qc.overallStatus"].value_counts(dropna=False))

    # 2. Alignment errors & reverse complements
    has_err = df["errors"].notna() & (df["errors"].astype(str).str.strip() != "")
    print(f"\nRows with alignment errors: {has_err.sum()}")
    print("\n=== REVERSE COMPLEMENT ===")
    print(df["isReverseComplement"].value_counts(dropna=False))

    # 3. Genotype & clade breakdown
    df["genotype"] = df["clade"].str.split("_").str[0]
    print("\n=== GENOTYPE COUNTS ===")
    print(df["genotype"].value_counts(dropna=False))

    print("\n=== CLADE BREAKDOWN (ALL) ===")
    print(df["clade"].value_counts(dropna=False))

    print("\n=== CLADE BREAKDOWN (GOOD + MEDIOCRE) ===")
    gm = df[df["qc.overallStatus"].isin(["good", "mediocre"])]
    print(gm["clade"].value_counts(dropna=False))

    # 4. Sub-QC failure reasons for 'bad' and 'mediocre'
    qc_cols = [
        "qc.missingData.status",
        "qc.mixedSites.status",
        "qc.privateMutations.status",
        "qc.snpClusters.status",
        "qc.frameShifts.status",
        "qc.stopCodons.status",
    ]

    for label in ["bad", "mediocre"]:
        sub = df[df["qc.overallStatus"] == label]
        print(f"\n=== SUB-QC STATUSES AMONG '{label}' (n={len(sub)}) ===")
        for col in qc_cols:
            print(f"  {col}: {sub[col].value_counts(dropna=False).to_dict()}")

    # 5. Coverage, missing data, frameshifts, stop codons by QC status
    print("\n=== COVERAGE & MISSING DATA STATS ===")
    for status in ["good", "mediocre", "bad"]:
        sub = df[df["qc.overallStatus"] == status]
        if len(sub) > 0:
            print(f"\nStatus={status} (n={len(sub)}):")
            print("  coverage:", sub["coverage"].describe()[["min", "50%", "mean", "max"]].to_dict())
            print("  totalMissing:", sub["totalMissing"].describe()[["min", "50%", "mean", "max"]].to_dict())
            print("  totalFrameShifts > 0:", (sub["totalFrameShifts"] > 0).sum())
            print("  totalStopCodons > 0:", (sub["qc.stopCodons.totalStopCodons"] > 0).sum())

    # 6. Inspect non-missingData 'bad' sequences
    print("\n=== BAD SEQUENCES WHERE missingData != 'bad' ===")
    bad_other = df[(df["qc.overallStatus"] == "bad") & (df["qc.missingData.status"] != "bad")]
    cols = [
        "seqName", "clade", "qc.overallScore", "qc.missingData.status",
        "qc.privateMutations.status", "qc.frameShifts.status",
        "totalMissing", "totalFrameShifts", "privateNucMutations.totalPrivateSubstitutions",
    ]
    print(bad_other[cols].to_string(index=False))

    # 7. Inspect the Genotype III sequence (1III_A.2)
    print("\n=== NON-GENOTYPE V SEQUENCES ===")
    non_5 = df[df["genotype"] != "1V"]
    print(non_5[["seqName", "clade", "qc.overallStatus", "coverage", "totalMissing", "totalSubstitutions"]].to_string(index=False))

    # Optional: cross-check PP234964 in the raw NCBI metadata
    try:
        meta = pd.read_csv(META_PATH)
        print("\n=== PP234964 IN RAW METADATA ===")
        print(meta[meta["Accession"].str.startswith("PP234964")].to_string(index=False))
    except FileNotFoundError:
        pass

output = buf.getvalue()
print(output, end="")
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
OUT_PATH.write_text(output, encoding="utf-8")