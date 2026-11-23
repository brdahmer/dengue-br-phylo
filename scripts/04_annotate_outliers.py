#!/usr/bin/env python3
"""
04_annotate_outliers.py -- explain TreeTime's clock outliers.

TreeTime reports outliers as bare accessions with a residual. That tells you
a tip is off the molecular clock; it does not tell you what those tips have
in common. This joins the outlier table back to the curated metadata and the
Nextclade lineage calls, so the question "is this a metadata error, a
submitter artefact, or a real lineage?" can be answered from evidence.

Reads  : <outdir>/outliers.tsv        (TreeTime)
         01_metadata/metadata_clean.csv
         02_qc/out/nextclade.tsv
Writes : <outdir>/outliers_annotated.tsv
Prints : a per-outlier table and a grouping summary

Usage
-----
  python3 scripts/04_annotate_outliers.py 05_timetree/clock_quick
  python3 scripts/04_annotate_outliers.py 05_timetree/full_stoch
"""

import csv, os, sys, collections

csv.field_size_limit(10 ** 9)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COLS = ["name", "given_date", "apparent_date", "residual", "direction",
        "collection_date", "date_precision", "uf", "municipality",
        "clade", "length", "ambig_frac", "nuc_completeness", "organization"]


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "05_timetree/clock_quick"
    outdir = outdir if os.path.isabs(outdir) else os.path.join(ROOT, outdir)
    ofile = os.path.join(outdir, "outliers.tsv")
    meta_f = os.path.join(ROOT, "01_metadata", "metadata_clean.csv")
    nc_f = os.path.join(ROOT, "02_qc", "out", "nextclade.tsv")
    for f in (ofile, meta_f, nc_f):
        if not os.path.exists(f):
            sys.exit("missing input: %s" % f)

    with open(meta_f, newline="") as fh:
        meta = {r["name"]: r for r in csv.DictReader(fh)}
    with open(nc_f, newline="") as fh:
        clade = {r["seqName"]: r["clade"] for r in csv.DictReader(fh, delimiter="\t")}

    with open(ofile, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        idx = {k: i for i, k in enumerate(h.strip() for h in header)}
        rows = [r for r in rd if r and r[0].strip()]

    recs = []
    for r in rows:
        name = r[0].strip()
        m = meta.get(name, {})
        try:
            resid = float(r[idx["residual"]])
        except (KeyError, ValueError):
            resid = float("nan")
        recs.append({
            "name": name,
            "given_date": r[idx["given_date"]] if "given_date" in idx else "",
            "apparent_date": r[idx["apparent_date"]] if "apparent_date" in idx else "",
            "residual": "%+.2f" % resid,
            # negative residual = less divergence than its date implies = below the line
            "direction": "below" if resid < 0 else "above",
            "collection_date": m.get("collection_date", "?"),
            "date_precision": m.get("date_precision", "?"),
            "uf": m.get("uf") or "-",
            "municipality": m.get("municipality") or "-",
            "clade": clade.get(name, "?"),
            "length": m.get("length", "?"),
            "ambig_frac": m.get("ambig_frac", "?"),
            "nuc_completeness": m.get("nuc_completeness", "?"),
            "organization": m.get("organization", "?"),
        })
    recs.sort(key=lambda x: (x["direction"], x["clade"], x["name"]))

    with open(os.path.join(outdir, "outliers_annotated.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(recs)

    print("=== %d clock outliers in %s ===\n" % (len(recs), os.path.basename(outdir)))
    hdr = "%-10s %-6s %-8s %-6s %-10s %-5s %-20s %-7s %s"
    print(hdr % ("name", "resid", "dir", "year", "clade", "uf", "municipality",
                 "ambig", "organization"))
    for r in recs:
        print(hdr % (r["name"], r["residual"], r["direction"],
                     (r["collection_date"] or "")[:4], r["clade"], r["uf"],
                     r["municipality"][:20], r["ambig_frac"],
                     r["organization"][:30]))

    print("\n=== what the outliers share ===")
    for field, label in (("clade", "lineage"), ("uf", "state"),
                         ("organization", "submitter"),
                         ("date_precision", "date precision"),
                         ("direction", "side of the line")):
        c = collections.Counter(r[field] for r in recs)
        print("  by %-16s %s" % (label,
              "  ".join("%s=%d" % (k[:30] or "-", v) for k, v in c.most_common())))

    print("\n=== groups (same lineage + state + side) ===")
    g = collections.defaultdict(list)
    for r in recs:
        g[(r["clade"], r["uf"], r["direction"])].append(r)
    for (cl, uf, d), members in sorted(g.items(), key=lambda x: -len(x[1])):
        munis = sorted({m["municipality"] for m in members})
        print("  %-10s %-4s %-6s n=%d  %s" % (cl, uf, d, len(members),
              ", ".join(munis)[:70]))
        if len(members) == 1:
            m = members[0]
            print("      %s: %s precision, %s, %s" %
                  (m["name"], m["date_precision"], m["nuc_completeness"],
                   m["organization"][:50]))

    print("\nwrote %s" % os.path.join(outdir, "outliers_annotated.tsv"))


if __name__ == "__main__":
    main()
