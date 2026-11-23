#!/usr/bin/env python3
"""
03_tempest_dates.py -- write a TempEst-readable sampling-date file.

TempEst's File > Import Dates reads a TAB-delimited, header-less table of
taxon name and date. It cannot parse:
  * commas as separators
  * a header row
  * TreeTime interval syntax such as [2020.0:2020.999]
  * mixed ISO dates and decimal years

So this converts every date to a plain decimal year:
  full date  YYYY-MM-DD -> year + (day_of_year - 0.5) / days_in_year
  year-month YYYY-MM    -> the 15th of that month, same formula
  year only  YYYY       -> year + 0.5 (interval midpoint)

Reads  : 01_metadata/metadata_clean.csv (collection_date, date_precision)
Filters: to the names present in 02_qc/dates_qc_passed.csv
Writes : 05_timetree/tempest_dates.txt   (tab-delimited, no header)
         05_timetree/tempest_dates_precision.tsv  (audit trail)

Note: TempEst treats every value as exact. The 56 year-only records are
pinned to mid-year here, which TreeTime models as an interval instead --
so TempEst's regression is a visual diagnostic, not the reported estimate.
"""

import csv, os, sys
from calendar import isleap, monthrange
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "01_metadata", "metadata_clean.csv")
KEEP = os.path.join(ROOT, "02_qc", "dates_qc_passed.csv")
OUTD = os.path.join(ROOT, "05_timetree")


def decimal_year(y, m, d):
    doy = date(y, m, d).timetuple().tm_yday
    return y + (doy - 0.5) / (366.0 if isleap(y) else 365.0)


def convert(raw):
    """Return (decimal_year, precision_label) or (None, reason)."""
    raw = (raw or "").strip()
    parts = raw.split("-")
    try:
        if len(parts) == 3:
            return decimal_year(int(parts[0]), int(parts[1]), int(parts[2])), "full"
        if len(parts) == 2:
            y, m = int(parts[0]), int(parts[1])
            day = min(15, monthrange(y, m)[1])
            return decimal_year(y, m, day), "month"
        if len(parts) == 1 and len(raw) == 4:
            return int(raw) + 0.5, "year"
    except ValueError:
        pass
    return None, "unparseable"


def main():
    for f in (META, KEEP):
        if not os.path.exists(f):
            sys.exit("missing input: %s" % f)
    os.makedirs(OUTD, exist_ok=True)

    with open(KEEP, newline="") as fh:
        keep = {r["name"] for r in csv.DictReader(fh)}
    with open(META, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["name"] in keep]

    out, audit, bad = [], [], []
    for r in rows:
        val, lab = convert(r["collection_date"])
        if val is None:
            bad.append(r["name"])
            continue
        out.append((r["name"], val))
        audit.append((r["name"], r["collection_date"], lab, "%.6f" % val))

    out.sort()
    dates_path = os.path.join(OUTD, "tempest_dates.txt")
    with open(dates_path, "w") as fh:
        for name, val in out:
            fh.write("%s\t%.6f\n" % (name, val))

    with open(os.path.join(OUTD, "tempest_dates_precision.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "collection_date", "precision", "decimal_year"])
        w.writerows(sorted(audit))

    n = len(out)
    counts = {}
    for _, _, lab, _ in audit:
        counts[lab] = counts.get(lab, 0) + 1
    span = (min(v for _, v in out), max(v for _, v in out))
    print("wrote %s" % dates_path)
    print("  taxa: %d of %d expected" % (n, len(keep)))
    print("  precision: " + ", ".join("%s=%d" % kv for kv in sorted(counts.items())))
    print("  span: %.3f - %.3f" % span)
    if bad:
        print("  UNPARSEABLE (omitted, will show as 0 in TempEst): %s" % ", ".join(bad))
    if n != len(keep):
        print("  WARNING: %d taxa in the tree have no date" % (len(keep) - n))


if __name__ == "__main__":
    main()
