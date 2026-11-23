#!/usr/bin/env python3
"""
01_curate.py -- curate NCBI Virus dengue downloads into analysis-ready inputs.

Reads  : 00_data/denv<N>_br.fasta, 00_data/denv<N>_br_metadata.csv
Writes : 01_metadata/curated.fasta
         01_metadata/metadata_clean.csv
         01_metadata/dates.csv            (name,date for TreeTime)
         01_metadata/attrition.md         (per-filter record counts)
         01_metadata/unresolved_geo.csv   (records needing manual review)

Design notes
------------
* Sequence IDs are the bare accession (version suffix stripped). NCBI Virus
  exports the Accession column with or without the ".N" version depending on
  the download; stripping it on both sides is the only join that is stable.
* Geo_Location arrives in six mutually incompatible formats. Position is NOT
  reliable -- "Goiania, GO" is City,UF but "Sao Paulo, Taubate" is State,City.
  So each token is resolved against a UF lookup and whichever side matches a
  state wins, rather than assuming an order.
* Records with unresolvable geography or imprecise dates are FLAGGED and kept,
  not silently dropped. Dropping them is a judgement call for the analyst;
  the counts are in attrition.md so the call is visible.
* Finer sequence QC is delegated to Nextclade (stage 2). This script removes
  only catastrophically ambiguous sequences (default >30% non-ACGT), because
  those degrade the alignment that Nextclade's QC is measured against.

Stdlib only, deliberately: `python3 scripts/01_curate.py` must work in a clean
checkout with no environment.

Usage
-----
  python3 scripts/01_curate.py                 # DENV-1, default thresholds
  python3 scripts/01_curate.py --serotype 2
  python3 scripts/01_curate.py --self-test     # run unit tests, touch no data
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import unicodedata
from collections import Counter, OrderedDict

# --------------------------------------------------------------------------
# Brazilian states: UF code -> canonical name
# --------------------------------------------------------------------------

UF_TO_NAME = OrderedDict([
    ("AC", "Acre"), ("AL", "Alagoas"), ("AP", "Amapa"), ("AM", "Amazonas"),
    ("BA", "Bahia"), ("CE", "Ceara"), ("DF", "Distrito Federal"),
    ("ES", "Espirito Santo"), ("GO", "Goias"), ("MA", "Maranhao"),
    ("MT", "Mato Grosso"), ("MS", "Mato Grosso do Sul"), ("MG", "Minas Gerais"),
    ("PA", "Para"), ("PB", "Paraiba"), ("PR", "Parana"), ("PE", "Pernambuco"),
    ("PI", "Piaui"), ("RJ", "Rio de Janeiro"), ("RN", "Rio Grande do Norte"),
    ("RS", "Rio Grande do Sul"), ("RO", "Rondonia"), ("RR", "Roraima"),
    ("SC", "Santa Catarina"), ("SP", "Sao Paulo"), ("SE", "Sergipe"),
    ("TO", "Tocantins"),
])

# Spelling variants seen in the wild, mapped to UF. Keys are fold()ed.
NAME_ALIASES = {
    "amapa": "AP", "ceara": "CE", "espirito santo": "ES", "goias": "GO",
    "goiaz": "GO", "maranhao": "MA", "para": "PA", "paraiba": "PB",
    "parana": "PR", "piaui": "PI", "rondonia": "RO", "sao paulo": "SP",
    "brasilia": "DF", "distrito federal": "DF",
}

NOISE_WORDS = {"state", "estado", "brazil", "brasil", "province", "uf"}


def fold(s: str) -> str:
    """Lowercase, strip accents, collapse whitespace. For matching only."""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"\s+", " ", s).strip().lower()


NAME_TO_UF = {fold(v): k for k, v in UF_TO_NAME.items()}
NAME_TO_UF.update(NAME_ALIASES)

# Brazilian name connectors stay lowercase when title-casing
CONNECTORS = {"de", "da", "do", "das", "dos", "e", "d"}


def split_camel(tok: str) -> str:
    """FeiraDeSantana -> Feira De Santana. Leaves all-caps UF codes alone."""
    return re.sub(r"(?<=[a-z])(?=[A-Z])", " ", tok)


def titlecase_br(s: str) -> str:
    """Title-case a place name, keeping Portuguese connectors lowercase."""
    words = s.split()
    out = []
    for i, w in enumerate(words):
        lw = fold(w)
        if i > 0 and lw in CONNECTORS:
            out.append(lw)
        else:
            out.append(w[:1].upper() + w[1:].lower() if w else w)
    return " ".join(out)


def clean_token(tok: str) -> str:
    """Normalise one Geo_Location token to a comparable place name."""
    tok = tok.replace(")", " ").replace("(", " ")
    tok = split_camel(tok)
    tok = re.sub(r"[_/]+", " ", tok)
    tok = re.sub(r"\s+", " ", tok).strip(" .,;:-")
    # drop noise words, but never reduce the token to nothing
    kept = [w for w in tok.split() if fold(w) not in NOISE_WORDS]
    if kept:
        tok = " ".join(kept)
    return tok.strip()


def match_state(tok: str):
    """Return UF code if tok names a state (by code or name), else None."""
    t = clean_token(tok)
    if not t:
        return None
    if re.fullmatch(r"[A-Za-z]{2}", t) and t.upper() in UF_TO_NAME:
        return t.upper()
    return NAME_TO_UF.get(fold(t))


def split_geo_tokens(body: str):
    """Split the post-'Brazil:' part of Geo_Location into candidate tokens."""
    # parenthetical form: "Tocantins state(Porto Nacional)"
    body = re.sub(r"\(", ";", body).replace(")", " ")
    parts = re.split(r"[;,]", body)
    return [p for p in (p.strip() for p in parts) if p]


def parse_geo(raw: str):
    """
    Resolve a Geo_Location string.

    Returns (uf, state_name, municipality, status) where status is one of
    'ok', 'state_only', 'country_only', 'unresolved'.
    """
    raw = (raw or "").strip()
    body = raw.split(":", 1)[1] if ":" in raw else ""
    tokens = split_geo_tokens(body)

    if not tokens:
        return (None, None, None, "country_only")

    hits = [(i, match_state(t)) for i, t in enumerate(tokens)]
    state_idx = [i for i, uf in hits if uf]

    if not state_idx:
        return (None, None, titlecase_br(clean_token(tokens[0])), "unresolved")

    # Both tokens name a state: prefer an explicit UF code as the state,
    # so "Sao Paulo, SP" -> state SP, municipality Sao Paulo.
    if len(state_idx) > 1:
        coded = [i for i in state_idx
                 if re.fullmatch(r"[A-Za-z]{2}", clean_token(tokens[i]))]
        chosen = coded[0] if coded else state_idx[0]
    else:
        chosen = state_idx[0]

    uf = dict(hits)[chosen]
    others = [t for i, t in enumerate(tokens) if i != chosen]
    muni = titlecase_br(clean_token(others[0])) if others else None

    # "Rio de Janeiro, Rio de Janeiro" -- capital shares the state name
    if muni is None and len(tokens) == 1:
        return (uf, UF_TO_NAME[uf], None, "state_only")
    return (uf, UF_TO_NAME[uf], muni, "ok")


# --------------------------------------------------------------------------
# Dates
# --------------------------------------------------------------------------

def parse_date(raw: str):
    """
    Returns (precision, treetime_value) where precision is
    'full' | 'month' | 'year' | 'missing'.

    TreeTime accepts ISO dates and, for imprecise records, a [lo:hi] range.
    Year-month is given the 15th; year-only becomes a one-year range so the
    uncertainty is modelled rather than invented.
    """
    raw = (raw or "").strip()
    if not raw:
        return ("missing", None)
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", raw):
        return ("full", raw)
    if re.fullmatch(r"\d{4}-\d{2}", raw):
        return ("month", raw + "-15")
    if re.fullmatch(r"\d{4}", raw):
        y = int(raw)
        return ("year", "[%d.0:%d.999]" % (y, y))
    return ("missing", None)


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def bare(acc: str) -> str:
    """Strip an NCBI version suffix: 'OZ570380.1' -> 'OZ570380'."""
    return re.sub(r"\.\d+$", "", (acc or "").strip())


def read_fasta(path):
    seqs = OrderedDict()
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = bare(line[1:].split()[0])
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def ambig_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    return sum(1 for c in seq if c not in "ACGT") / len(seq)


# --------------------------------------------------------------------------
# Attrition log
# --------------------------------------------------------------------------

class Attrition:
    """Per-filter record counts.

    Removals are reported as a share of that step's own input; retention as a
    share of the starting total, so both the local cost of a filter and the
    cumulative survival are readable from one row.
    """

    def __init__(self):
        self.steps = []
        self.baseline = None
        self.baseline_label = None

    def add(self, label, kept, removed=0, note=""):
        self.steps.append((label, kept, removed, note))

    def set_baseline(self, n, label):
        """Anchor retention percentages here rather than at the raw download.

        The download is all four serotypes, so measuring DENV-1 retention
        against it would count the other serotypes as attrition.
        """
        self.baseline = n
        self.baseline_label = label

    def render(self, title):
        start = (self.steps[0][1] + self.steps[0][2]) if self.steps else 0
        base = self.baseline or start
        anchor = self.baseline_label or ("%d records downloaded" % start)
        out = ["# %s" % title, "",
               "Record counts after each curation step. Generated by",
               "`scripts/01_curate.py`; do not edit by hand.", "",
               "Removed is given as a share of that step's own input.",
               "Retained is a share of the %s." % anchor, "",
               "| Step | Removed | Retained |",
               "|---|---:|---:|"]
        for label, kept, removed, note in self.steps:
            lab = label + ((" — " + note) if note else "")
            step_in = kept + removed
            rem = ("%d/%d (%.1f%%)" % (removed, step_in, 100.0 * removed / step_in)
                   if removed and step_in else "—")
            ret = ("%d/%d (%.1f%%)" % (kept, base, 100.0 * kept / base)
                   if base and kept <= base else "%d" % kept)
            out.append("| %s | %s | %s |" % (lab, rem, ret))
        return "\n".join(out) + "\n"


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def curate(root, serotype, max_ambig, drop_unresolved_geo, drop_year_only,
           stem=None, outdir_name="01_metadata"):
    stem = stem or ("denv%d_br" % serotype)
    fa_in = os.path.join(root, "00_data", "%s.fasta" % stem)
    csv_in = os.path.join(root, "00_data", "%s_metadata.csv" % stem)
    outdir = os.path.join(root, outdir_name)
    os.makedirs(outdir, exist_ok=True)

    for p in (fa_in, csv_in):
        if not os.path.exists(p):
            sys.exit("missing input: %s" % p)

    seqs = read_fasta(fa_in)
    with open(csv_in, newline="") as fh:
        rows = list(csv.DictReader(fh))

    att = Attrition()
    att.add("Downloaded from NCBI Virus", len(rows), 0,
            "%d FASTA records" % len(seqs))

    # --- join -------------------------------------------------------------
    meta = {}
    dup_acc = 0
    for r in rows:
        k = bare(r.get("Accession", ""))
        if not k:
            continue
        if k in meta:
            dup_acc += 1
            continue
        meta[k] = r
    joined = [k for k in seqs if k in meta]
    att.add("Joined FASTA to metadata on accession",
            len(joined), len(seqs) - len(joined),
            "version suffix stripped both sides")

    # --- organism / serotype sanity --------------------------------------
    want = re.compile(r"type\s*%s\b|type\s*%s\b" %
                      (serotype, "IV" if serotype == 4 else
                       "III" if serotype == 3 else
                       "II" if serotype == 2 else "I"), re.I)
    keep = [k for k in joined if want.search(meta[k].get("Organism_Name", ""))]
    att.add("Organism name matches DENV-%d" % serotype,
            len(keep), len(joined) - len(keep),
            "scoping step, not attrition: the download covers all four serotypes")
    att.set_baseline(len(keep), "%d DENV-%d records retrieved" % (len(keep), serotype))

    # --- excluded record types -------------------------------------------
    # Lab-adapted material violates clock assumptions. Older GenBank records
    # often name the culture system in the title, isolate or source field.
    bad_re = re.compile(
        r"patent|synthetic|construct|vaccine|infectious clone|chimeric|"
        r"passag|recombinant|C6/36|C6-36|LLC-MK2|suckling mouse|mouse brain|"
        r"vero cell|cell cultur|cell-cultur", re.I)
    def is_bad(k):
        r = meta[k]
        blob = " ".join([r.get("GenBank_Title", ""), r.get("Isolate_Strain", ""),
                         r.get("Tissue_Specimen_Source", "")])
        return bool(bad_re.search(blob))
    before = len(keep)
    keep = [k for k in keep if not is_bad(k)]
    att.add("Removed patent / synthetic / clone / passaged / cell-cultured records",
            len(keep), before - len(keep))

    # --- country -----------------------------------------------------------
    before = len(keep)
    keep = [k for k in keep if fold(meta[k].get("Country", "")) == "brazil"]
    att.add("Country is Brazil", len(keep), before - len(keep))

    # --- catastrophic ambiguity -------------------------------------------
    amb = {k: ambig_fraction(seqs[k].upper()) for k in keep}
    before = len(keep)
    keep = [k for k in keep if amb[k] <= max_ambig]
    att.add("Ambiguity filter", len(keep), before - len(keep),
            ">%.0f%% non-ACGT removed; finer QC left to Nextclade"
            % (max_ambig * 100))

    # --- dates -------------------------------------------------------------
    dates = {k: parse_date(meta[k].get("Collection_Date", "")) for k in keep}
    before = len(keep)
    keep = [k for k in keep if dates[k][0] != "missing"]
    att.add("Collection date parseable", len(keep), before - len(keep))

    prec = Counter(dates[k][0] for k in keep)
    if drop_year_only:
        before = len(keep)
        keep = [k for k in keep if dates[k][0] != "year"]
        att.add("Dropped year-only collection dates",
                len(keep), before - len(keep))
    else:
        att.add("Year-only collection dates retained as TreeTime ranges",
                len(keep), 0, "%d records" % prec.get("year", 0))

    # --- geography ---------------------------------------------------------
    geo = {k: parse_geo(meta[k].get("Geo_Location", "")) for k in keep}
    unresolved = [k for k in keep if geo[k][3] in ("unresolved", "country_only")]
    if drop_unresolved_geo:
        before = len(keep)
        keep = [k for k in keep if k not in set(unresolved)]
        att.add("Dropped records without a resolvable state",
                len(keep), before - len(keep))
    else:
        att.add("Records without a resolvable state retained and flagged",
                len(keep), 0, "%d records" % len(unresolved))

    # --- duplicates --------------------------------------------------------
    # identical sequence strings: keep the first accession seen
    seen_seq, before = {}, len(keep)
    dedup = []
    for k in keep:
        sig = seqs[k].upper().replace("-", "")
        if sig in seen_seq:
            continue
        seen_seq[sig] = k
        dedup.append(k)
    keep = dedup
    att.add("Collapsed identical sequences", len(keep), before - len(keep))

    # repeated isolate names: keep the least ambiguous, then the longest
    by_iso = {}
    for k in keep:
        iso = (meta[k].get("Isolate_Strain") or "").strip()
        by_iso.setdefault(iso or "__blank__%s" % k, []).append(k)
    before = len(keep)
    chosen = []
    for iso, ks in by_iso.items():
        if len(ks) == 1:
            chosen.append(ks[0])
        else:
            chosen.append(sorted(ks, key=lambda k: (amb[k], -len(seqs[k])))[0])
    keep = [k for k in keep if k in set(chosen)]
    att.add("Collapsed repeated isolate names", len(keep), before - len(keep))

    # --- outputs -----------------------------------------------------------
    keep.sort()
    write_fasta(os.path.join(outdir, "curated.fasta"),
                OrderedDict((k, seqs[k].upper()) for k in keep))

    cols = ["name", "accession_versioned", "collection_date", "date_precision",
            "treetime_date", "year", "uf", "state", "municipality",
            "geo_status", "geo_raw", "host", "length", "ambig_frac",
            "nuc_completeness", "bioproject", "organization", "submitters"]
    with open(os.path.join(outdir, "metadata_clean.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for k in keep:
            r = meta[k]
            uf, state, muni, status = geo[k]
            p, tt = dates[k]
            w.writerow([
                k, (r.get("Accession") or "").strip(),
                (r.get("Collection_Date") or "").strip(), p, tt,
                (r.get("Collection_Date") or "")[:4],
                uf or "", state or "", muni or "", status,
                (r.get("Geo_Location") or "").strip(),
                (r.get("Host") or "").strip(), len(seqs[k]),
                "%.4f" % amb[k], (r.get("Nuc_Completeness") or "").strip(),
                (r.get("BioProject") or "").strip(),
                (r.get("Organization") or "").strip(),
                (r.get("Submitters") or "").strip(),
            ])

    with open(os.path.join(outdir, "dates.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "date"])
        for k in keep:
            w.writerow([k, dates[k][1]])

    still_unresolved = [k for k in keep
                        if geo[k][3] in ("unresolved", "country_only")]
    with open(os.path.join(outdir, "unresolved_geo.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "geo_raw", "geo_status", "guessed_municipality"])
        for k in still_unresolved:
            w.writerow([k, (meta[k].get("Geo_Location") or "").strip(),
                        geo[k][3], geo[k][2] or ""])

    att.add("**Final curated set**", len(keep), 0, "written to `curated.fasta`")
    with open(os.path.join(outdir, "attrition.md"), "w") as fh:
        fh.write(att.render("Curation attrition — DENV-%d, Brazil" % serotype))

    # --- console summary ---------------------------------------------------
    print(att.render("Curation attrition — DENV-%d, Brazil" % serotype))
    print("Date precision: " + ", ".join("%s=%d" % kv for kv in prec.items()))
    st = Counter(geo[k][0] or "(none)" for k in keep)
    print("States represented: %d" % len([s for s in st if s != "(none)"]))
    print("  " + ", ".join("%s=%d" % kv for kv in st.most_common(12)))
    print("Municipalities resolved: %d"
          % len({geo[k][2] for k in keep if geo[k][2]}))
    print("Geo status: " + ", ".join(
        "%s=%d" % kv for kv in Counter(geo[k][3] for k in keep).items()))
    print("Wrote 4 files + attrition.md to 01_metadata/")


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------

GEO_CASES = [
    # raw                                        uf    municipality      status
    ("Brazil: Goiania, GO",                      "GO", "Goiania",        "ok"),
    ("Brazil: Anapolis, GO",                     "GO", "Anapolis",       "ok"),
    ("Brazil: Cuiaba, MT",                       "MT", "Cuiaba",         "ok"),
    ("Brazil: Chapadao do Sul, Mato Grosso do Sul", "MS", "Chapadao do Sul", "ok"),
    ("Brazil: Salvador, BA",                     "BA", "Salvador",       "ok"),
    # order reversed: State, City
    ("Brazil: Sao Paulo, Taubate",               "SP", "Taubate",        "ok"),
    ("Brazil: Sao Paulo, Ribeirao Preto",        "SP", "Ribeirao Preto", "ok"),
    # semicolon form
    ("Brazil:Rio Grande do Sul; Canoas",         "RS", "Canoas",         "ok"),
    # CamelCase, no spaces
    ("Brazil: FeiraDeSantana, BA",               "BA", "Feira de Santana", "ok"),
    ("Brazil: AparecidaDeGoiania, GO",           "GO", "Aparecida de Goiania", "ok"),
    ("Brazil: AltoParaisoDeGoias, GO",           "GO", "Alto Paraiso de Goias", "ok"),
    ("Brazil: RibasDoRioPardo, MS",              "MS", "Ribas do Rio Pardo", "ok"),
    # parenthetical
    ("Brazil: Tocantins state(Porto Nacional)",  "TO", "Porto Nacional", "ok"),
    # UF code only
    ("Brazil: PB",                               "PB", None,       "state_only"),
    ("Brazil: SP",                               "SP", None,       "state_only"),
    # state name only
    ("Brazil: Mato Grosso",                      "MT", None,       "state_only"),
    ("Brazil: Goias",                            "GO", None,       "state_only"),
    ("Brazil: Alagoas",                          "AL", None,       "state_only"),
    ("Brazil: Distrito Federal",                 "DF", None,       "state_only"),
    # capital shares the state name
    ("Brazil: Rio de Janeiro, Rio de Janeiro",   "RJ", "Rio de Janeiro", "ok"),
    # both tokens resolve; UF code should win as the state
    ("Brazil: Sao Paulo, SP",                    "SP", "Sao Paulo",      "ok"),
    # accented input (none in the current file, but cheap insurance)
    ("Brazil: São Paulo, Taubaté",     "SP", "Taubate",        "ok"),
    ("Brazil: Goiânia, GO",                 "GO", "Goiania",        "ok"),
    # no subdivision at all
    ("Brazil",                                   None, None,  "country_only"),
    ("Brazil:",                                  None, None,  "country_only"),
    # genuinely unresolvable municipality
    ("Brazil: Nowhereville",                     None, "Nowhereville", "unresolved"),
    # must not substring-match a longer state name
    ("Brazil: Mato Grosso do Sul",               "MS", None,       "state_only"),
]

DATE_CASES = [
    ("2024-05-17", "full",  "2024-05-17"),
    ("2019-03",    "month", "2019-03-15"),
    ("2021",       "year",  "[2021.0:2021.999]"),
    ("",           "missing", None),
    ("   ",        "missing", None),
    ("not a date", "missing", None),
]

EXCLUDE_CASES = [
    ("DENV-1 isolate BR/90, passage 3 in C6/36 cells", True),
    ("dengue virus 1 strain X, suckling mouse brain passage", True),
    ("infectious clone of DENV-1", True),
    ("DENV-1 vaccine candidate construct", True),
    ("dengue virus type 1 isolate sample-77337 genome assembly", False),
    ("dengue virus type 1 strain BR/RJ/2024, complete genome", False),
]

ACC_CASES = [("OZ570380.1", "OZ570380"), ("MK517773", "MK517773"),
             ("PP269455.12", "PP269455"), ("", "")]


def self_test():
    fails = []
    for raw, uf, muni, status in GEO_CASES:
        got = parse_geo(raw)
        exp = (uf, UF_TO_NAME[uf] if uf else None, muni, status)
        # compare municipality case-insensitively on accent-folded form
        ok = (got[0] == exp[0] and got[3] == exp[3] and
              fold(got[2] or "") == fold(exp[2] or ""))
        if not ok:
            fails.append("GEO  %-46r\n       got %r\n       want %r"
                         % (raw, got, exp))
    for raw, prec, val in DATE_CASES:
        got = parse_date(raw)
        if got != (prec, val):
            fails.append("DATE %-46r got %r want %r" % (raw, got, (prec, val)))
    for raw, exp in ACC_CASES:
        if bare(raw) != exp:
            fails.append("ACC  %-46r got %r want %r" % (raw, bare(raw), exp))
    bad_re = re.compile(
        r"patent|synthetic|construct|vaccine|infectious clone|chimeric|"
        r"passag|recombinant|C6/36|C6-36|LLC-MK2|suckling mouse|mouse brain|"
        r"vero cell|cell cultur|cell-cultur", re.I)
    for text, want in EXCLUDE_CASES:
        got = bool(bad_re.search(text))
        if got != want:
            fails.append("EXCL %-46r got %r want %r" % (text, got, want))
    # ambiguity maths
    if abs(ambig_fraction("ACGTNNNN") - 0.5) > 1e-9:
        fails.append("AMBIG ambig_fraction('ACGTNNNN') != 0.5")
    if ambig_fraction("ACGTACGT") != 0.0:
        fails.append("AMBIG clean sequence should be 0.0")

    total = (len(GEO_CASES) + len(DATE_CASES) + len(ACC_CASES)
             + len(EXCLUDE_CASES) + 2)
    if fails:
        print("SELF-TEST FAILED: %d/%d" % (len(fails), total))
        for f in fails:
            print("  " + f)
        return 1
    print("SELF-TEST PASSED: %d/%d cases" % (total, total))
    print("  geo formats: %d, dates: %d, accessions: %d, exclusions: %d"
          % (len(GEO_CASES), len(DATE_CASES), len(ACC_CASES),
             len(EXCLUDE_CASES)))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--serotype", type=int, default=1, choices=[1, 2, 3, 4])
    ap.add_argument("--max-ambig", type=float, default=0.30,
                    help="drop sequences with more than this fraction of "
                         "non-ACGT bases (default 0.30)")
    ap.add_argument("--drop-unresolved-geo", action="store_true",
                    help="drop instead of flag records with no resolvable state")
    ap.add_argument("--drop-year-only", action="store_true",
                    help="drop instead of range-encode year-only dates")
    ap.add_argument("--stem", default=None,
                    help="input basename in 00_data/ (default denv<N>_br)")
    ap.add_argument("--outdir", default="01_metadata",
                    help="output directory under project root")
    ap.add_argument("--self-test", action="store_true",
                    help="run unit tests against fabricated rows and exit")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(self_test())
    curate(a.root, a.serotype, a.max_ambig,
           a.drop_unresolved_geo, a.drop_year_only, a.stem, a.outdir)


if __name__ == "__main__":
    main()
