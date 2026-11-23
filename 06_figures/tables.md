# Report Tables

## Table 1: Record Counts After Each Curation Step

| Step | Removed | Retained |
|---|---:|---:|
| Downloaded from NCBI Virus — 4261 FASTA records | — | 4261 |
| Joined FASTA to metadata on accession — version suffix stripped both sides | — | 4261 |
| Organism name matches DENV-1 — scoping step, not attrition: the download covers all four serotypes | 1595/4261 (37.4%) | 2666/2666 (100.0%) |
| Ambiguity filter — >30% non-ACGT removed; finer QC left to Nextclade | 1/2666 (0.0%) | 2665/2666 (100.0%) |
| Collection date parseable | 1/2665 (0.0%) | 2664/2666 (99.9%) |
| Year-only collection dates retained as TreeTime ranges — 56 records | — | 2664/2666 (99.9%) |
| Records without a resolvable state retained and flagged — 119 records | — | 2664/2666 (99.9%) |
| Collapsed identical sequences | 22/2664 (0.8%) | 2642/2666 (99.1%) |
| Collapsed repeated isolate names | 12/2642 (0.5%) | 2630/2666 (98.6%) |
| **Final curated set** (`curated.fasta`) | — | 2630/2666 (98.6%) |

## Table 2: Lineage Summary

| Lineage | Sequences (n) | Share (%) | Year Range | States (n) |
|---|---:|---:|---|---:|
| `1III_A.2` | 1 | 0.04% | 2019 | 0 |
| `1V` | 5 | 0.20% | 1986–1989 | 1 |
| `1V_A` | 60 | 2.45% | 2000–2023 | 6 |
| `1V_A.1` | 30 | 1.23% | 2013–2022 | 5 |
| `1V_A.2` | 10 | 0.41% | 2021–2022 | 2 |
| `1V_B` | 1 | 0.04% | 2019 | 0 |
| `1V_C` | 1 | 0.04% | 2020 | 1 |
| `1V_D.1.1` | 507 | 20.72% | 2020–2024 | 13 |
| `1V_D.1.4` | 90 | 3.68% | 2021 | 2 |
| `1V_D.2` | 14 | 0.57% | 2021–2024 | 4 |
| `1V_E` | 9 | 0.37% | 2010–2017 | 1 |
| `1V_E.1` | 1572 | 64.24% | 2015–2025 | 21 |
| `1V_E.2` | 21 | 0.86% | 2020–2022 | 2 |
| `1V_F` | 40 | 1.63% | 2008–2024 | 11 |
| `1V_F.1` | 38 | 1.55% | 2019–2022 | 2 |
| `1V_F.2` | 32 | 1.31% | 2021–2022 | 2 |
| `1V_H` | 16 | 0.65% | 1996–2020 | 3 |

## Table 3: Clock Outliers Flagged by TreeTime (`--clock-filter 4`)

| Accession | Lineage | State | Stated Date |
|---|---|---|---|
| `PP234964` | `1III_A.2` | nan | 2019-07-01 |
| `NODE_0000023` | `nan` | nan | nan |
| `PV205900` | `1V_A` | ES | 2023-05-02 |
| `NODE_0000024` | `nan` | nan | nan |
| `NODE_0000025` | `nan` | nan | nan |
| `PV205870` | `1V_A` | ES | 2023-01-17 |
| `PV205875` | `1V_A` | ES | 2023-01-16 |
| `NODE_0000026` | `nan` | nan | nan |
| `PV205911` | `1V_A` | ES | 2023-05-05 |
| `NODE_0000027` | `nan` | nan | nan |
| `PV205883` | `1V_A` | ES | 2023-04-24 |
| `PV205889` | `1V_A` | ES | 2023-04-24 |
| `OR258460` | `1V_A.1` | MS | 2022-04-06 |
| `PV682577` | `1V_H` | nan | 2020 |
| `MW377899` | `1V_F` | AP | 2014-12-10 |
| `OR518260` | `1V_E.1` | TO | 2022-02-10 |
| `NODE_0001704` | `nan` | nan | nan |
| `OR258487` | `1V_E.1` | MS | 2022-04-06 |
| `OR258539` | `1V_E.1` | MS | 2022-04-16 |
