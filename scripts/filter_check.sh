#!/bin/bash

awk -F',' 'NR==FNR {bad[$1]; next} FNR==1 || !($1 in bad)' \
    02_qc/qc_bad.txt 01_metadata/dates.csv > 02_qc/dates_qc_passed.csv

wc -l 02_qc/dates_qc_passed.csv

diff <(grep '^>' 02_qc/aligned_qc_passed.fasta | tr -d '>' | sort) \
     <(tail -n +2 02_qc/dates_qc_passed.csv | cut -d, -f1 | sort) && echo "names match"