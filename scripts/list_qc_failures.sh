#!/bin/bash

cd ~/Documents/denv-project

awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++){if($i=="seqName")s=i; if($i=="qc.overallStatus")q=i}; next}
            $q=="bad" {print $s}' 02_qc/out/nextclade.tsv > 02_qc/qc_bad.txt

wc -l 02_qc/qc_bad.txt