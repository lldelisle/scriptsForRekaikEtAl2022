#!/bin/bash

pathWithCHiC="/scratch/ldelisle/GEO_Hocine/toGEO/CHiC/"

cd $pathWithCHiC

# Merge replicates
# Use HiCExplorer version 3.7.2
my_sizes=${pathWithCHiC}/../../mm10.chrom.sizes
for f in *_rep3*.cool; do
  if [ ! -e ${f/rep3/merge} ]; then
    tmpdir=$(mktemp -d)
    echo $tmpdir
    hicConvertFormat --matrices $f --inputFormat cool --load_raw_values --outputFormat cool -o ${tmpdir}/$f
    echo "First matrix extracted"
    hicConvertFormat --matrices ${f/rep3/rep2} --inputFormat cool --load_raw_values --outputFormat cool -o ${tmpdir}/${f/rep3/rep2}
    echo "Second matrix extracted"
    hicConvertFormat --matrices ${f/_rep3/} --inputFormat cool --load_raw_values --outputFormat cool -o ${tmpdir}/${f/_rep3/_tmp}
    echo "Third matrix extracted"
    # Reorder bins for the rep1:
    bin=$(echo $f | awk -F "_" '{gsub("kb", "", $5);print $5}')
    cooler dump --join ${tmpdir}/${f/_rep3/_tmp}  | cooler load --format bg2 "${my_sizes}:${bin}000" - ${tmpdir}/${f/_rep3/}
    echo "Third matrix reordered"
    hicSumMatrices --matrices ${tmpdir}/$f ${tmpdir}/${f/rep3/rep2} ${tmpdir}/${f/_rep3/} -o ${f/rep3/merge}
    cooler balance --cis-only -f ${f/rep3/merge}
  fi
done
