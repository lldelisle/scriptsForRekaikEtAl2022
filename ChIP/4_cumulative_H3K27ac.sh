#!/bin/bash

pathWithChIP="/scratch/ldelisle/GEO_Hocine/toGEO/ChIP/"
pathWithPlots="/scratch/ldelisle/GEO_Hocine/plots/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"

mkdir -p $pathWithPlots

python ${pathWithGitHub}/scripts/CumulativeBigWig.py \
    --bigwig ${pathWithChIP}/*_H3K27ac_reptc_Normalized.bigwig \
    --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
    --output ${pathWithPlots}/wt_H3K27ac_reptc_cum.bedgraph
