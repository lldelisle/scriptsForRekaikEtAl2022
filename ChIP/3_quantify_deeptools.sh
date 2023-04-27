#!/bin/bash

# deeptools version 3.0.0

pathWithChIP="/scratch/ldelisle/GEO_Hocine/toGEO/ChIP/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"

mkdir -p /scratch/ldelisle/GEO_Hocine/quantifs/
cd /scratch/ldelisle/GEO_Hocine/quantifs/

# Figure 1D
multiBigwigSummary BED-file -b ${pathWithChIP}/*_H3K27ac_reptc_Normalized.bigwig \
    -o Fig1D.npz --BED ${pathWithGitHub}/annotations/HoxDdiv10.bed \
    --outRawCounts Fig1D.txt

# Figure Extended Data 1
cat ${pathWithGitHub}/annotations/TADs.bed \
    ${pathWithGitHub}/annotations/subTADs.bed \
    ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_TDOM.bed \
    ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_CDOM.bed | cut -f1-4 > regions_ExtD1.bed

multiBigwigSummary BED-file -b ${pathWithChIP}/*_H3K27ac_reptc_Normalized.bigwig \
    -o ExtDFig1.npz --BED regions_ExtD1.bed \
    --outRawCounts ExtDFig1.txt

# Supplementary Figure 1B
multiBigwigSummary BED-file -b ${pathWithChIP}/*_PolII_Normalized.bigwig \
    -o SupFig1B.npz --BED ${pathWithGitHub}/annotations/HoxDdiv10.bed \
    --outRawCounts SupFig1B.txt

# Extended Data Figure 2A (only use rep1)
cp ${pathWithGitHub}/annotations/RAD21_HoxD_intCBS_regions.bed regions_ExtD2A.bed
multiBigwigSummary BED-file -b $(ls ${pathWithChIP}/wt_*_RAD21_rep1_Normalized.bigwig | grep -v wt_48h_) \
    -o ExtDFig2A.npz --BED regions_ExtD2A.bed \
    --outRawCounts ExtDFig2A.txt

# Extended Data Figure 2B
cat ${pathWithGitHub}/annotations/RAD21_HoxD_CBS_regions.bed | grep -v [TC]D \
    > regions_ExtD2B.bed
multiBigwigSummary BED-file -b ${pathWithChIP}/wt_*_RAD21*_Normalized.bigwig \
    -o ExtDFig2B.npz --BED regions_ExtD2B.bed \
    --outRawCounts ExtDFig2B.txt

# Extended Data Figure 3B
grep [TC]D ${pathWithGitHub}/annotations/RAD21_HoxD_CBS_regions.bed > regions_ExtD3B.bed
multiBigwigSummary BED-file -b ${pathWithChIP}/wt_*_RAD21*_Normalized.bigwig \
    -o ExtDFig3B.npz --BED regions_ExtD3B.bed \
    --outRawCounts ExtDFig3B.txt

# Extended Data Figure 8C
multiBigwigSummary BED-file -b ${pathWithChIP}/*96h_RAD21*_Normalized.bigwig \
    -o ExtDFig8C.npz --BED regions_ExtD2B.bed \
    --outRawCounts ExtDFig8C.txt

mkdir -p ${pathWithGitHub}/ChIP/quantifs/
cp Fig*.txt ${pathWithGitHub}/ChIP/quantifs/
cp regions*.bed ${pathWithGitHub}/ChIP/quantifs/

# Plots
Rscript ${pathWithGitHub}/scripts/plot_ChIP_quantif.R $PWD/ ${pathWithGitHub}/ChIP/plots/
