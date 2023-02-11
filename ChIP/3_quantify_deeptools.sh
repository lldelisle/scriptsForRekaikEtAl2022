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

# Figure S1
cat ${pathWithGitHub}/annotations/TADs.bed \
    ${pathWithGitHub}/annotations/subTADs.bed \
    ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_TDOM.bed \
    ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_CDOM.bed | cut -f1-4 > regions_S1.bed

multiBigwigSummary BED-file -b ${pathWithChIP}/*_H3K27ac_reptc_Normalized.bigwig \
    -o FigS1.npz --BED regions_S1.bed \
    --outRawCounts FigS1.txt

# Figure S2B
multiBigwigSummary BED-file -b ${pathWithChIP}/*_PolII_Normalized.bigwig \
    -o FigS2B.npz --BED ${pathWithGitHub}/annotations/HoxDdiv10.bed \
    --outRawCounts FigS2B.txt

# Figure S4A (only use rep1)
cp ${pathWithGitHub}/annotations/RAD21_HoxD_intCBS_regions.bed regions_S4A.bed
multiBigwigSummary BED-file -b $(ls ${pathWithChIP}/wt_*_RAD21_rep1_Normalized.bigwig | grep -v wt_48h_) \
    -o FigS4A.npz --BED regions_S4A.bed \
    --outRawCounts FigS4A.txt

# Figure S4B
cat ${pathWithGitHub}/annotations/RAD21_HoxD_CBS_regions.bed | grep -v [TC]D \
    > regions_S4B.bed
multiBigwigSummary BED-file -b ${pathWithChIP}/wt_*_RAD21*_Normalized.bigwig \
    -o FigS4B.npz --BED regions_S4B.bed \
    --outRawCounts FigS4B.txt

# Figure S5B
grep [TC]D ${pathWithGitHub}/annotations/RAD21_HoxD_CBS_regions.bed > regions_S5B.bed
multiBigwigSummary BED-file -b ${pathWithChIP}/wt_*_RAD21*_Normalized.bigwig \
    -o FigS5B.npz --BED regions_S5B.bed \
    --outRawCounts FigS5B.txt

# Figure S11C
multiBigwigSummary BED-file -b ${pathWithChIP}/*96h_RAD21*_Normalized.bigwig \
    -o FigS11C.npz --BED regions_S4B.bed \
    --outRawCounts FigS11C.txt

mkdir -p ${pathWithGitHub}/ChIP/quantifs/
cp Fig*.txt ${pathWithGitHub}/ChIP/quantifs/
cp regions*.bed ${pathWithGitHub}/ChIP/quantifs/

# Plots
Rscript ${pathWithGitHub}/scripts/plot_ChIP_quantif.R $PWD/ ${pathWithGitHub}/ChIP/plots/
