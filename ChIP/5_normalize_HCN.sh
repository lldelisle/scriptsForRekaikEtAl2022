#!/bin/bash

# For H3K27ac and H3K27me3,
# in order to reduce the impact of variations in gastruloid growth speed, 
# the profile of mutant and corresponding wild-type (Figure 7 and S10)
# were corrected using the time-course experiment 
# and the profile on HoxA, HoxB and HoxC clusters as a guide.

pathWithChIP="/scratch/ldelisle/GEO_Hocine/toGEO/ChIP/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"
workingDirectory="/scratch/ldelisle/GEO_Hocine/ChIP_HCN/"

mkdir -p $workingDirectory
cd $workingDirectory

# First we generate a matrix where a profile at each hour is extrapolated from
# existing profiles (simple linear regression)
# At the same time a film is performed

# Here is the ini file:
ini_file="Film1.ini"
echo "[spacer]
[x-axis]
where=top

[spacer]
height = 0.5

[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
height = 0.4

[H3K27ac]
file=current.bedgraph
title = H3K27ac
height = 2
color = #55cf21
min_value = 0
max_value = 200
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true

[spacer]
height = 0.5

[genes]
file= ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 1.5
title = genes
fontsize = 10
prefered_name = gene_name
style = flybase
display = interleaved
fontstyle = italic
" > ${ini_file}

# First for acethylation
bigwigs=""
times=""
for time in {72..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_H3K27ac_reptc_Normalized.bigwig ]; then
        bigwigs="${bigwigs} ${pathWithChIP}/wt_${time}h_H3K27ac_reptc_Normalized.bigwig"
        times="${times} $time"
    fi
done

python ${pathWithGitHub}/scripts/BigWigMatrixGen.py \
    --bigwig $bigwigs --timepoints $times \
    --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
    --newTimeInterval 1 --prefix "H3K27ac_" \
    --ini ${ini_file}

mkdir test_film
convert *__[4-9]*.png *__1*.png test_film/Movie1.gif
ffmpeg -i test_film/Movie1.gif -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p" test_film/Movie1.mp4

# For the real movie provided with the publication,
# each png frame has been modified on ImageJ and then exported as avi,
# then converted to mp4 with ffmpeg

# Then for trimethylation
# No film is generated
bigwigs=""
times=""
for time in {72..144..24}; do
    if [ -e ${pathWithChIP}/wt_${time}h_H3K27me3_reptc_Normalized.bigwig ]; then
        bigwigs="${bigwigs} ${pathWithChIP}/wt_${time}h_H3K27me3_reptc_Normalized.bigwig"
        times="${times} $time"
    fi
done

python ${pathWithGitHub}/scripts/BigWigMatrixGen.py \
    --bigwig $bigwigs --timepoints $times \
    --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
    --newTimeInterval 1 --prefix "H3K27me3_"

# Then for each profile we find the time-point
# which seems to be the closer from the profile
# We use HoxA
timePoints='timePoints_H3K27ac.txt'
if [ -e $timePoints ]; then
    rm $timePoints
fi
for bw in ${pathWithChIP}/*h_H3K27ac_rep[1-2]_Normalized.bigwig; do
    echo -n "$(basename $bw _Normalized.bigwig) " >> $timePoints
    python ${pathWithGitHub}/scripts/TimePointFinder.py \
        --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
        --bigwig ${bw} \
        --matrix H3K27ac_HoxAll_1h_timepoint.txt --normChr chr6 >> ${timePoints}
done

timePoints='timePoints_H3K27me3.txt'
if [ -e $timePoints ]; then
    rm $timePoints
fi
for bw in ${pathWithChIP}/*h_H3K27me3_rep[1-2]_Normalized.bigwig; do
    echo -n "$(basename $bw _Normalized.bigwig) " >> $timePoints
    python ${pathWithGitHub}/scripts/TimePointFinder.py \
        --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
        --bigwig ${bw} \
        --matrix H3K27me3_HoxAll_1h_timepoint.txt --normChr chr6 >> ${timePoints}
done

# Compute the normalization:
while read l; do
    sample=$(echo $l | awk '{print $1}')
    computed=$(echo $l | awk '{print $2}')
    target=$(echo $sample | awk -F "_" '{gsub("h", "", $2); print $2}')
    echo $l
    echo "$sample $computed $target"
    python ${pathWithGitHub}/scripts/HoxClustersNormalization.py \
        --computedTimePoint $computed --targetTimePoint $target \
        --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
        --bigwig ${pathWithChIP}/${sample}_Normalized.bigwig \
        --matrix H3K27ac_HoxAll_1h_timepoint.txt \
        --normWindowSize 50
done < timePoints_H3K27ac.txt

while read l; do
    sample=$(echo $l | awk '{print $1}')
    computed=$(echo $l | awk '{print $2}')
    target=$(echo $sample | awk -F "_" '{gsub("h", "", $2); print $2}')
    echo $l
    echo "$sample $computed $target"
    python ${pathWithGitHub}/scripts/HoxClustersNormalization.py \
        --computedTimePoint $computed --targetTimePoint $target \
        --bed ${pathWithGitHub}/annotations/HoxAll_10bin.bed \
        --bigwig ${pathWithChIP}/${sample}_Normalized.bigwig \
        --matrix H3K27me3_HoxAll_1h_timepoint.txt \
        --normWindowSize 50
done < timePoints_H3K27me3.txt


# We perform average of replicates:
for rep1 in ${pathWithChIP}/*rep1_Normalized_HCN*.bedgraph; do
    join <(awk '{print $1":"$2"-"$3"\t"$4}' ${rep1}) \
        <(awk '{print $1":"$2"-"$3"\t"$4}' ${rep1/_rep1/_rep2}) \
        | awk '{gsub(":", "\t", $1); gsub("-","\t", $1); print $1"\t"($2+$3)/2}' \
        > ${rep1/_rep1/_bothrep}
done

# We perform substraction of control:
for mut in $(ls $pathWithChIP/*both*HCN* | grep -v wt); do
    wt=$pathWithChIP/$(basename $mut | awk -F "_" -v OFS="_" '{$1="wt"; print}')
    join <(awk '{print $1":"$2"-"$3"\t"$4}' ${mut}) \
        <(awk '{print $1":"$2"-"$3"\t"$4}' ${wt}) \
        | awk '{gsub(":", "\t", $1); gsub("-","\t", $1); print $1"\t"($2-$3)}' \
        > ${mut/.bedgraph/_minuswt.bedgraph}
done
