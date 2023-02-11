#!/bin/bash

pathWithCHiC="/scratch/ldelisle/GEO_Hocine/toGEO/CHiC/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"

mkdir -p /scratch/ldelisle/GEO_Hocine/cHi-C_movie/
cd /scratch/ldelisle/GEO_Hocine/cHi-C_movie/

# First we generate a matrix where a matrix at each hour is extrapolated from
# existing matrices (simple linear regression)
# At the same time a film is performed

# Here is the ini file:
ini_file="Film2.ini"
echo "[spacer]
[x-axis]
where=top

[spacer]
height = 0.5

[wt_CHiC]
file = current.cool
depth = 1300000
show_masked_bins = false

[CBS1_highlight]
file = ${pathWithGitHub}/annotations/CBS1_highlight.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dotted
color = white
line_width = 2

[spacer]

[wt_168h_CTCF]
file = /scratch/ldelisle/GEO_Hocine/toGEO/ChIP//wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig

[spacer]
height = 0.1

[annotations]
file = ${pathWithGitHub}/annotations/TADs.bed
height = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb
labels = false

[labels]
file = ${pathWithGitHub}/annotations/TADs_label.bed
overlay_previous = yes
color = none
border_color = none
display = collapsed
fontsize = 20

[spacer]
height = 0.2

[annotations2]
file = ${pathWithGitHub}/annotations/subTADs.bed
height = 0.2
display = collapsed
border_color = white
color = bed_rgb
labels = false

[labels]
file = ${pathWithGitHub}/annotations/subTADs_label.bed
color = none
border_color = none
display = collapsed
fontsize = 20
height = 1.5
" > ${ini_file}

# Use 5kb wt
cools=""
times=""
for time in {48..168..24}; do
    if [ -e ${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool ]; then
        cools="${cools} ${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool"
        times="${times} $time"
    elif [ -e ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool ]; then
        cools="${cools} ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool"
        times="${times} $time"
    else
        echo "Did not find wt_${time}h_CHiC_5kb_Balanced_mm10.cool"
    fi
done

python ${pathWithGitHub}/scripts/CoolMatrixGen.py \
    --cool $cools --timepoints $times \
    --region 'chr2:73,000,000-77,000,000' \
    --newTimeInterval 1 --prefix "wt_5kb_" \
    --ini ${ini_file} --plottedRegion "chr2:74,179,681-75,653,876"

mkdir test_film
convert *__[4-9]*.png *__1*.png test_film/Movie2.gif
ffmpeg -i test_film/Movie2.gif -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p" test_film/Movie2.mp4

# For the real movie provided with the publication,
# each png frame has been modified on ImageJ and then exported as avi,
# then converted to mp4 with ffmpeg
