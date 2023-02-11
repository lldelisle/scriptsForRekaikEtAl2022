# In the original version of the figures, pgt version 3.6 was used and images were modified with illustrator
# In order to produce figures closer with fewer modifications with illustrator I used 3.7
pathWithChIP="/scratch/ldelisle/GEO_Hocine/toGEO/ChIP/"
pathWithCHiC="/scratch/ldelisle/GEO_Hocine/toGEO/CHiC/"
pathWithv4C="/scratch/ldelisle/GEO_Hocine/v4C/"
pathWithHiChIP="/scratch/ldelisle/GEO_Hocine/toGEO/HiChIP/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"
pathWithPlots="/scratch/ldelisle/GEO_Hocine/plots/"
mkdir -p ${pathWithPlots}
cd ${pathWithPlots}

wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -nc
gunzip -k Mus_musculus.GRCm38.102.gtf.gz
transcripts=$(cat ${pathWithGitHub}/annotations/selected_transcripts.txt)
if [ -e ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf ]; then
    rm ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
fi
for tr in $transcripts; do
    grep "$tr" Mus_musculus.GRCm38.102.gtf >> ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
done

cat ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf | awk -F "\t" -v OFS="\t" '
$3=="transcript"{
    i=match($9, "gene_name")
    j=match(substr($9,i), ";")
    gene_name=substr($9,i + 14, j - 14 - 2)
    pos=$4
    if (gene_name == "d3") {
        pos=74746000
    } else if (gene_name == "d4"){
        pos=74727000
    } else if (gene_name == "2"){
        gene_name=substr($9,i + 11, j - 11 - 2)
    }
    print $1"\t"pos"\t"pos+1"\t"gene_name
}' > Hox_left.bed

# Get the script to switch annotations:
wget "https://raw.githubusercontent.com/lldelisle/scriptsForRodriguezCarballoEtAl2020/45972dc006649c9bc018d2a0411f40f662d2b8a4/scripts/shiftAnnotFunctions.R"
wget "https://raw.githubusercontent.com/lldelisle/scriptsForRodriguezCarballoEtAl2020/45972dc006649c9bc018d2a0411f40f662d2b8a4/scripts/shiftAnnot_splitIfOV.R"

# Write the "BR" file to shift the annotations in the mutant genome
echo -e "genome\tbr1\tbr2\ttg1" > br.txt
# a 900bp cassette was inserted between bases mm10:74,716,919 and 74,717,063
echo -e "Ins(2xCBS-d4d8)\t74716920\t74717063\t900" >> br.txt
# Del(d1-d4) mutant, the deletion includes bases between mm10:74,719,771 and 74,765,994
echo -e "Del(d1-d4)\t74719771\t74765994\t0" >> br.txt
# Del(sub-TAD1) chr2 contains the deletion of mm10:74768588 to 75133800. 
echo -e "Del(sub-TAD1)\t74768588\t75133800\t0" >> br.txt

# Shift exons and CDS (and add 'chr' to be UCSC style and be able to shift gtf)
awk -v OFS="\t" '$3=="exon"||$3=="CDS"{print "chr"$0}' ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf > temp.gtf
Rscript shiftAnnot_splitIfOV.R br.txt temp.gtf 1 4 5 Hox_single_transcripts.gtf

# Shift other annotations
Rscript shiftAnnot_splitIfOV.R br.txt ${pathWithGitHub}/annotations/TADs.bed 1 2 3 TADs.bed
Rscript shiftAnnot_splitIfOV.R br.txt ${pathWithGitHub}/annotations/subTADs.bed 1 2 3 subTADs.bed
Rscript shiftAnnot_splitIfOV.R br.txt ${pathWithGitHub}/annotations/TADs_label.bed 1 2 3 TADs_label.bed
Rscript shiftAnnot_splitIfOV.R br.txt ${pathWithGitHub}/annotations/subTADs_label.bed 1 2 3 subTADs_label.bed
Rscript shiftAnnot_splitIfOV.R br.txt ${pathWithGitHub}/annotations/CS3840.bed 1 2 3 CS3840.bed

# Figure 1 core:
ini_file="Fig1core.ini"
echo "" > ${ini_file}
for time in {48..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_H3K27ac_reptc_Normalized.bigwig ]; then
        echo "[wt_${time}h_H3k27ac_reptc]
file = ${pathWithChIP}/wt_${time}h_H3K27ac_reptc_Normalized.bigwig
title = ${time}h_H3k27ac
height = 3
color = #55cf21
min_value = 0
max_value = 150 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done


# Figure 1B
ini_file="Fig1B.ini"
echo "[x-axis]
where = top

[spacer]
height = 0.5
" > ${ini_file}
cat Fig1core.ini >> ${ini_file}
echo "[HoxD highlight]
file = ${pathWithGitHub}/annotations/HoxD.bed
type = vhighlight
color = #00aeef
alpha = 0.1

[C-DOM highlight]
file = ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_CDOM.bed
type = vhighlight
color = #ef4036
alpha = 0.1

[T-DOM highlight]
file = ${pathWithGitHub}/annotations/H3K27ac_HoxD_locations_TDOM.bed
type = vhighlight
color = #8cc63f
alpha = 0.1

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
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,058,812-75,485,391 -o ${ini_file/.ini/.pdf}

# Figure 1C
ini_file="Fig1C.ini"
cp Fig1core.ini ${ini_file}
sed -i 's/max_value = 150/max_value = 200/' ${ini_file}
echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 2
title = genes 102 single transcripts
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.025

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,667,374-74,767,842 -o ${ini_file/.ini/.pdf}

# Figure 1D bottom
ini_file="Fig1Dbottom.ini"
echo "[bins]
file = ${pathWithGitHub}/annotations/HoxDdiv10.bed
height = 1
display = interleaved
color = black
labels_in_margin = true

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 2
title = genes 102 single transcripts
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" > ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,667,374-74,767,842 -o ${ini_file/.ini/.pdf}


# Figure 2A
ini_file="Fig2A.ini"

echo "[CTCF]
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.5
file_type = bed
labels = false
color = bed_rgb

[wt_72h_CTCF]
file = ${pathWithChIP}/wt_72h_CTCF.bigwig
title = 72h CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 69
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]

[72h-168h H3K27ac]
file = ${pathWithPlots}/wt_H3K27ac_reptc_cum.bedgraph
title = 72h-168h_H3k27ac
height = 3
color = #55cf21
min_value = 0
show_data_range = true

[spacer]
height = 0.1

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" > ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,650,000-74,767,842 -o ${ini_file/.ini/.pdf} --width 25

# Figure 2B
ini_file="Fig2B.ini"
echo "[CTCF]
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.5
file_type = bed
labels = false
color = bed_rgb

[spacer]
height = 0.2
" > ${ini_file}
for time in {72..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_RAD21_rep1_Normalized.bigwig ]; then
        echo "[wt_${time}h_RAD21_rep1]
file = ${pathWithChIP}/wt_${time}h_RAD21_rep1_Normalized.bigwig
title = ${time}h_RAD21_rep1
height = 2
color = #f562d5
min_value = 0
max_value = 500 
number_of_bins = 2000
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done
echo "[spacer]
height = 0.2

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 1
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,667,374-74,767,842 -o ${ini_file/.ini/.pdf}

# Figure 2C
ini_file="Fig2C.ini"
echo "[CTCF]
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.5
file_type = bed
labels = false
color = bed_rgb

[spacer]
height = 0.2
" > ${ini_file}
for time in {72..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_NIPBL_Normalized.bigwig ]; then
        echo "[wt_${time}h_NIPBL]
file = ${pathWithChIP}/wt_${time}h_NIPBL_Normalized.bigwig
title = ${time}h_NIPBL
height = 2
color = #204D48
min_value = 0
max_value = 160
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done
echo "[spacer]
height = 0.2

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 1
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,667,374-74,767,842 -o ${ini_file/.ini/.pdf}

# Figure 3A
ini_file="Fig3A.ini"
echo "[wt_48h_CHiC]
file = ${pathWithCHiC}/wt_48h_CHiC_merge_5kb_Balanced_mm10.cool
title = wt 48h CHiC 5kb balanced
depth = 1000000
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

pgt --tracks ${ini_file} --region chr2:73,900,037-75,621,560 -o ${ini_file/.ini/.pdf}
vmin48=0.0002769132964053131
vmax48=0.019123641300184402

# Figure 3B
ini_file="Fig3B.ini"
echo "" > ${ini_file}
for time in {72..144..24}; do
    if [ -e ${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool
    elif [ -e ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool
    else
        cool_file=""
    fi
    if [ ! -z $cool_file ]; then
        echo "[wt_${time}h_CHiC]
file = ${cool_file}
title = ${time}h_CHiC_5kb_Balanced
depth = 600000
show_masked_bins = false

[CBS1_highlight]
file = ${pathWithGitHub}/annotations/CBS1_highlight.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dotted
color = white
line_width = 1
" >> ${ini_file}
        if [ "$time" = "72" ]; then
            echo "[zoom]
file = ${pathWithGitHub}/annotations/HoxDvsCS3840.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dashed
color = white
line_width = 1
" >> ${ini_file}
        fi
echo "[spacer]
height = 0.1
" >> ${ini_file}
    fi
done
echo "[wt_168h_CTCF]
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
orientation = inverted
file_type = bigwig
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf} --width 20


ini_file="Fig3B_notused.ini"
echo "" > ${ini_file}
for time in 168; do
    if [ -e ${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_merge_5kb_Balanced_mm10.cool
    elif [ -e ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool
    else
        cool_file=""
    fi
    if [ ! -z $cool_file ]; then
        echo "[wt_${time}h_CHiC]
file = ${cool_file}
title = ${time}h_CHiC_5kb_Balanced
depth = 600000
show_masked_bins = false

[CBS1_highlight]
file = ${pathWithGitHub}/annotations/CBS1_highlight.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dotted
color = white
line_width = 1
" >> ${ini_file}
        if [ "$time" = "72" ]; then
            echo "[zoom]
file = ${pathWithGitHub}/annotations/HoxDvsCS3840.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dashed
color = white
line_width = 1
" >> ${ini_file}
        fi
echo "[spacer]
height = 0.1
" >> ${ini_file}
    fi
done
echo "[wt_168h_CTCF]
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
orientation = inverted
file_type = bigwig
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf} --width 20

times=(48 72 96 120 144 168)
vmins=($vmin48 0.001095951913515006 0.0011560076281525809 0.0011909313825654738 0.0012646813965259626 0.0011796867803342272)
vmaxs=($vmax48 0.01642239782107665 0.018603166415317677 0.017884795102223634 0.015623469750531742 0.016004993782163884)

# Figure 3C
ini_file="Fig3C.ini"
echo "" > ${ini_file}
for i in "${!times[@]}"; do
    time=${times[$i]}
    vmin=${vmins[$i]}
    vmax=${vmaxs[$i]}
    if [ -e ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool ]; then
        echo "[wt_${time}h_CHiC]
file = ${pathWithCHiC}/wt_${time}h_CHiC_5kb_Balanced_mm10.cool
title = ${time}h_CHiC_5kb_Balanced
file_type = hic_matrix_square
show_masked_bins = false
region2 = chr2:75105000-75190000
orientation = inverted
min_value = $vmin
max_value = $vmax
height = 3.5
colormap = RdYlBu_r

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
orientation = inverted
file_type = bigwig

[spacer]
height = 2
" >> ${ini_file}
    fi
done
echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,637,423-74,765,000 -o ${ini_file/.ini/.pdf} --width 10

sed "s/RdYlBu_r/Blues_r/g" ${ini_file} > ${ini_file/.ini/_forquantif.ini}
pgt --tracks ${ini_file/.ini/_forquantif.ini} --region chr2:74,637,423-74,765,000 -o ${ini_file/.ini/_forquantif.pdf} --width 10

# The CTCF track on the vertical was added manually with illustrator.
# The histogram was done with Plot profile on imageJ version 2.1.0/1.53c

# Figure 3D
ini_file="Fig3D.ini"
echo "[wt_120h_HiChIP]
file = ${pathWithHiChIP}/wt_120h_H3K27ac_HiChIP_10kb_Raw_mm10.cool
title = 120h HiChIP H3K27ac
colormap = hot
depth = 1000000
show_masked_bins = false

[spacer]
height = 0.1

[wt_120h_H3K27ac]
file = ${pathWithChIP}/wt_120h_H3K27ac_reptc_Normalized.bigwig
title = 120h_H3k27ac
height = 1.5
color = #55cf21
min_value = 0
max_value = 150
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

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

[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig
" > ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,190,037-75,621,560 -o ${ini_file/.ini/.pdf}

# Figure 5ACE
panels=('A' 'C' 'E')
genotypes=('Del(CBS1)' 'Del(CBS1-2)' 'Del(CBS2)')
for i in "${!panels[@]}"; do
    ini_file="Fig5${panels[$i]}.ini"
    geno=${genotypes[$i]}
    echo "[CTCF]
arrowhead_fraction = 0.006
display = collapsed
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.3
file_type = bed
labels = false
color = bed_rgb
line_width = 0.05

[spacer]

[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[${geno}_96h_CTCF]
file = ${pathWithChIP}/${geno}_96h_CTCF.bigwig
title = CTCF ${geno} 96h
height = 1.5
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" > ${ini_file}
    pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 10
done

# Figure 5G
ini_file="Fig5G.ini"
geno="Ins(2xCBS-d4d8)"
echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[${geno}_96h_CTCF]
file = ${pathWithChIP}/${geno}_144h_CTCF.bigwig
title = CTCF ${geno} 144h
height = 1.5
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]

[genes]
file = ${geno}/${geno}_vSplit_Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts shifted
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 10

# Figure 6C
ini_file="Fig6C.ini"
geno="Del(d1-d4)"
time=96
echo "[${geno}_${time}h_CHiC]
file = ${pathWithCHiC}/${geno}_${time}h_CHiC_5kb_Balanced_${geno}.cool
title = ${geno}_CHiC_5kb_Balanced
depth = 500000
show_masked_bins = false
min_value = 0

[HoxD-CS3840]
file = ${pathWithGitHub}/annotations/HoxDvsCS3840_Del(d1-d4).bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dotted
color = white
line_width = 1

[annotations]
file = ${geno}/${geno}_vSplit_TADs.bed
height = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb
labels = false

[labels]
file = ${geno}/${geno}_vSplit_TADs_label.bed
overlay_previous = yes
color = none
border_color = none
display = collapsed
fontsize = 20


[spacer]
height = 0.05

[annotations2]
file = ${geno}/${geno}_vSplit_subTADs.bed
height = 0.2
display = collapsed
border_color = white
color = bed_rgb
labels = false

[genes]
file = ${geno}/${geno}_vSplit_Hox_single_transcripts.gtf
overlay_previous = yes
title = genes 102 single transcripts shifted
style = flybase
prefered_name = gene_name
display = collapsed
labels = false

[labels]
file = ${geno}/${geno}_vSplit_subTADs_label.bed
color = none
border_color = none
display = collapsed
fontsize = 20
height = 1.5

[CS38-40]
file = ${geno}/${geno}_vSplit_CS3840.bed
overlay_previous = yes
color = bed_rgb
" > ${ini_file}

# Check coo
pgt --tracks ${ini_file} --region chr2:74,635,000-75,165,000 -o ${ini_file/.ini/.pdf} --width 20


# Figure 6F
ini_file="Fig6F.ini"
geno="Del(sub-TAD1)"
time=96
echo "[${geno}_${time}h_CHiC]
file = ${pathWithCHiC}/${geno}_${time}h_CHiC_5kb_Balanced_${geno}.cool
title = ${geno}_CHiC_5kb_Balanced
depth = 500000
show_masked_bins = false
min_value = 0

[annotations]
file = ${geno}/${geno}_vSplit_TADs.bed
height = 1
display = collapsed
color = bed_rgb
border_color = bed_rgb
labels = false

[labels]
file = ${geno}/${geno}_vSplit_TADs_label.bed
overlay_previous = yes
color = none
border_color = none
display = collapsed
fontsize = 20


[spacer]
height = 0.05

[annotations2]
file = ${geno}/${geno}_vSplit_subTADs.bed
height = 0.2
display = collapsed
border_color = white
color = bed_rgb
labels = false

[genes]
file = ${geno}/${geno}_vSplit_Hox_single_transcripts.gtf
overlay_previous = yes
title = genes 102 single transcripts shifted
style = flybase
prefered_name = gene_name
display = collapsed
labels = false

[labels]
file = ${geno}/${geno}_vSplit_subTADs_label.bed
color = none
border_color = none
display = collapsed
fontsize = 20
height = 1.5

[CS38-40]
file = ${geno}/${geno}_vSplit_CS3840.bed
overlay_previous = yes
color = bed_rgb
" > ${ini_file}

# Check coo
pgt --tracks ${ini_file} --region chr2:74,635,000-75,165,000 -o ${ini_file/.ini/.pdf} --width 20


# Figure 7A
ini_file=Fig7A.ini
echo "[CTCF]
arrowhead_fraction = 0.006
display = collapsed
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.3
file_type = bed
labels = false
color = bed_rgb
line_width = 0.05


[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
height = 0.1

[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF Del(CBS1-5) 96h
height = 1.5
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
height = 0.1

[Del(CBS1-5)_120h_RAD21]
file = ${pathWithChIP}/Del(CBS1-5)_120h_RAD21.bigwig
title = RAD21 Del(CBS1-5) 120h
height = 1.5
color = #b849a0
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 10


# Figure 7D
ini_file="Fig7D.ini"
echo "" > ${ini_file}
for time in {48..168..24}; do
    if [ -e ${pathWithCHiC}/Del\(CBS1-5\)_${time}h_CHiC_2kb_Balanced_mm10.cool ]; then
        echo "[Del(CBS1-5)_${time}h_CHiC]
file = ${pathWithCHiC}/Del(CBS1-5)_${time}h_CHiC_2kb_Balanced_mm10.cool
title = ${time}h_CHiC_2kb_Balanced
show_masked_bins = false
depth = 150000

[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF
height = 2
color = #ca6a02
min_value = 0
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig

[spacer]
height = 0.5
" >> ${ini_file}
    fi
done
echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 1
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.01

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf}


# Figure S2A
ini_file="FigS2A.ini"
echo "" > ${ini_file}
for time in {72..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_PolII_Normalized.bigwig ]; then
        echo "[wt_${time}h_PolII]
file = ${pathWithChIP}/wt_${time}h_PolII_Normalized.bigwig
title = ${time}h_PolII
height = 1.5
color = #0394fc
min_value = 0
max_value = 200 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done

echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 25

# Figure S2D
ini_file="FigS2D.ini"
echo "" > ${ini_file}
for time in {72..144..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_pSer2PolII_Normalized.bigwig ]; then
        echo "[wt_${time}h_pSer2PolII]
file = ${pathWithChIP}/wt_${time}h_pSer2PolII_Normalized.bigwig
title = ${time}h_pSer2PolII
height = 1.5
color = #1a6985
min_value = 0
max_value = 100 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done
echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 25

# Figure S3
ini_file="FigS3.ini"
echo "[x-axis]
where = top
" > ${ini_file}
for time in {72..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_PolII_Normalized.bigwig ]; then
        echo "[wt_${time}h_PolII]
file = ${pathWithChIP}/wt_${time}h_PolII_Normalized.bigwig
title = ${time}h_PolII
height = 1.5
color = #0394fc
min_value = 0
max_value = 200 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done

echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16

[spacer]
height = 2
" >> ${ini_file}
for time in {72..144..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_pSer2PolII_Normalized.bigwig ]; then
        echo "[wt_${time}h_pSer2PolII]
file = ${pathWithChIP}/wt_${time}h_pSer2PolII_Normalized.bigwig
title = ${time}h_pSer2PolII
height = 1.5
color = #1a6985
min_value = 0
max_value = 100 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done
echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --BED ${pathWithGitHub}/annotations/HoxAB_plot_regions.bed -o ${ini_file/.ini/.pdf} --width 25

# Figure S4A top
ini_file="FigS4Atop.ini"
echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16

[CTCF]
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.5
file_type = bed
labels = false
color = bed_rgb

[quantif]
file = ${pathWithGitHub}/annotations/RAD21_HoxD_intCBS_regions.bed
display = collapsed
labels = false
border_color = none
" > ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf}



# Figure S5A
ini_file="FigS5A.ini"
echo "[CTCF]
file = ${pathWithGitHub}/annotations/RAD21_HoxD_CBS_regions.bed
title = CTCF-binding sites
height = 0.5
file_type = bed
labels = false
color = bed_rgb
display = collapsed

[wt_72h_CTCF]
file = ${pathWithChIP}/wt_72h_CTCF.bigwig
title = 72h CTCF
height = 2
color = #fc8403
min_value = 0
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" > ${ini_file}

for time in {48..168..12}; do
    if [ -e ${pathWithChIP}/wt_${time}h_RAD21_rep1_Normalized.bigwig ]; then
        echo "[wt_${time}h_RAD21_rep1]
file = ${pathWithChIP}/wt_${time}h_RAD21_rep1_Normalized.bigwig
title = ${time}h_RAD21_rep1
height = 2
color = #f562d5
min_value = 0
max_value = 600 
number_of_bins = 2000
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done
echo "[annotations]
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
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:73,900,037-75,621,560 -o ${ini_file/.ini/.pdf}

# Figure S6Atop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep CBS > vp_CBS.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep -P "sub-TAD1$" > v4C_subTAD1.bed

ini_file="FigS6Atop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_CBS.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_subTAD1.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}

# Figure S6Btop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep CBS > vp_CBS.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "CBS" > v4C_CBSsubTAD1.bed

ini_file="FigS6Btop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_CBS.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_CBSsubTAD1.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}

# Figure S6Ctop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep CBS > vp_CBS.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "CBS-CS38-40" > v4C_CBSCS3840.bed

ini_file="FigS6Ctop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_CBS.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_CBSCS3840.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}

# Figure S6Dtop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep "d1-d4" > vp_d1d4.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep -P "sub-TAD1$" > v4C_subTAD1.bed

ini_file="FigS6Dtop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_d1d4.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_subTAD1.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}

# Figure S6Etop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep "d9" > vp_d9.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep -P "sub-TAD1$" > v4C_subTAD1.bed

ini_file="FigS6Etop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_d9.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_subTAD1.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}


# Figure S7A
ini_file="FigS7A.ini"
echo "" > ${ini_file}
for time in {48..168..24}; do
    if [ -e ${pathWithCHiC}/wt_${time}h_CHiC_merge_2kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_merge_2kb_Balanced_mm10.cool
    elif [ -e ${pathWithCHiC}/wt_${time}h_CHiC_2kb_Balanced_mm10.cool ]; then
        cool_file=${pathWithCHiC}/wt_${time}h_CHiC_2kb_Balanced_mm10.cool
    else
        cool_file=""
    fi
    if [ ! -z $cool_file ]; then
        echo "[wt_${time}h_CHiC]
file = ${cool_file}
title = ${time}h_CHiC_2kb_Balanced
show_masked_bins = false
depth = 150000

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
orientation = inverted
file_type = bigwig

[spacer]
height = 0.5
" >> ${ini_file}
    fi
done
echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf}

# Figure S7B
ini_file="FigS7B.ini"
echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
" > ${ini_file}
for time in 48 96; do
    echo "[wt_${time}h_CHiC_d9]
file = ${pathWithv4C}/wt_${time}h_CHiC_merge_d9.bedgraph
title = wt_${time}h_v4C_d9
height = 4
number_of_bins = 1000
nans_to_zeros = true
summary_method = mean
show_data_range = true

[spacer]
height = 0.5
" >> ${ini_file}
done

echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

# Check coo
pgt --tracks ${ini_file} --region chr2:74,616,423-74,812,095 -o ${ini_file/.ini/.pdf} --width 20

# Get the first 2 CTCF:
cat ${pathWithGitHub}/annotations/CTCF.bed | tail -n 2 > CBS12.bed

# Figure S7C
ini_file="FigS7C.ini"
echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
" > ${ini_file}
for time in {72..144..24}; do
    if [ -e ${pathWithChIP}/wt_${time}h_H3K27me3_reptc_Normalized.bigwig ]; then
        echo "[wt_${time}h_H3K27me3_reptc]
file = ${pathWithChIP}/wt_${time}h_H3K27me3_reptc_Normalized.bigwig
title = ${time}h_H3K27me3_reptc
height = 4
color = #873e23
min_value = 0
max_value = 220 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done

echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16

[vlines]
file = CBS12.bed
type = vlines
" >> ${ini_file}

# Check coo
pgt --tracks ${ini_file} --region chr2:74,616,423-74,812,095 -o ${ini_file/.ini/.pdf} --width 20

# Figure S8top
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep -v CBS | grep -v "d1-d4" > vp_genes.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "Evx2-Hoxd12" > v4C_Evx2d12.bed

ini_file="FigS8top.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_genes.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_Evx2d12.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf}

# Figure S9ABC
ini_file="FigS9ABC.ini"
echo "[x-axis]
where=top

[spacer]

[CTCF]
arrowhead_fraction = 0.006
display = collapsed
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.3
file_type = bed
labels = false
color = bed_rgb
line_width = 0.05

[wt_72h_CTCF]
file = ${pathWithChIP}/wt_72h_CTCF.bigwig
title = 72h CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 69
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
" > ${ini_file}
for time in 72 96 144; do
    echo "[wt_${time}h_H3k27ac_reptc]
file = ${pathWithChIP}/wt_${time}h_H3K27ac_reptc_Normalized.bigwig
title = ${time}h_H3k27ac
height = 3
color = #55cf21
min_value = 0
max_value = 250 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
done
echo "[spacer]
" >> ${ini_file}
for time in 72 96 144; do
    echo "[wt_${time}h_H3K27me3_reptc]
file = ${pathWithChIP}/wt_${time}h_H3K27me3_reptc_Normalized.bigwig
title = ${time}h_H3K27me3_reptc
height = 4
color = #873e23
min_value = 0
max_value = 220 
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
done
echo "[spacer]
" >> ${ini_file}
for time in 72 96 144; do
    echo "[wt_${time}h_RAD21_rep1]
file = ${pathWithChIP}/wt_${time}h_RAD21_rep1_Normalized.bigwig
title = ${time}h_RAD21_rep1
height = 2
color = #f562d5
min_value = 0
max_value = 500 
number_of_bins = 2000
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bigwig
" >> ${ini_file}
done
echo "[spacer]
" >> ${ini_file}
for time in 72 96 144; do
    echo "[wt_${time}h_NIPBL]
file = ${pathWithChIP}/wt_${time}h_NIPBL_Normalized.bigwig
title = ${time}h_NIPBL
height = 2
color = #204D48
min_value = 0
max_value = 160
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
done
echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 7
title = genes 102 single transcripts
fontsize = 10
style = flybase
gene rows = 10
prefered_name = gene_name
" >> ${ini_file}

pyGenomeTracks --tracks ${ini_file} --region chr6:52147123-52277553 --outFileName FigS9A.pdf
pyGenomeTracks --tracks ${ini_file} --region chr11:96183548-96376004 --outFileName FigS9B.pdf
pyGenomeTracks --tracks ${ini_file} --region chr15:102,919,005-103,068,700 --outFileName FigS9C.pdf

# Figure S9D
ini_file="FigS9D.ini"
echo "[x-axis]
where = top

[spacer]

[wt_120h_HiChIP]
file = ${pathWithHiChIP}/wt_120h_H3K27ac_HiChIP_10kb_Raw_mm10.cool
title = 120h HiChIP H3K27ac
colormap = hot
depth = 1200000
show_masked_bins = false

[spacer]
height = 0.1

[wt_120h_H3K27ac]
file = ${pathWithChIP}/wt_120h_H3K27ac_reptc_Normalized.bigwig
title = 120h_H3k27ac
height = 1.5
color = #55cf21
min_value = 0
max_value = 150
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 2000
orientation = inverted
show_data_range = true
file_type = bigwig

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 7
title = genes 102 single transcripts
fontsize = 10
style = flybase
gene rows = 10
prefered_name = gene_name
" > ${ini_file}

pyGenomeTracks --tracks ${ini_file} --region chr6:51,119,395-53,072,988 --outFileName FigS9D_HoxA.pdf
pyGenomeTracks --tracks ${ini_file} --region chr11:95,634,641-97,036,681 --outFileName FigS9D_HoxB.pdf
pyGenomeTracks --tracks ${ini_file} --region chr15:102,417,766-103,467,108 --outFileName FigS9D_HoxC.pdf

# Figure S10F
ini_file="FigS10F.ini"
echo "" > ${ini_file}
for time in {72..120..24}; do
    if [ -e ${pathWithChIP}/wt_${time}h_CDX2_Normalized.bigwig ]; then
        echo "[wt_${time}h_CDX2]
file = ${pathWithChIP}/wt_${time}h_CDX2_Normalized.bigwig
title = ${time}h_CDX2
height = 3
color = #3F4716
min_value = 0
max_value = 20
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig
" >> ${ini_file}
    fi
done

echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 1
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" >> ${ini_file}

pgt --tracks ${ini_file} --region chr2:74,667,374-74,767,842 -o ${ini_file/.ini/.pdf} --width 30

# Figure S11ABCD
fignames=('AB' 'DE')
mutants=('Del(CBS1)' 'Del(CBS1-2)')
for i in "${!fignames[@]}"; do
    mutant=${mutants[$i]}
    ini_file="FigS11${fignames[$i]}.ini"
    echo "[CTCF]
arrowhead_fraction = 0.006
display = collapsed
file = ${pathWithGitHub}/annotations/CTCF.bed
title = CTCF-binding sites
height = 0.3
file_type = bed
labels = false
color = bed_rgb
line_width = 0.05

[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[${mutant}_96h_CTCF]
file = ${pathWithChIP}/${mutant}_96h_CTCF.bigwig
title = CTCF ${mutant} 96h
height = 1.5
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]
height = 1
" > ${ini_file}
    for time in 96 120; do
        echo "[wt_${time}h_H3K27ac]
file = ${pathWithChIP}/wt_${time}h_H3K27ac_bothrep_Normalized_HCN_${time}h.bedgraph
title = wt_${time}h_HCN
height = 1.5
color = #55cf21
min_value = 0
max_value = 380
number_of_bins = 800
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bedgraph

[spacer]
height = 0.1

[${mutant}_${time}h_H3K27ac]
file = ${pathWithChIP}/${mutant}_${time}h_H3K27ac_bothrep_Normalized_HCN_${time}h.bedgraph
title = ${mutant}_${time}h_HCN
height = 1.5
color = #22530d
min_value = 0
max_value = 380
number_of_bins = 800
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bedgraph

[spacer]
height = 0.1

[${mutant}minuswt_${time}h_H3K27ac]
file = ${pathWithChIP}/${mutant}_${time}h_H3K27ac_bothrep_Normalized_HCN_${time}h_minuswt.bedgraph
title = ${mutant}vswt_${time}h_HCN
height = 1.5
color = #22530d
negative_color = #55cf21
number_of_bins = 800
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bedgraph

[spacer]
height = 0.5
" >> ${ini_file}
    done
    echo "[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 12
" >> ${ini_file}

    pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 10
done

# Figure S12A
ini_file="FigS12A.ini"
geno="Del(CBS4)"
echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 1.5
color = #fc8403
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[${geno}_120h_CTCF]
file = ${pathWithChIP}/${geno}_120h_CTCF.bigwig
title = CTCF ${geno} 120h
height = 1.5
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 16
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 15

# Figure S14Atop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep CBS | grep -v "CBS[36]" > vp_CBS_pos.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "CBS-CS38-40" > v4C_CBSCS3840.bed

ini_file="FigS14Atop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_CBS_pos.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_CBSCS3840.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}

# Figure S14Btop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep -v CBS | grep -v "d1-d4" > vp_genes.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "CBS-CS38-40" > v4C_CBSCS3840.bed

ini_file="FigS14Btop.ini"
echo "[wt_168h_CTCF]
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

[vp]
file = vp_genes.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_CBSCS3840.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,645,050-75,200,352 -o ${ini_file/.ini/.pdf}


# Figure S15A
ini_file="FigS15A.ini"
echo "[wt_96h_H3K27me3]
file = ${pathWithChIP}/wt_96h_H3K27me3_bothrep_Normalized_HCN_96h.bedgraph
title = wt_96h_HCN
height = 2
color = #873e23
min_value = 0
max_value = 180
number_of_bins = 800
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bedgraph

[spacer]
height = 0.1

[Del(CBS1-5)_96h_H3K27me3]
file = ${pathWithChIP}/Del(CBS1-5)_96h_H3K27me3_bothrep_Normalized_HCN_96h.bedgraph
title = Del(CBS1-5)_96h_HCN
height = 2
color = #59210c
min_value = 0
max_value = 180
number_of_bins = 800
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bedgraph

[spacer]
height = 0.1

[Del(CBS1-5)minuswt_96h_H3K27me3]
file = ${pathWithChIP}/Del(CBS1-5)_96h_H3K27me3_bothrep_Normalized_HCN_96h_minuswt.bedgraph
title = Del(CBS1-5)vswt_96h_HCN
height = 2
color = #59210c
negative_color = #873e23
min_value = -75
max_value = 75
number_of_bins = 800
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bedgraph

[spacer]

[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 12
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,646,423-74,780,095 -o ${ini_file/.ini/.pdf} --width 15

# Figure S15Btop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep -v "CBS" | grep -v "d1-d4" | grep -v "d8" > vp_genes_S15B.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "Hoxd3-4" > v4C_Hoxd34.bed

ini_file="FigS15Btop.ini"
echo "[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF Del(CBS1-5) 96h
height = 2
color = #ca6a02
min_value = 0
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
file_type = bigwig

[vp]
file = vp_genes_S15B.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_Hoxd34.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,636,423-74,780,095 -o ${ini_file/.ini/.pdf}

# Figure S16A
ini_file="FigS16A.ini"
echo "" > ${ini_file}
for time in {48..168..24}; do
    if [ ${time} = "48" ]; then
        depth=1500000
    else
        depth=1000000
    fi
    if [ -e ${pathWithCHiC}/Del\(CBS1-5\)_${time}h_CHiC_5kb_Balanced_mm10.cool ]; then
        echo "[Del(CBS1-5)_${time}h_CHiC]
file = ${pathWithCHiC}/Del(CBS1-5)_${time}h_CHiC_5kb_Balanced_mm10.cool
title = ${time}h_CHiC_5kb_Balanced
show_masked_bins = false
depth = ${depth}

[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF
height = 2
color = #ca6a02
min_value = 0
number_of_bins = 2000
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig

[spacer]
height = 0.5
" >> ${ini_file}
    fi
    if [ ${time} = "48" ]; then
        echo "[annotations]
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
" >> ${ini_file}
    fi
done
pgt --tracks ${ini_file} --region chr2:73,900,037-75,621,560 -o ${ini_file/.ini/.pdf}

# Figure S16B
genos=('' 'wt' 'Del(CBS1-5)')
for i in 1 2; do
    geno=${genos[$i]}
    for time in 48 96; do
        ini_file=FigS16B${i}-${time}h-1.ini
        echo "" > ${ini_file}
        if [ -e ${pathWithCHiC}/${geno}_${time}h_CHiC_merge_5kb_Balanced_mm10.cool ]; then
            cool_file=${pathWithCHiC}/${geno}_${time}h_CHiC_merge_5kb_Balanced_mm10.cool
        else
            cool_file=${pathWithCHiC}/${geno}_${time}h_CHiC_5kb_Balanced_mm10.cool
        fi
        echo "[${geno}_${time}h_CHiC]
file = ${cool_file}
title = ${geno} ${time}h_CHiC_5kb_Balanced
show_masked_bins = false
depth = 500000

[HoxDCS3840]
file = ${pathWithGitHub}/annotations/HoxDvsCS3840_tight.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dashed
color = white
line_width = 0.5

[line]
file = ${pathWithGitHub}/annotations/FigS16B${i}_${time}.bedpe
file_type = links
links_type = loops
overlay_previous = share-y
line_style = dashed
color = white
line_width = 0.5

[spacer]
height = 0.1
" >> ${ini_file}
        if [ "${geno}" = "Del(CBS1-5)" ]; then
            echo "[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF
height = 2
color = #ca6a02
min_value = 0
max_value = 59.6
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig
" >> ${ini_file}
        else
            echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig
" >> ${ini_file}
        fi
        echo "[vlines]
file = ${pathWithGitHub}/annotations/FigS16B${i}_${time}.bedpe
type = vlines" >> ${ini_file}
        pgt --tracks ${ini_file} --region chr2:74,590,000-74,985,000 -o ${ini_file/.ini/.pdf} --width 10 --dpi 250
    done
done

# Plot the vlines on top of genes
genos=('' 'wt' 'Del(CBS1-5)')
for i in 1 2; do
    geno=${genos[$i]}
    for time in 48 96; do
        ini_file=FigS16B${i}-${time}h-2.ini
        if [ "${geno}" = "Del(CBS1-5)" ]; then
            echo "[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF
height = 2
color = #ca6a02
min_value = 0
max_value = 59.6
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig
" > ${ini_file}
        else
            echo "[wt_168h_CTCF]
file = ${pathWithChIP}/wt_168h_CTCF.bigwig
title = CTCF
height = 2
color = #fc8403
min_value = 0
max_value = 100
number_of_bins = 500
nans_to_zeros = true
summary_method = mean
show_data_range = true
orientation = inverted
file_type = bigwig
" > ${ini_file}
        fi
        echo "[genes]
file = ${pathWithGitHub}/annotations/Hox_single_transcripts.gtf
height = 0.5
title = genes 102 single transcripts
fontsize = 10
style = flybase
prefered_name = gene_name
display = collapsed
labels = false
arrowhead_fraction = 0.017

[labels]
file = Hox_left.bed
color = none
border_color = none
fontstyle = oblique
display = collapsed
height = 1
fontsize = 6

[vlines]
file = ${pathWithGitHub}/annotations/FigS16B${i}_${time}.bedpe
type = vlines
" >> ${ini_file}
        pgt --tracks ${ini_file} --region chr2:74,667,374-74,750,000 -o ${ini_file/.ini/.pdf} --width 10 --dpi 250
    done
done

# Figure S16Ctop
cat ${pathWithGitHub}/annotations/v4C_vp.bed | grep -v "CBS" | grep -v "d1-d4" | grep -v "d8" > vp_genes_S15B.bed
cat ${pathWithGitHub}/annotations/v4C_quantified.bed | grep "5p-sub-TAD1" > v4C_5psubTAD1.bed

ini_file="FigS16Ctop.ini"
echo "[Del(CBS1-5)_96h_CTCF]
file = ${pathWithChIP}/Del(CBS1-5)_96h_CTCF.bigwig
title = CTCF Del(CBS1-5) 96h
height = 2
color = #ca6a02
min_value = 0
number_of_bins = 800
nans_to_zeros = true
summary_method = max
show_data_range = true
file_type = bigwig

[vp]
file = vp_genes_S15B.bed
height = 1
color = blue
display = collapsed
border_color = none
labels = false

[quantifs]
file = v4C_5psubTAD1.bed
overlay_previous = yes
display = collapsed
color = orange
border_color = none
labels = false
" > ${ini_file}
pgt --tracks ${ini_file} --region chr2:74,625,050-75,180,352 -o ${ini_file/.ini/.pdf}
