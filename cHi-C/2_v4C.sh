#!/bin/bash

pathWithCHiC="/scratch/ldelisle/GEO_Hocine/toGEO/CHiC/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"

mkdir -p /scratch/ldelisle/GEO_Hocine/v4C/
cd /scratch/ldelisle/GEO_Hocine/v4C/

# Get the python script from Amandio et al. 2021:
wget https://raw.githubusercontent.com/lldelisle/scriptsForAmandioEtAl2021/103d950ddc8b92012295a728c58c95b19cecec68/scripts/fromFragFileAndValidPairsToVirtualCaptureC.py -nc
# As well as the digester file:
wget https://raw.githubusercontent.com/lldelisle/scriptsForAmandioEtAl2021/c286f9028d20277af0d6f4ab45539e7c066ce21d/cHi-C/Digester_File_chr2.txt -nc

pathWithVP=${pathWithGitHub}/annotations/v4C_vp.bed

while read l; do
    sample=$(echo $l | awk '{print $1}')
    vp=$(echo $l | awk '{print $2}')
    output=${sample}_${vp}.bedgraph
    if [ -e $output ]; then
        echo "$output already exists"
    else
        if [ ! -e ${sample}.validPair ]; then
            if [ -e ${pathWithCHiC}/${sample}_validPairs_mm10.txt.gz ]; then
                zcat ${pathWithCHiC}/${sample}_validPairs_mm10.txt.gz > ${sample}.validPair
            else
                geno=$(echo $sample | awk -F "_" '{print $1}')
                zcat ${pathWithCHiC}/${sample}_validPairs_${geno}.txt.gz > ${sample}.validPair
            fi
        fi
        coo=$(cat $pathWithVP | awk -v vp=$vp '$4==vp{print $1":"($2+$3)/2}')
        vpReg=$(cat $pathWithVP | awk -v vp=$vp '$4==vp{print ($3-$2)/1000}')
        python fromFragFileAndValidPairsToVirtualCaptureC.py \
            --validPair ${sample}.validPair \
            --colChr1 3 --colChr2 7 --colFrag1 5 --colFrag2 9 --lineToSkipInValidPair 0 \
            --fragmentFile Digester_File_chr2.txt \
            --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 \
            --lineToSkipInFragmentFile 0 --viewpointRegionInKb ${vpReg} \
            --viewpointCoo ${coo} --output ${output}
    fi
done < ${pathWithGitHub}/cHi-C/v4C_table.txt

# Plot FigS6 S7 S9 S10
Rscript ${pathWithGitHub}/scripts/plot_v4C_quantif.R $PWD/ ${pathWithGitHub}
