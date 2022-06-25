#!/bin/bash

pathWithChIPoutputs="/scratch/ldelisle/GEO_Hocine/ChIP/"
pathWithGitHub="/home/ldelisle/softwares/scriptsForRekaikEtAl2022/"
pathWithTable="${pathWithGitHub}/ChIP/chip_norm_groups.txt"
pathWithPythonScript="${pathWithGitHub}/scripts/BigWigScaling.py"

groups=$(cat $pathWithTable | awk -F "\t" 'NR>1 && $2!=""{print $2}' | sort | uniq)
for group in $groups; do
    controlSample=$(cat $pathWithTable | awk -F "\t" -v g=$group 'NR>1 && $2==g && $3=="X"{print $1}')
    otherSamples=$(cat $pathWithTable | awk -F "\t" -v g=$group 'NR>1 && $2==g && $3!="X"{print $1}')
    otherBW=""
    otherNP=""
    for os in $otherSamples; do
        otherBW="$otherBW ${pathWithChIPoutputs}/${os}.bigwig"
        otherNP="$otherNP ${pathWithChIPoutputs}/${os}.narrowPeak"
    done
    python $pathWithPythonScript \
        --referenceBW ${pathWithChIPoutputs}/${controlSample}.bigwig \
        --referenceNP ${pathWithChIPoutputs}/${controlSample}.narrowPeak \
        --otherBW ${otherBW} --otherNP ${otherNP} --prefixReports group${group}_
done
