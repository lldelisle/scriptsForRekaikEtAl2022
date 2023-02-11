# Using homer version 4.10
pathWithChIP="/scratch/ldelisle/GEO_Hocine/toGEO/ChIP/"
gunzip -c ${pathWithChIP}/wt_72h_CDX2.narrowPeak.gz > wt_72h_CDX2.narrowPeak
findMotifsGenome.pl wt_72h_CDX2.narrowPeak mm10 output_CDX2 -size given -len 8,10,12
