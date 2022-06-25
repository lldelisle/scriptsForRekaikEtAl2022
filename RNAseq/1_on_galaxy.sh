# Here are the command lines launched by galaxy when using the workflow Galaxy-Workflow-RNAseq_PE_Cufflinks.ga

# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16

cutadapt -j ${GALAXY_SLOTS:-1} -a 'TruSeqR1'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'TruSeq'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' --output='out1.fq.gz' --paired-output='out2.fq.gz' 'sample_R1.fq.gz' 'sample_R2.fq.gz' > report.txt

# toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.7a
# command_version:
STAR --runThreadN ${GALAXY_SLOTS:-4} --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 99 --sjdbGTFfile 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz.gtf' --readFilesIn 'cutadapt of sample_R1.fastq.gz' 'cutadapt of sample_R2.fastq.gz' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate '' --quantMode GeneCounts --outSAMattributes NH HI AS nM --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outBAMsortingThreadN ${GALAXY_SLOTS:-4} --limitBAMsortRAM $((${GALAXY_MEMORY_MB:-0}*1000000)) 
 samtools view -b -o 'RNA STAR on cutadapt of sample.bam' Aligned.sortedByCoord.out.bam

# toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.4.1
# command_version:
cp '/data/galaxy/galaxy/jobs/091/91052/configs/tmpzpq2o7bp' 'Filter on data 81: JSON filter rules.txt' 
 ln -s 'RNA STAR on cutadapt of sample.bam' localbam.bam 
 ln -s '/data/galaxy/data/_metadata_files/012/metadata_12499.dat' localbam.bam.bai 
 cat '/data/galaxy/galaxy/jobs/091/91052/configs/tmpzpq2o7bp' 
 bamtools filter -script '/data/galaxy/galaxy/jobs/091/91052/configs/tmpzpq2o7bp' -in localbam.bam -out 'uniquely mapped of RNA STAR on cutadapt of sample.bam'

# toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/cufflinks/2.2.1.3
# command_version:cufflinks v2.2.1

cufflinks -q --no-update-check 'uniquely mapped of RNA STAR on cutadapt of sample.bam' --num-threads "${GALAXY_SLOTS:-4}" -G 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz.gtf' 2> stderr 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/d080005cffe1/cufflinks/mass.py' stderr 'None' "transcripts.gtf"
