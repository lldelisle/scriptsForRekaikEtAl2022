# Here are the command lines launched by galaxy when using
# the workflow Galaxy-Workflow-imported__Hi-C_fromFastqToValidPairs_mm10_DpnII.ga followed by Galaxy-Workflow-Hi-C_fromValidPairAndChrSizeTo2kbAnd5kbAnd10kbMatrixInHoxDRegion.ga

# genome can be mm10 or Del(d1-d4) or Del(sub-TAD1)

# toolshed.g2.bx.psu.edu/repos/bgruening/hicup_digester/hicup_digester/0.6.1.0
# command_version:HiCUP v0.6.1

hicup_digester --re1 '^GATC' --genome 'genome' genome.fasta

# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16

cutadapt -j ${GALAXY_SLOTS:-1} -a 'TruSeqR1'='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'TruSeqR2'='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' --output='out1.fq.gz' --paired-output='out2.fq.gz' 'sample_R1.fq.gz' 'sample_R2.fq.gz' > report.txt

# toolshed.g2.bx.psu.edu/repos/bgruening/hicup_hicup/hicup_hicup/0.6.1.0
# command_version:HiCUP v0.6.1
 BOWTIE_PATH_BASH="$(which bowtie2)" 
 bowtie2-build "genome.fasta" genome 
 ln -s "genome.fasta" genome.fa 
 hicup_digester --re1 '^GATC' --genome 'genome' genome.fasta 
 mv *Digest_* digest_file.txt 
 hicup --zip --threads ${GALAXY_SLOTS:-1} --digest digest_file.txt --index 'genome' --bowtie2 $BOWTIE_PATH_BASH --keep cutadapt_R1.fastq.gz cutadapt_R2.fastq.gz  
 ls -al

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/fromHicupToJuicebox/0.0.2
# This tool can be downloaded from https://testtoolshed.g2.bx.psu.edu/repository/download?repository_id=be5040251cd4afb7&changeset_revision=44365a4feb3b&file_type=gz
# command_version:
python /data/galaxy/galaxy/var/shed_tools/testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/552eb7782435/from_hicup_to_juicebox/fromHicupToJuicebox.py --fragmentFile digest_file.txt --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --methodForFrag hicup --useMid --output validPairs.txt hicup.bam
 bash /data/galaxy/galaxy/var/shed_tools/testtoolshed.g2.bx.psu.edu/repos/lldelisle/from_hicup_to_juicebox/552eb7782435/from_hicup_to_juicebox/switchAndSort.sh validPairs.txt 'sample_validPairs file with midFrag positions.tabular'

# Filter1
# command_version:
python '/data/galaxy/galaxy/server/tools/stats/filtering.py' 'sample_validPairs file with midFrag positions.tabular' 'both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabular' '/data/galaxy/galaxy/jobs/082/82279/configs/tmp5nfsi_r8' 11 "str,int,str,int,int,int,str,int,int,int,int" 0

# Filter1
# command_version:
python '/data/galaxy/galaxy/server/tools/stats/filtering.py' 'both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabular' 'both pairs in captured region of both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabular' '/data/galaxy/galaxy/jobs/082/82280/configs/tmpkh5xk6ys' 11 "str,int,str,int,int,int,str,int,int,int,int" 0

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_csort_pairix/0.0.1
# command_version:cooler, version 0.7.4
cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o 'sorted and indexed contact list of both pairs in captured region of both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabix' 'both pairs in captured region of both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabular' genome.chrom.sizes.tabular

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_makebins/0.0.1
# command_version:cooler, version 0.7.4

cooler makebins -o bins.5kb_genome.bed genome.chrom.sizes.tabular 5000
cooler makebins -o bins.2kb_genome.bed genome.chrom.sizes.tabular 2000

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_cload_pairix/0.0.1
# command_version:cooler, version 0.7.4

cooler cload tabix --assembly mm10 -c2 7 -p2 8 bins.5kb_genome.bed 'sorted and indexed contact list of both pairs in captured region of both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabix' 'cool_file_with_matrices_of_5kb_Raw.cool'
cooler cload tabix --assembly mm10 -c2 7 -p2 8 bins.2kb_genome.bed 'sorted and indexed contact list of both pairs in captured region of both pairs above MAPQ30 of sample_validPairs file with midFrag positions.tabix' 'cool_file_with_matrices_of_2kb_Raw.cool'

# testtoolshed.g2.bx.psu.edu/repos/lldelisle/cooler/cooler_balance/0.0.1
# command_version:cooler, version 0.7.4

cp 'cool_file_with_matrices_of_5kb_Raw.cool' 'cool_file_with_matrices_of_5kb_Balanced.cool'
# For all genomes except Del(sub-TAD1)
cooler balance --cis-only -f 'cool_file_with_matrices_of_5kb_Balanced.cool'
# For Del(sub-TAD1)
cooler balance --ignore-diags 5 --cis-only -f 'cool_file_with_matrices_of_5kb_Balanced.cool'

# For mm10 genome
cp 'cool_file_with_matrices_of_2kb_Raw.cool' 'cool_file_with_matrices_of_2kb_Balanced.cool'
cooler balance --cis-only -f 'cool_file_with_matrices_of_2kb_Balanced.cool'
