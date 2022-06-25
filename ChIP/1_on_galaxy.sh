# Here are the command lines launched by galaxy when using
# the workflow Galaxy-Workflow-inported__ChIPseqOrChipMentationPairEnd_cutadaptBoth_withMAPQandCPfilter_newMacs2_v20190602.ga

# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.1
# command_version:1.16

cutadapt -j ${GALAXY_SLOTS:-1} -a 'Please use: For R1: - For Nextera: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC - For TrueSeq: GATCGGAAGAGCACACGTCTGAACTCCAGTCAC or AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC '='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'Please use: For R2: - For Nextera: CTGTCTCTTATACACATCTGACGCTGCCGACGA - For TruSeq: GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT or AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' --output='out1.gz' --paired-output='out2.gz' 'sample_R1.fq.gz' 'sample_R2.fq.gz' > report.txt

# toolshed.g2.bx.psu.edu/repos/devteam/bowtie2/bowtie2/2.3.4.1
# command_version:/data/galaxy/galaxy/var/dependencies/_conda/envs/mulled-v1-65d5efe4f1b69ab7166d1a5a5616adebe902133ea3e4c189d87d7de2e21ddc17/bin/bowtie2-align-s version 2.3.4.1
# 64-bit
# Built on default-df05fd51-3d07-4109-abba-6883676f3ae8
# Mon Jun 25 23:07:08 UTC 2018
# Compiler: gcc version 4.8.2 20140120 (Red Hat 4.8.2-15) (GCC) 
# Options: -O3 -m64 -msse2 -funroll-loops -g3  -DBOOST_MATH_DISABLE_FLOAT128 -m64 -fPIC -std=c++98 -DPOPCNT_CAPABILITY -DWITH_TBB -DNO_SPINLOCK -DWITH_QUEUELOCK=1
# Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

set -o | grep -q pipefail 
 set -o pipefail;   ln -f -s 'cutadapt of sample_R1.fastqsanger.gz' input_f.fastq.gz 
  ln -f -s 'cutadapt of sample_R2.fastqsanger.gz' input_r.fastq.gz 
   bowtie2  -p ${GALAXY_SLOTS:-4}  -x '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/bowtie2_index/mm10_UCSC/mm10_UCSC'   -1 'input_f.fastq.gz' -2 'input_r.fastq.gz'                      2> 'mapping stats of cutadapt of sample .txt'  | samtools sort -@${GALAXY_SLOTS:-2} -O bam -o 'mapping of cutadapt of sample.bam'

# toolshed.g2.bx.psu.edu/repos/devteam/samtool_filter2/samtool_filter2/1.1.2
# command_version:
ln -s 'mapping of cutadapt of sample.bam' input.bam 
 ln -s '/data/galaxy/data/_metadata_files/012/metadata_12490.dat' input.bai 
 samtools view -o 'filtered bam of mapping of cutadapt of sample.bam' -h   -b  -q 30 -f 0x2 input.bam 2>&1

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.1.1.20160309.3
# command_version:macs2 2.1.1.20160309

export PYTHON_EGG_CACHE=`pwd` 
 (macs2 callpeak --name 'MACS2' -t 'filtered bam of mapping of cutadapt of sample.bam' --format BAMPE --gsize '1870000000' --call-summits --bdg 2>&1 > macs2_stderr) 
 cp MACS2_peaks.xls 'peaks xls of filtered bam of mapping of cutadapt of sample.tabular' 
 ( count=`ls -1 MACS2* 2>/dev/null | wc -l`; if [ $count != 0 ]; then mkdir '/data/galaxy/galaxy/jobs/090/90933/working/dataset_252695_files' 
 cp -r MACS2* '/data/galaxy/galaxy/jobs/090/90933/working/dataset_252695_files' 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/macs2/01cded2297b7/macs2/dir2html.py' '/data/galaxy/galaxy/jobs/090/90933/working/dataset_252695_files' macs2_stderr > 'MACS2 callpeak on data 46 (html report).html'; fi; ) 
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# This awk command remove lines with no coverage and add a header
# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.0
# command_version:GNU Awk 4.1.3, API: 1.1

awk --sandbox -v FS='	' -v OFS='	' --re-interval -f "/data/galaxy/galaxy/jobs/090/90935/configs/tmp31n6vy4r" "macs2treatment coverage of filtered bam of mapping of cutadapt of sample.bedgraph" > "macs2treatment coverage of filtered bam of mapping of cutadapt of sample with header.bedgraph"

# wig_to_bigWig
# command_version:
grep -v "^track" 'macs2treatment coverage of filtered bam of mapping of cutadapt of sample with header.bedgraph' | wigToBigWig stdin '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len' 'bigWig macs2treatment coverage of filtered bam of mapping of cutadapt of sample with header.bigwig' -clip 2>&1 || echo "Error running wigToBigWig." >&2
