#!/usr/bin/sh
#cd $PBS_O_WORKDIR

# 1. paired-end fastq files -> bam 
# java -jar <path to Drop-seq tools>/3rdParty/picard/picard.jar FastqToSam F1=<read 1> F2=<read 2> O=<output unaligned bam> SM=<library #>

# 2. Drop-seq alignment tool pipeline
# <path to Drop-seq tools>/Drop-seq_alignment.sh -g <path to STAR index> -r <Drop-seq metadata; fasta> -n <expected cell count; Bumsik expects 2000 cells> -d <path to Drop-seq tools> -s <STAR alinger in full path> -o <path to output directory> -t <path to temporary output directory; create 'tmp' inside the -o path and set it here> -p <unaligned bam>

# 3. Aligend read summarization (per barcode using TAG=XC)
# <path to Drop-seq tools>/BAMTagHistogram I=<aligned bam> O=<output read count> TAG=XC
