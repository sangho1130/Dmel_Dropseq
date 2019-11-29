#!/usr/bin/sh
cd $PBS_O_WORKDIR

### 96 h AEL lib1 ###
java -jar /home/sangho/programs/src/Drop-seq_tools-1.13/3rdParty/picard/picard.jar FastqToSam F1=lib1/DropSeq_Normal_LG_96hAEL_lib1_R1.fastq F2=lib1/DropSeq_Normal_LG_96hAEL_lib1_R2.fastq O=lib1/DropSeq_Normal_LG_96hAEL_lib1.bam SM=lib1
/home/sangho/programs/src/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g STARindex-withLincrnas/ -r metafiles/Drosophila_melanogaster.BDGP6.withLincrnas.fasta -n 2000 -d /home/sangho/programs/src/Drop-seq_tools-1.13/ \
	-s /home/sangho/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR -o lib1/ -t lib1/tmp/ -p lib1/DropSeq_Normal_LG_96hAEL_lib1.bam
/home/sangho/programs/src/Drop-seq_tools-1.13/BAMTagHistogram I=lib1/error_detected.bam O=lib1/lib1_cell_readCounts.txt TAG=XC

### 96 h AEL lib2 ###
java -jar /home/sangho/programs/src/Drop-seq_tools-1.13/3rdParty/picard/picard.jar FastqToSam F1=lib2/DropSeq_Normal_LG_96hAEL_lib2_R1.fastq F2=lib2/DropSeq_Normal_LG_96hAEL_lib2_R2.fastq O=lib2/DropSeq_Normal_LG_96hAEL_lib2.bam SM=lib2
/home/sangho/programs/src/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g STARindex-withLincrnas/ -r metafiles/Drosophila_melanogaster.BDGP6.withLincrnas.fasta -n 4000 -d /home/sangho/programs/src/Drop-seq_tools-1.13/ \
        -s /home/sangho/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR -o lib2/ -t lib2/tmp/ -p lib2/DropSeq_Normal_LG_96hAEL_lib2.bam
/home/sangho/programs/src/Drop-seq_tools-1.13/BAMTagHistogram I=lib2/error_detected.bam O=lib2/lib2_cell_readCounts.txt TAG=XC

### 96 h AEL lib3 ###
java -jar /home/sangho/programs/src/Drop-seq_tools-1.13/3rdParty/picard/picard.jar FastqToSam F1=lib3/DropSeq_Normal_LG_96hAEL_lib3_R1.fastq F2=lib3/DropSeq_Normal_LG_96hAEL_lib3_R2.fastq O=lib3/DropSeq_Normal_LG_96hAEL_lib3.bam SM=lib3
/home/sangho/programs/src/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g STARindex-withLincrnas/ -r metafiles/Drosophila_melanogaster.BDGP6.withLincrnas.fasta -n 3000 -d /home/sangho/programs/src/Drop-seq_tools-1.13/ \
        -s /home/sangho/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR -o lib3/ -t lib3/tmp/ -p lib3/DropSeq_Normal_LG_96hAEL_lib3.bam
/home/sangho/programs/src/Drop-seq_tools-1.13/BAMTagHistogram I=lib3/error_detected.bam O=lib3/lib3_cell_readCounts.txt TAG=XC

### 96 h AEL lib4 ###
java -jar /home/sangho/programs/src/Drop-seq_tools-1.13/3rdParty/picard/picard.jar FastqToSam F1=lib4/DropSeq_Normal_LG_96hAEL_lib4_R1.fastq F2=lib4/DropSeq_Normal_LG_96hAEL_lib4_R2.fastq O=lib4/DropSeq_Normal_LG_96hAEL_lib4.bam SM=lib4
/home/sangho/programs/src/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g STARindex-withLincrnas/ -r metafiles/Drosophila_melanogaster.BDGP6.withLincrnas.fasta -n 4500 -d /home/sangho/programs/src/Drop-seq_tools-1.13/ \
        -s /home/sangho/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR -o lib4/ -t lib4/tmp/ -p lib4/DropSeq_Normal_LG_96hAEL_lib4.bam
/home/sangho/programs/src/Drop-seq_tools-1.13/BAMTagHistogram I=lib4/error_detected.bam O=lib4/lib4_cell_readCounts.txt TAG=XC

### 96 h AEL lib5 ###
java -jar /home/sangho/programs/src/Drop-seq_tools-1.13/3rdParty/picard/picard.jar FastqToSam F1=lib5/DropSeq_Normal_LG_96hAEL_lib5_R1.fastq F2=lib5/DropSeq_Normal_LG_96hAEL_lib5_R2.fastq O=lib5/DropSeq_Normal_LG_96hAEL_lib5.bam SM=lib4
/home/sangho/programs/src/Drop-seq_tools-1.13/Drop-seq_alignment.sh -g STARindex-withLincrnas/ -r metafiles/Drosophila_melanogaster.BDGP6.withLincrnas.fasta -n 4000 -d /home/sangho/programs/src/Drop-seq_tools-1.13/ \
        -s /home/sangho/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR -o lib5/ -t lib5/tmp/ -p lib5/DropSeq_Normal_LG_96hAEL_lib5.bam
/home/sangho/programs/src/Drop-seq_tools-1.13/BAMTagHistogram I=lib5/error_detected.bam O=lib5/lib5_cell_readCounts.txt TAG=XC

