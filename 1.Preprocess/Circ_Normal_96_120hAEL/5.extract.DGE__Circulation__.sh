#!/usr/bin/sh
cd $PBS_O_WORKDIR


# Circulation Normal
### Circ 96 h AEL lib1 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib1/error_detected.bam O=lib1/error_detected.expr.txt SUMMARY=lib1/error_detected.summary.txt MIN_NUM_READS_PER_CELL=10000
### Circ 96 h AEL lib2 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib2/error_detected.bam O=lib2/error_detected.expr.txt SUMMARY=lib2/error_detected.summary.txt MIN_NUM_READS_PER_CELL=25000
### Circ 96 h AEL lib3 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib3/error_detected.bam O=lib3/error_detected.expr.txt SUMMARY=lib3/error_detected.summary.txt MIN_NUM_READS_PER_CELL=10000


# Circulation Normal
### Circ 120 h AEL lib1 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib1/error_detected.bam O=lib1/error_detected.expr.txt SUMMARY=lib1/error_detected.summary.txt MIN_NUM_READS_PER_CELL=5000
### Circ 120 h AEL lib2 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib2/error_detected.bam O=lib2/error_detected.expr.txt SUMMARY=lib2/error_detected.summary.txt MIN_NUM_READS_PER_CELL=10000
### Circ 120 h AEL lib3 ###
/home/sangho/programs/src/Drop-seq_tools-1.13/DigitalExpression I=lib3/error_detected.bam O=lib3/error_detected.expr.txt SUMMARY=lib3/error_detected.summary.txt MIN_NUM_READS_PER_CELL=10000

