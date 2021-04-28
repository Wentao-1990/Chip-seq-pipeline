# Chip-seq-pipeline
Simply the pipeline for Chipseq data
Pipeline for Chip-seq data analysis


Prerequirements:
1. cutadapt 3.4
2. trim_galore 0.6.6
3. bowtie2 2.4.2
4. macs2 2.2.7.1
4. R 4.0.3
5. ChIPseeker 1.26.2
6. deeptools 3.5.1

##pipeline for chip-seq analysis###
Overview of pipeline:
1. quality control(QC)
filter low-quality reads and remove adapters 
2. index genome
 bwa index -p genome.fasta
3. alignment

4. alignment filter
 	a. alignment qulity >1 && sort according to the coordinate
b. mark and remove pcr-duplication
5. assesment of the alignment 
a. sample correaltaion/ PCA analysis
b. enrichment efficity of antibody
6. Peakcalling
7.Peak annotation

Command line
1. QC and alingment are concentrated in one script:
perl chipseq_processing.pl -read1 20060FL-04-02-04_S51_L004_R1_001.fastq.gz,20060FL-04-02-05_S52_L004_R1_001.fastq.gz,20060FL-04-02-06_S53_L004_R1_001.fastq.gz -read2 20060FL-04-02-04_S51_L004_R2_001.fastq.gz,20060FL-04-02-05_S52_L004_R2_001.fastq.gz,20060FL-04-02-06_S53_L004_R2_001.fastq.gz -sample o3input_1,o3input_2,o3input_3 -rawdata 1 -genome ../../../../genomes/hg19/hg19

-read1: sample1_R1, if multiple samples were provides, 
-read2: sample1_R2
-genome: indexed genome 
-rawdata: 1 indocates raw data; 0 indicates clean data, which will do the alignment directly.
2. Alignment results filter
sh sam

##peak calling##
nohup macs2 callpeak -t ../../molm13_hoxa9/o3hoxa9_3new/o3hoxa9_3new.bowtiealignment.fianl_sort_by_pos.bam -c ../../molm13_input/o3hoxa9_3/o3hoxa9_3.bowtiealignment.fianl_sort_by_pos.bam -B --outdir macs2 -g hs -f BAM -n moml13_rep3 2>runmacserr.log&