#!/bin/sh

##input and output##

insam=$1
outbam=`basename $insam .sam`.bam
sortoutbam=`basename $insam .sam`.sort.bam
filtbam=`basename $insam .sam`.unique.bam
sortbam=`basename $insam .sam`.unique_sort.bam
finalbam=`basename $insam .sam`.final.bam
finalsortbam=`basename $insam .sam`.final_sort_by_pos.bam
bwfile=`basename $finalsortbam .bowtiealignment.final_sort_by_pos.bam`.bw
picardlog=`basename $insam .sam`.rmdup.log


#transfor sam to bam##
samtools view -b -S $insam -o $outbam


##sort raw bam##
samtools sort -@4 -T ./tmp -o $sortoutbam $outbam

##index raw-alignment  bam##
samtools index $sortoutbam


##statistics flagstat##


##Alignment filter##
samtools view -F 4 -q 1 -u $outbam -o $filtbam

picard SortSam -SO coordinate --INPUT $filtbam  --OUTPUT $sortbam --TMP_DIR `pwd`

##mark and remove duplicates##

picard MarkDuplicates --INPUT $sortbam   --OUTPUT $finalbam --METRICS_FILE $picardlog --CREATE_INDEX false --REMOVE_DUPLICATES true --TMP_DIR `pwd`


##samtools sort by position###
samtools sort -@4 -T ./tmp -o $finalsortbam $finalbam

##samtools index finalsortbam##
samtools index $finalsortbam

##transfer bam to bigWig##
bamCoverage --normalizeUsing RPKM -b $finalsortbam -o $bwfile
