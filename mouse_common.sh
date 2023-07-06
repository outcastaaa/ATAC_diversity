#!/usr/bin bash
# This script is for both pari-end and single-end sequence.
# the process of part of post-alignment.

# rmdup
parallel -k -j 96 'picard MarkDuplicates -I ./{}.sort.bam -O ../rmdup/{}.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE ../rmdup/{}.log'
# index
samtools index -@ 96 ../rmdup/{}.rmdup.bam
samtools flagstat -@ 96 ../rmdup/{}.rmdup.bam > ../rmdup/${sample}.rmdup.stat

# rm chrM etal
samtools view -h -f 2 -q 30 ../rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 96 -O bam  -o ../filter/{}.filter.bam
samtools index -@ 96 ../filter/{}.filter.bam
samtools flagstat -@ 96 ../filter/{}.filter.bam > ../filter/{}.filter.stat
