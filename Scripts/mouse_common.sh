#!/usr/bin bash
# This script is for both pari-end and single-end sequence.
# the process of part of post-alignment.

# rmdup
# parallel -k -j 48 'picard MarkDuplicates -I ./{}.sort.bam -O ../rmdup/{}.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE ../rmdup/{}.log'
# index
# samtools index -@ 48 ./{}.rmdup.bam
# samtools flagstat -@ 48 ./{}.rmdup.bam > ./{}.rmdup.stat

# # rm chrM etal
# samtools view -h -F 1804 -q 30 -f 2 ./{}.rmdup.bam | grep -v  chrM | samtools sort -@ 48 -O bam  -o ../filter/{}.filter.bam
# samtools index -@ 48 ../filter/{}.filter.bam
# samtools flagstat -@ 48 ../filter/{}.filter.bam > ../filter/{}.filter.stat


# Blacklist filtering
# cd /mnt/xuruizhi/ATAC_brain/mouse/filter
# echo {} 
#   echo "{}.filter.bam"
#   bedtools intersect -wa -a {}.filter.bam -b ../blklist/mm10.blacklist.bed | \
#   wc -l  > ../blklist/{}.intersect.list

#   bedtools intersect -v -a {}.filter.bam -b ../blklist/mm10.blacklist.bed > ../final/{}.final.bam
#   samtools index ../final/{}.final.bam
#   samtools flagstat ../final/{}.final.bam > ../final/{}.final.stat


# call peaks
# cd /mnt/xuruizhi/ATAC_brain/mouse/final
# bam2bed
# bedtools bamtobed -i {}.final.bam > ../bed/{}.bed
# # shift
cat ../bed/{}.bed | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/{}.Tn5.bed

# # call peaks
macs2 callpeak  -g mm --shift -75 --extsize 150 --nomodel --nolambda --keep-dup all -n {} -t ../Tn5_shift/{}.Tn5.bed --outdir ../peaks/