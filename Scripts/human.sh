#!/usr/bin bash
# This script is for pari-end sequence.


# local
# sra2fq.sh
fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/human/sequence /mnt/xuruizhi/ATAC_brain/human/sra/{}/{}.sra

# fastqc
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/human/fastqc /mnt/xuruizhi/ATAC_brain/human/sequence/{}_1.fastq.gz
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/human/fastqc /mnt/xuruizhi/ATAC_brain/human/sequence/{}_2.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o /mnt/xuruizhi/ATAC_brain/human/trim  {}_1.fastq.gz  {}_2.fastq.gz

# fatsqc_again
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/human/fastqc_again /mnt/xuruizhi/ATAC_brain/human/trim/{}_1_val_1.fq.gz
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/human/fastqc_again /mnt/xuruizhi/ATAC_brain/human/trim/{}_2_val_2.fq.gz

# hpcc
# align比对尽量不要拆多节点去跑，有可能造成数据读取错误
bowtie2  -p 20 -x /scratch/wangq/xrz/ATAC_brain/human/genome/GRCh38_noalt_as/GRCh38_noalt_as  --very-sensitive -X 2000 -1 /scratch/wangq/xrz/ATAC_brain/human/trim/{}_1_val_1.fq.gz -2 /scratch/wangq/xrz/ATAC_brain/human/trim/{}_2_val_2.fq.gz -S /scratch/wangq/xrz/ATAC_brain/human/align/{}.sam

# sort_transfertobam_index 
samtools sort -@ 20 /scratch/wangq/xrz/ATAC_brain/human/align/{}.sam > /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam
samtools flagstat  -@ 8 /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam > /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.raw.stat

# rmdup
# cd /scratch/wangq/xrz/ATAC_brain/human/sort_bam
parallel -k -j 20 'picard MarkDuplicates -I /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam -O /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.log'
# index
# cd /scratch/wangq/xrz/ATAC_brain/human/rmdup
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam
samtools flagstat -@ 8 /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam > /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.stat

# rm chrM etal
samtools view -h -f 2 -F 1804 -q 30  /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 20 -O bam  -o /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam
samtools flagstat -@ 8 /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam > /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.stat


# local
# rm blklist
bedtools intersect -wa -a {}.filter.bam -b ../blklist/hg38.blacklist.bed | \
wc -l  > ../blklist/{}.intersect.list
bedtools intersect -v -a {}.filter.bam -b ../blklist/hg38.blacklist.bed > ../final/{}.final.bam
samtools index -@ 6 ../final/{}.final.bam
samtools flagstat -@ 6 ../final/{}.final.bam > ../final/{}.final.stat

# shift bed
# cd /mnt/xuruizhi/ATAC_brain/human/final
bedtools bamtobed -i {}.final.bam > ../bed/{}.bed
cat ../bed/{}.bed | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/{}.Tn5.bed

# cd /mnt/xuruizhi/ATAC_brain/human/final
bedtools bamtobed -i {}.final.bam | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/{}.Tn5.bed

macs2 callpeak  -g hs --shift -75 --extsize 150 --nomodel --nolambda --keep-dup all -n {} -t ../Tn5_shift/{}.Tn5.bed --outdir ../peaks/