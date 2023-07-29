#!/usr/bin bash
# This script is for pari-end sequence.

# 本地
# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/mouse/sequence /mnt/xuruizhi/ATAC_brain/mouse/sra/{}/{}.sra

# # fastqc
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_1.fastq.gz
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_2.fastq.gz

# # trim
# trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired \
# -o /mnt/xuruizhi/ATAC_brain/mouse/trim  {}_1.fastq.gz  {}_2.fastq.gz

# # fatsqc_again
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_1_val_1.fq.gz
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_2_val_2.fq.gz

# 超算
# align
bowtie2  -p 96 -x /scratch/wangq/xrz/ATAC_brain/mouse/genome/mm10 --very-sensitive -X 2000 -1 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_1_val_1.fq.gz -2 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_2_val_2.fq.gz -S /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam


# sort_transfertobam_index 
samtools sort -@ 96 /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam > /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools index -@ 96 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools flagstat  -@ 96 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam > /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.raw.stat

# rmdup
# cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam
parallel -k -j 48 'picard MarkDuplicates -I /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam -O /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.log'
# index
# cd /scratch/wangq/xrz/ATAC_brain/mouse/rmdup
samtools index -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.rmdup.bam
samtools flagstat -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.rmdup.bam > /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.rmdup.stat

# rm chrM etal
samtools view -h -f 2 -F 1804 -q 30  /scratch/wangq/xrz/ATAC_brain/mouse/rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 48 -O bam  -o /scratch/wangq/xrz/ATAC_brain/mouse/filter/{}.filter.bam
samtools index -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/filter/{}.filter.bam
samtools flagstat -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/filter/{}.filter.bam > /scratch/wangq/xrz/ATAC_brain/mouse/filter/{}.filter.stat


# 本地
# rm blklist
  bedtools intersect -wa -a {}.filter.bam -b ../blklist/mm10.blacklist.bed | \
  wc -l  > ../blklist/{}.intersect.list

parallel -k -j 6 "
  echo {} 
  echo "{}.filter.bam"
  bedtools intersect -wa -a {}.filter.bam -b ../blklist/mm10.blacklist.bed | \
  wc -l  > ../blklist/{}.intersect.list

  bedtools intersect -v -a {}.filter.bam -b ../blklist/mm10.blacklist.bed > ../final/{}.final.bam
  samtools index ../final/{}.final.bam
  samtools flagstat ../final/{}.final.bam > ../final/{}.final.stat
"
