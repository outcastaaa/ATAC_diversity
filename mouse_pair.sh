#!/usr/bin bash
# This script is for pari-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/mouse/sequence /mnt/xuruizhi/ATAC_brain/mouse/sra/{}/{}.sra

# # fastqc
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_1.fastq.gz
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_2.fastq.gz

# # trim
# trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired \
# -o /mnt/xuruizhi/ATAC_brain/mouse/trim  {}_1.fastq.gz  {}_2.fastq.gz

# # fatsqc_again
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_1.fastq.gz
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_2.fastq.gz

# align
bowtie2  -p 48 -x /scratch/wangq/xrz/ATAC_brain/mouse/genome/mm10 \
--very-sensitive -X 2000 -1 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_1_trimmed.fq.gz \
-2 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_2_trimmed.fq.gz \
-S /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam

# sort_transfertobam_index 
samtools sort -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam \
> /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools index -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools flagstat  -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam \
> /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.raw.stat