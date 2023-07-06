#!/usr/bin bash
# This script is for single-end sequence.

# sra2fq.sh
fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/mouse/sequence \
/mnt/xuruizhi/ATAC_brain/mouse/sra/{}/{}.sra

# fastqc
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 \
  -o /mnt/xuruizhi/ATAC_brain/mouse/trim {}.fastq.gz

# fatsqc_again
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}.fastq.gz