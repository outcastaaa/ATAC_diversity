#!/usr/bin bash

# 转格式，sra2fq
fastq-dump --gzip --split-3 -O /mnt/xuruizhi/RNA_brain/mouse/sequence /mnt/xuruizhi/RNA_brain/mouse/sra/{}/{}.sra