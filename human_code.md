- [0. 安装下载](#0-安装下载)
- [1. 建目录](#1-建目录)
- [2. 下载数据](#2-下载数据)
- [3. 比对前质控](#3-比对前质控)
- [4. 比对](#4-比对)
- [5. post-alignment](#5-post-alignment)
- [6. call peak](#6-call-peaks)
- [7. Quality cheak](#7-quality-check)
- [8. Visualization](#8-visualization)
- [9. 相同组织rep间consensus peak](#9-寻找rep间consensus-peak) 
- [10. 每个脑区独有peak](#10-每个脑区独有peak)
- [11. 每个脑区共有peak](#11-每个脑区共有peak)  
- [12. 对应RNA-seq数据](#12-对应rna-seq数据)


# 0. 安装下载
```bash
# /mnt/d/perl/perl_scripts/download_srr.pl
cpan Getopt::Long 
cpan File::Basename
cpan Digest::MD5
cpan IPC::System::Simple
```
# 1. 建目录
```bash
cd /mnt/xuruizhi # 挂载到NAS
mkdir -p /mnt/xuruizhi/ATAC_brain/human
```

# 2. 下载数据
1. 测序数据
```bash
# human 68个+26个
cd /mnt/xuruizhi/ATAC_brain/human
vim HUMAN.list 
SRR21163180
SRR21163181
SRR21163184
SRR21163185
SRR21163186
SRR21163187
SRR21163190
SRR21163191
SRR21163196
SRR21163197
SRR21163203
SRR21163204
SRR21163207
SRR21163208
SRR21163209
SRR21163210
SRR21163213
SRR21163214
SRR21163216
SRR21163217
SRR21163218
SRR21163219
SRR21163220
SRR21163221
SRR21163226
SRR21163227
SRR21163228
SRR21163229
SRR21163232
SRR21163233
SRR21163234
SRR21163235
SRR21163240
SRR21163241
SRR21163249
SRR21163250
SRR21163254
SRR21163255
SRR21163256
SRR21163257
SRR21163267
SRR21163268
SRR21163293
SRR21163294
SRR21163298
SRR21163299
SRR21163304
SRR21163305
SRR21163320
SRR21163321
SRR21163322
SRR21163323
SRR21163335
SRR21163336
SRR21163337
SRR21163338
SRR21163343
SRR21163344
SRR21163347
SRR21163348
SRR21163349
SRR21163350
SRR21163365
SRR21163366
SRR21163367
SRR21163368
SRR21163376
SRR21163377

# PAC+PVC有必要再下载
SRR21163168
SRR21163169
SRR21163174
SRR21163175
SRR21163182
SRR21163183
SRR21163192
SRR21163193
SRR21163198
SRR21163199
SRR21163326
SRR21163327
SRR21163328
SRR21163329
SRR21163362
SRR21163363

# 加上HAB +PSC 
SRR21163263
SRR21163264
SRR21163311
SRR21163312
SRR21163339
SRR21163340
SRR21163357
SRR21163358
SRR21163371
SRR21163372


mkdir -p ./sra
cd ./sra
cp /mnt/d/perl/perl_scripts/download_srr.pl ./
cp /mnt/xuruizhi/ATAC_brain/human/HUMAN.list ./

# 批量下载
cat new_HUMAN.list | parallel -k -j 8 "
  echo {} >> ./download.log
  perl download_srr.pl --output-dir . --srr {} >> ./download.log 2>&1
"
cat download.log | grep " downloaded successfully"
```

# 3. sra2fa+qc+trim+qc
```bash
cd /mnt/xuruizhi/ATAC_brain/human
mkdir -p sequence
mkdir -p fastqc
mkdir -p trim
mkdir -p fastqc_again

cd sequence 
# 因为磁盘空间不够，分多部分处理
vim 1.list
# PSM
SRR21163180
SRR21163181
SRR21163186
SRR21163187
SRR21163293
SRR21163294
# VLPFC
SRR21163184
SRR21163185
SRR21163203
SRR21163204
SRR21163207
SRR21163208
SRR21163320
SRR21163321
SRR21163337
SRR21163338
# CERE
SRR21163190
SRR21163191
SRR21163209
SRR21163210
SRR21163216
SRR21163217
SRR21163365
SRR21163366


vim 4.list
# OFC
SRR21163196
SRR21163197
SRR21163218
SRR21163219
SRR21163220
SRR21163221
SRR21163298
SRR21163299
SRR21163347
SRR21163348
# NACC
SRR21163213
SRR21163214
SRR21163240
SRR21163241
SRR21163343
SRR21163344
# CN
SRR21163226
SRR21163227
SRR21163249
SRR21163250
SRR21163256
SRR21163257
SRR21163367
SRR21163368


vim 5.list
# AMY
SRR21163228
SRR21163229
SRR21163234
SRR21163235
SRR21163335
SRR21163336
# DLPFC
SRR21163232
SRR21163233
SRR21163254
SRR21163255
SRR21163322
SRR21163323
SRR21163349
SRR21163350
SRR21163376
SRR21163377
# HIPP
SRR21163267
SRR21163268
SRR21163304
SRR21163305



vim 2.list
SRR21163168
SRR21163169
SRR21163174
SRR21163175
SRR21163182
SRR21163183
SRR21163192
SRR21163193
SRR21163198
SRR21163199
SRR21163287
SRR21163288
SRR21163326
SRR21163327
SRR21163328
SRR21163329
SRR21163362
SRR21163363


vim 3.list
SRR21163263
SRR21163264
SRR21163311
SRR21163312
SRR21163339
SRR21163340
SRR21163357
SRR21163358
SRR21163371
SRR21163372



vim human.sh
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




cat 3.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ${id}_qc_trim.sh
done
cat 3.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_qc_trim.sh >> ./sra2qc_trim_fastqc.log 2>&1"
# 其他list同理
cat 1.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ${id}_qc_trim.sh
done
cat 1.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_qc_trim.sh >> ./sra2qc_trim_fastqc.log 2>&1"

cd ../fastqc
multiqc .
cd  ../fastqc_again
multiqc .
```
# 4. align_filter

1. genome: 后续用bowtie2比对，在bowtie2官网下载已经建立好的基因组索引文件。   
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/genome
mkdir -p /mnt/d/ATAC_brain/human/genome
cd /mnt/d/ATAC_brain/human/genome
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip 
unzip GRCh38_noalt_as.zip

cp -r ./* /mnt/xuruizhi/ATAC_brain/human/genome/
```
2. 批量处理
```bash
# 转入超算
rsync -avP /mnt/xuruizhi/ATAC_brain/human/sequence \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -avP /mnt/xuruizhi/ATAC_brain/human/trim \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -avP /mnt/xuruizhi/ATAC_brain/human/sra \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -avP /mnt/xuruizhi/ATAC_brain/human/fastqc \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -avP /mnt/xuruizhi/ATAC_brain/human/fastqc_again \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -avP /mnt/xuruizhi/ATAC_brain/human/genome \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/


cd /scratch/wangq/xrz/ATAC_brain/human
mkdir -p align 
mkdir -p sort_bam
mkdir -p filter
mkdir -p rmdup
mkdir -p blklist

cd align
cp ../sequence/*.list .

```
* align
```bash
# align
vim human.sh
#!/usr/bin bash
# This script is for pari-end sequence.

# hpcc
# align
bowtie2  -p 20 -x /scratch/wangq/xrz/ATAC_brain/human/genome/GRCh38_noalt_as/GRCh38_noalt_as --very-sensitive -X 2000 -1 /scratch/wangq/xrz/ATAC_brain/human/trim/{}_1_val_1.fq.gz -2 /scratch/wangq/xrz/ATAC_brain/human/trim/{}_2_val_2.fq.gz -S /scratch/wangq/xrz/ATAC_brain/human/align/{}.sam

# sort_transfertobam_index 
samtools sort -@ 20 /scratch/wangq/xrz/ATAC_brain/human/align/{}.sam > /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam
samtools flagstat  -@ 8 /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam > /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.raw.stat


cat 3.list  | while read id
do 
  sed "s/{}/${id}/g" human.sh > ${id}_align.sh
  bsub -q mpi -n 24  "
  bash  ${id}_align.sh >> ./align.log 2>&1"
done
```
* rmdup
```bash
# rmdup
parallel -k -j 20 --no-run-if-empty --linebuffer "\
picard MarkDuplicates -I /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam -O /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.log" :::: 3.list
```
* rm chrM etal
```bash
vim human.sh
#!/usr/bin bash
# This script is for pari-end sequence.

# hpcc
# index
# cd /scratch/wangq/xrz/ATAC_brain/human/rmdup
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam
samtools flagstat -@ 8 /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam > /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.stat

# rm chrM etal
samtools view -h -f 2 -F 1804 -q 30  /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 20 -O bam  -o /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam
samtools index -@ 8 /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam
samtools flagstat -@ 8 /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.bam > /scratch/wangq/xrz/ATAC_brain/human/filter/{}.filter.stat


cd rmdup
cp ../sequence/3.list .
cat 3.list  | while read id
do 
  sed "s/{}/${id}/g" human.sh > ${id}_filter.sh
  bsub -q mpi -n 24  "
  bash  ${id}_filter.sh >> ./align.log 2>&1"
done
```
* rm blklist

```bash
# 转入本地，使用beyond compare需要使用md5检查后再操作

# 下载hg38的blklist
mkdir -p /mnt/xuruizhi/ATAC_brain/human/blklist
mkdir -p /mnt/xuruizhi/ATAC_brain/human/final
cd /mnt/xuruizhi/ATAC_brain/human/blklist
wget https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gzip -dc hg38.blacklist.bed.gz > hg38.blacklist.bed
rm hg38.blacklist.bed.gz


cd ../filter
cp ../sequence/3.list .

vim human.sh
# local
# rm blklist
bedtools intersect -wa -a {}.filter.bam -b ../blklist/hg38.blacklist.bed | \
wc -l  > ../blklist/{}.intersect.list
bedtools intersect -v -a {}.filter.bam -b ../blklist/hg38.blacklist.bed > ../final/{}.final.bam
samtools index -@ 6 ../final/{}.final.bam
samtools flagstat -@ 6 ../final/{}.final.bam > ../final/{}.final.stat


cat 3.list  | while read id
do 
  sed "s/{}/${id}/g" human.sh > ${id}_final.sh
  bash  ${id}_final.sh >> ./final.log 2>&1
done
```

# 5. 合并neuron和non-neuron

samtools merge要求是对排序后的bam文件进行合并，生成和现在顺序一样的bam文件。需要注意的是，bam文件需要有.bai索引。
```bash
samtools merge [options] -o out.bam [options] in1.bam ... inN.bam
# 将对应文件合并且用脑区命名
```


