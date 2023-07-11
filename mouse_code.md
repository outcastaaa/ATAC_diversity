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


- [10. 不同组织可重复peak](#10-不同组织可重复peak)  
- [9. diffbind](#9-使用diffbind做主成分分析)


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
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse
```

# 2. 下载数据
1. 测序数据
```bash
# mouse 23个，后三个是单端测序SRR13443549,SRR13443553,SRR13443554
cd /mnt/xuruizhi/ATAC_brain/mouse
cat >MOUSE.list <<EOF
SRR11179779    0e6ac756b030cc84f42de5bce499f5fd
SRR11179780    d1766b149fe0182ae872bf79a3decd5a
SRR11179781    b93250c0ce5affeb603fa071c78fcce0
SRR13049359    3f280eaf565605981db8d3062bf6a6a3
SRR13049360    70677096eb974366da7152f36fa0c5b0
SRR13049361    119bfddf4301e49835661bf24693199c
SRR13049362    674996a273c74bbb3d52cfc025e155da
SRR13049363    a54bb8dad247d4e7c114f12fd1b8791a
SRR13049364    bf1dc5b5a62011235fdf854f860c98b8
SRR14362271    8667509648f566de1e8b3b3e3638b3e0
SRR14362272    bb509a5db71ff478426c773aa2c9770d
SRR14362275    5eb17b41b9ac7c2caa4967752c416c98
SRR14362276    075c3926a6a49c15ce85ac9194cc5370
SRR14362281    c2cee7a6f408c5db7c8a81119e2ee43b
SRR14362282    33a21b0f50030658b5a2d48f6573c1c3
SRR14614715×
SRR3595211    fd2862221c839f1da0bf667704eb63c1
SRR3595212    21076454a9730edede42d98f09aa10bc
SRR3595213    c277a3891a38d31a0beb94105937e45d
SRR3595214    b2ca2eee6d312f08831d9e27657dda37
SRR13443549.   0e3302401f4581bddc071b41c11f939a
SRR13443553.   861aa5a62152623cb54d34fd12ede271
SRR13443554.   5db94c5616aaeba4e3970b8f487fd4c1
EOF 

mkdir -p ./sra
cd ./sra
cp /mnt/d/perl/perl_scripts/download_srr.pl ./
cp /mnt/xuruizhi/ATAC_brain/mouse/MOUSE.list ./
# 单样本举例
perl download_srr.pl --output-dir . --srr SRR11179779
# 批量下载
cat MOUSE.list | parallel -k -j 6 "
  echo {} >> ./download.log
  perl download_srr.pl --output-dir . --srr {} >> ./download.log 2>&1
"
cat download.log | grep " downloaded successfully"
```

```bash
# 在本地转换格式
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/sequence
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
cat >MOUSE.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049360
SRR13049361
SRR13049362
SRR13049363
SRR13049364
SRR14362271
SRR14362272
SRR14362275
SRR14362276
SRR14362281
SRR14362282
SRR3595211  
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF

fqdir=/mnt/xuruizhi/ATAC_brain/mouse/sequence
cat MOUSE.list | while read id
do
  sed "s/{}/${id}/g" sra2fq.sh > ${id}_sra2fq.sh
done

cat MOUSE.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_sra2fq.sh"


rsync -av /mnt/xuruizhi/ATAC_brain/mouse/sequence \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/mouse/
# find /scratch/wangq/xrz/ATAC_brain/mouse/sra -type d -name "SRR*"\
# -exec sh -c 'find {} -type f -name "*.sra" -exec cp -r {} ../ \;' \;此代码有错
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
md5sum *.fastq.gz
# e7c21902b5f6228c04bfefc4c7605d3d  SRR11179779_1.fastq.gz
# a17f0de8bb0baf6d9612e6cca491420c  SRR11179779_2.fastq.gz
# 642e3990eca76e041138f0d1dc12f9a1  SRR11179780_1.fastq.gz
# b3cbea2262204638e15b2fb4727f4337  SRR11179780_2.fastq.gz
# b8721636cf9ec8ac56502f5a430fbf62  SRR11179781_1.fastq.gz
# be248c9c98d6457ad4602d83efc85d29  SRR11179781_2.fastq.gz
# 8146a8f87792beaafa66da54106523d8  SRR13049359_1.fastq.gz
# 10bf9afbee98da9265ec19de38c8a2bc  SRR13049359_2.fastq.gz
# 7caf7e99e7948298e6a63e826933537a  SRR13049360_1.fastq.gz
# 30b5cfef2978503998369c9d44674720  SRR13049360_2.fastq.gz
# 4149b45bd8e854d5b93a2a13e62b9529  SRR13049361_1.fastq.gz
# deba00dd254e75092ca3ae2efeabaff6  SRR13049361_2.fastq.gz
# 990cf81d4e9a11d1e8cfb5cf4795db4c  SRR13049362_1.fastq.gz
# 0bba68428f864feb3d0052e6e8024c84  SRR13049362_2.fastq.gz
# 0fefd64e285d298dbf7c9516bf5c15e2  SRR13049363_1.fastq.gz
# b8fbd54a787d43b16b62a9db59c0a5c5  SRR13049363_2.fastq.gz
# 6729a52e1ec94690ac3134e5d0db71f4  SRR13049364_1.fastq.gz
# e878afe75d759c37bc3c9aeb08601d96  SRR13049364_2.fastq.gz
# b2276393d4aff553f4062996c88db88b  SRR13443549.fastq.gz
# 0639ec23552d6ff5f6fd03fdc2c563de  SRR13443553.fastq.gz
# 838ffe4696d2ca4e1545a533bf8dc450  SRR13443554.fastq.gz
# e69b3869609a59aaebaa143904c7edf3  SRR14362271_1.fastq.gz
# a9a1d6d99b5335ebf7f65efe89782a0c  SRR14362271_2.fastq.gz
# 6260ef37666965814587a39792f61e07  SRR14362272_1.fastq.gz
# b2897843281890540530dbfb5e2ece16  SRR14362272_2.fastq.gz
# ff4995efc2f931cc4b35cc56d5ab7324  SRR14362275_1.fastq.gz
# 46a5c6651d1b7e3153eb2976794ff7cc  SRR14362275_2.fastq.gz
# 06202182d3e3c9c32ccc957a440c7176  SRR14362276_1.fastq.gz
# 4c938feea42532eea6aee709f8ca46df  SRR14362276_2.fastq.gz
# c923c0b7ed80e4b8a396e6ab2b33e6ca  SRR14362281_1.fastq.gz
# a3496dffea0ab3eecb2e9253953e302f  SRR14362281_2.fastq.gz
# afc7182bedb6d4a9b93308bba235fa7b  SRR14362282_1.fastq.gz
# e183b449c01a574d6c5fb4861e753230  SRR14362282_2.fastq.gz
# d9c2fff6da973c66e93a9c747c4ce587  SRR3595211_1.fastq.gz
# 2d0cde1499d18225c4f65746a0832b4a  SRR3595211_2.fastq.gz
# 22c9367d651462327e25a01b11cf8bce  SRR3595212_1.fastq.gz
# ba0e915bbb587b562f5cfe4723521893  SRR3595212_2.fastq.gz
# 5a0cb0d75c67239c40031b5046f954b7  SRR3595213_1.fastq.gz
# e568383fa7f65d5ef8d34395012a944b  SRR3595213_2.fastq.gz
# d2704eabde9487fac8b21617afc4ad71  SRR3595214_1.fastq.gz
# 7acc0d714ac4ef99e6b5c37d24321335  SRR3595214_2.fastq.gz
cd /scratch/wangq/xrz/ATAC_brain/mouse/sequence
md5sum *.fastq.gz
# 核对没有问题
```
2. 基因组数据

# 3. 比对前质控
* 本地质控，先不要执行  
```bash 
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/fastqc
# 查看测序质量
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc *.gz
cd /mnt/xuruizhi/ATAC_brain/mouse/fastqc
multiqc .

# 质控
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/trim
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
# 双端测序
ls *_1.fastq.gz | sed 's/_1.fastq.gz//g'  >pair.list
trim_dir=/mnt/xuruizhi/ATAC_brain/mouse/trim
cat pair.list | parallel --no-run-if-empty --linebuffer -k -j 8 ' \
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired \
-o /mnt/xuruizhi/ATAC_brain/mouse/trim  {}_1.fastq.gz  {}_2.fastq.gz'
# 单端测序
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
cat >single.list <<EOF
SRR13443549
SRR13443553
SRR13443554
EOF
cat single.list | while read id
do 
  trim_galore --phred33 --length 35 -e 0.1 --stringency 3 \
  -o /mnt/xuruizhi/ATAC_brain/mouse/trim ${id}.fastq.gz
done
```
* 本地循环，执行
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/trim
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again

# 双端
ls *_1.fastq.gz | sed 's/_1.fastq.gz//g'  >pair.list
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
# 循环脚本
cat >mouse_pair.sh <<EOF
#!/usr/bin bash
# This script is for pari-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/mouse/sequence /mnt/xuruizhi/ATAC_brain/mouse/sra/{}/{}.sra

# fastqc
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_1.fastq.gz
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}_2.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired \
-o /mnt/xuruizhi/ATAC_brain/mouse/trim  {}_1.fastq.gz  {}_2.fastq.gz

# fatsqc_again
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_1_val_1.fq.gz
fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_2_val_2.fq.gz
EOF

cat pair.list | while read id
do
  sed "s/{}/${id}/g" mouse_pair.sh > ${id}_pair_trim.sh
done
cat pair.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_pair_trim.sh >> ./trim_fastqc.log 2>&1"


# 单端
cd /mnt/xuruizhi/ATAC_brain/mouse/sequence
cat >single.list <<EOF
SRR13443549
SRR13443553
SRR13443554
EOF

# 循环脚本
cat >mouse_single.sh <<EOF
#!/usr/bin bash
# This script is for single-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/mouse/sequence /mnt/xuruizhi/ATAC_brain/mouse/sra/{}/{}.sra

# fastqc
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc /mnt/xuruizhi/ATAC_brain/mouse/sequence/{}.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 -o /mnt/xuruizhi/ATAC_brain/mouse/trim {}.fastq.gz

# fatsqc_again
# fastqc -t 6 -o /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again /mnt/xuruizhi/ATAC_brain/mouse/trim/{}_trimmed.fq.gz
EOF

cat single.list | while read id
do
  sed "s/{}/${id}/g" mouse_single.sh > ${id}_single_trim.sh
  bash ${id}_single_trim.sh >> ./trim_fastqc.log 2>&1
done

cd  /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again
multiqc .

# 因为CG含量质控不合格，删除了
SRR11179779_1_val_1
SRR11179779_2_val_2
SRR13049360_1_val_1
SRR13049360_2_val_2
SRR13049361_1_val_1
SRR13049361_2_val_2
```
```bash
# 传到超算
# sequence后来的补充内容没有上传到超算
rsync -av /mnt/xuruizhi/ATAC_brain/mouse/trim \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/mouse/
rsync -av /mnt/xuruizhi/ATAC_brain/mouse/fastqc \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/mouse/
rsync -av /mnt/xuruizhi/ATAC_brain/mouse/fastqc_again \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/mouse/
```


# 4. 比对

1. 参考基因组  
```bash
# 小鼠 mm10 的 bowtie2 的 index 已经建立过，在超算中处理
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/genome
cp /mnt/d/atac/genome/* /mnt/xuruizhi/ATAC_brain/mouse/genome/
# 进入超算
mkdir -p /scratch/wangq/xrz/ATAC_brain/mouse/genome
cp /share/home/wangq/xuruizhi/brain/brain/genome/mouse/* /scratch/wangq/xrz/ATAC_brain/mouse/genome/

cd /scratch/wangq/xrz/ATAC_brain/mouse/trim
vim pair.list
SRR11179780
SRR11179781
SRR13049359
SRR13049362
SRR13049363
SRR13049364
SRR14362271
SRR14362272
SRR14362275
SRR14362276
SRR14362281
SRR14362282
SRR3595211
SRR3595212
SRR3595213
SRR3595214

vim single.list
SRR13443549
SRR13443553
SRR13443554
```

2. 比对
```bash
mkdir -p /scratch/wangq/xrz/ATAC_brain/mouse/align/
# 循环 
cd /scratch/wangq/xrz/ATAC_brain/mouse/trim

# 双端
bowtie2  -p 48 -x /scratch/wangq/xrz/ATAC_brain/mouse/genome/mm10 \
--very-sensitive -X 2000 -1 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_1_trimmed.fq.gz \
-2 /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_2_trimmed.fq.gz \
-S /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam

# 单端
bowtie2  -p 48 -x /scratch/wangq/xrz/ATAC_brain/mouse/genome/mm10 \
--very-sensitive -X 2000 -U /scratch/wangq/xrz/ATAC_brain/mouse/trim/{}_trimmed.fq.gz \
-S /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam
```

3. sort_transfertobam_index  
```bash
mkdir -p /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam
# 双端与单端一致
samtools sort -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/align/{}.sam \
> /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools index -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam
samtools flagstat  -@ 48 /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.sort.bam \
> /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam/{}.raw.stat
```
4. 大批量处理
```bash
cd /scratch/wangq/xrz/ATAC_brain/mouse/trim
# 双端 
cat pair.list  | while read id
do 
  sed "s/{}/${id}/g" mouse_pair.sh > ${id}_pair_align.sh
  bsub -q mpi -n 96 -o ../align "
  bash  ${id}_pair_align.sh >> ../align/align.log 2>&1"
done
# 单端
cat single.list  | while read id
do 
  sed "s/{}/${id}/g" mouse_single.sh > ${id}_single_align.sh;
  bsub -q largemem -n 96 -o ../align "
  bash ${id}_single_align.sh >> ../align/align_single.log 2>&1"
done
```

4. 传到本地
```bash
# 通过beyond compare传输
# 太大了，先不传输
```
5. 比对结果
```bash
# 比对率都可以在90以上
cat SRR14362282.raw.stat
92865004 + 0 in total (QC-passed reads + QC-failed reads)
92865004 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
91508313 + 0 mapped (98.54% : N/A)
91508313 + 0 primary mapped (98.54% : N/A)
92865004 + 0 paired in sequencing
46432502 + 0 read1
46432502 + 0 read2
90382150 + 0 properly paired (97.33% : N/A)
90816946 + 0 with itself and mate mapped
691367 + 0 singletons (0.74% : N/A)
29134 + 0 with mate mapped to a different chr
8432 + 0 with mate mapped to a different chr (mapQ>=5)
```
# 5. post-alignment
```bash
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam
# 下面list删掉了质控不好的
cat >MOUSE.list <<EOF
SRR11179780
SRR11179781
SRR13049359
SRR13049362
SRR13049363
SRR13049364
SRR14362271
SRR14362272
SRR14362275
SRR14362276
SRR14362281
SRR14362282
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF
```
1. remove PCR-duplicate reads
```bash
mkdir -p /scratch/wangq/xrz/ATAC_brain/mouse/rmdup
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam

vim mouse_common.sh
#!/usr/bin bash

# rmdup
parallel -k -j 96 'picard MarkDuplicates -I ./{}.sort.bam \
-O ../rmdup/{}.rmdup.bam \
-REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT \
-METRICS_FILE ../rmdup/{}.log'

# index
samtools index -@ 96 ../rmdup/{}.rmdup.bam
samtools flagstat -@ 96 ../rmdup/{}.rmdup.bam > ../rmdup/${sample}.rmdup.stat
```


2. remove bad quality reads and chrM reads
```bash
mkdir -p /scratch/wangq/xrz/ATAC_brain/mouse/filter
cd /scratch/wangq/xrz/ATAC_brain/mouse/rmdup
cp ../trim/*.list ./
# 双端
vim mouse_pair.sh
samtools view -h -f 2 -F 1804 -q 30 ../rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 48 -O bam  -o ../filter/{}.filter.bam
samtools index -@ 48 ../filter/{}.filter.bam
samtools flagstat -@ 48 ../filter/{}.filter.bam > ../filter/{}.filter.stat

# 单端
samtools view -h -F 1804 -q 30 ../rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 48 -O bam  -o ../filter/{}.filter.bam
samtools index -@ 48 ../filter/{}.filter.bam
samtools flagstat -@ 48 ../filter/{}.filter.bam > ../filter/{}.filter.stat
```


3. 大批量处理
```bash
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam
# 因为picard提交任务报错，直接跑
picard MarkDuplicates -I ./SRR3595211.sort.bam -O ../rmdup/SRR3595211.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE ../rmdup/SRR3595211.log


cd /scratch/wangq/xrz/ATAC_brain/mouse/rmdup
cat pair.list | while read id
do
  sed "s/{}/${id}/g" mouse_pair.sh > ${id}_filter.sh
  bsub -q mpi -n 48 -o ../filter "
  bash ${id}_filter.sh >> ../filter/filter.log 2>&1"
done
cat single.list | while read id
do
  sed "s/{}/${id}/g" mouse_single.sh > ${id}_filter.sh
  bsub -q mpi -n 48 -o ../filter "
  bash ${id}_filter.sh >> ../filter/filter.log 2>&1"
done
```


4. Blacklist filtering

```bash
# beyond compare传输到本地
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/blklist
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/final
# 下载对应物种的 blacklist.bed文件
wc -l  mm10.blacklist.bed #164
cp /mnt/d/ATAC/blklist/mm10.blacklist.bed /mnt/xuruizhi/ATAC_brain/mouse/blklist

cd /mnt/xuruizhi/ATAC_brain/mouse/filter
cat >MOUSE.list <<EOF
SRR11179780
SRR11179781
SRR13049359
SRR13049362
SRR13049363
SRR13049364
SRR14362271
SRR14362272
SRR14362275
SRR14362276
SRR14362281
SRR14362282
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF

# 取交集看bam文件和blacklist有多少重合部分
# 凡是bam中含有blacklist都删除
cat MOUSE.list | parallel -k -j 6 "
  echo {} 
  echo "{}.filter.bam"
  bedtools intersect -wa -a {}.filter.bam -b ../blklist/mm10.blacklist.bed | \
  wc -l  > ../blklist/{}.intersect.list

  bedtools intersect -v -a {}.filter.bam -b ../blklist/mm10.blacklist.bed > ../final/{}.final.bam
  samtools index ../final/{}.final.bam
  samtools flagstat ../final/{}.final.bam > ../final/{}.final.stat
"
```
5. 结果解读：  
```bash
# 原比对文件数据
cat SRR14362282.raw.stat
92865004 + 0 in total (QC-passed reads + QC-failed reads)
92865004 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
91508313 + 0 mapped (98.54% : N/A)
91508313 + 0 primary mapped (98.54% : N/A)
92865004 + 0 paired in sequencing
46432502 + 0 read1
46432502 + 0 read2
90382150 + 0 properly paired (97.33% : N/A)
90816946 + 0 with itself and mate mapped
691367 + 0 singletons (0.74% : N/A)
29134 + 0 with mate mapped to a different chr
8432 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除PCR重复+低质量+chrM后数据
# rmdup
72786670 + 0 in total (QC-passed reads + QC-failed reads)
72786670 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
71429979 + 0 mapped (98.14% : N/A)
71429979 + 0 primary mapped (98.14% : N/A)
72786670 + 0 paired in sequencing
36359425 + 0 read1
36427245 + 0 read2
70627278 + 0 properly paired (97.03% : N/A)
70971626 + 0 with itself and mate mapped
458353 + 0 singletons (0.63% : N/A)
25462 + 0 with mate mapped to a different chr
8250 + 0 with mate mapped to a different chr (mapQ>=5)
# filter
61465170 + 0 in total (QC-passed reads + QC-failed reads)
61465170 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
61465170 + 0 mapped (100.00% : N/A)
61465170 + 0 primary mapped (100.00% : N/A)
61465170 + 0 paired in sequencing
30732585 + 0 read1
30732585 + 0 read2
61465170 + 0 properly paired (100.00% : N/A)
61465170 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除blacklist后数据
61430470 + 0 in total (QC-passed reads + QC-failed reads)
61430470 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
61430470 + 0 mapped (100.00% : N/A)
61430470 + 0 primary mapped (100.00% : N/A)
61430470 + 0 paired in sequencing
30715234 + 0 read1
30715236 + 0 read2
61430470 + 0 properly paired (100.00% : N/A)
61430470 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
到这一步，比对文件已经过滤完成。   

# 6. Call peaks 
1. bamtobed
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/bed
cd /mnt/xuruizhi/ATAC_brain/mouse/final
bedtools bamtobed -i {}.final.bam > ../bed/{}.bed
```

2. shift reads
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/Tn5_shift
cd /mnt/xuruizhi/ATAC_brain/mouse/bed

cat ../final/MOUSE.list | while read id;
do 
  echo $id 
  cat ${id}.bed | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/${id}.Tn5.bed
done
```
3. Call peaks 
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/peaks
cd /mnt/xuruizhi/ATAC_brain/mouse/Tn5_shift
macs2 callpeak  -g mm \
  --shift -75 --extsize 150 --nomodel \
  --nolambda --keep-dup all \
  -n {} -t ./{}.Tn5.bed \
  --outdir ../peaks/
```
4. 批量处理
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/bed
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/Tn5_shift
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/peaks

cd /mnt/xuruizhi/ATAC_brain/mouse/final
cp /mnt/xuruizhi/ATAC_brain/mouse/filter/MOUSE.list ../final/
cat MOUSE.list | while read id
do
  sed "s/{}/${id}/g" mouse_common.sh > ../peaks/${id}_callpeaks.sh
done
cat MOUSE.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash ../peaks/{}_callpeaks.sh >> ../peaks/peaks.log 2>&1"
```

# 7. Quality check  
判断ATAC-seq是否合格的几个[Current Standards](https://www.encodeproject.org/atac-seq/)  

* 实验重复 Experiments should have two or more biological replicates. √
* 数据量 Each replicate should have 25 million non-duplicate, non-mitochondrial aligned reads for single-end sequencing and 50 million for paired-ended sequencing (i.e. 25 million fragments, regardless of sequencing run type). 
```bash
cat *.final.stat | grep 'read1'
# SRR11179780.final.stat
# 55174171 + 0 in total (QC-passed reads + QC-failed reads)
# SRR11179781.final.stat    reads不够
# 35179699 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13049359.final.stat    reads不够
# 11652432 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13049362.final.stat
# 57619495 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13049363.final.stat
# 64637521 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13049364.final.stat    reads不够
# 43270057 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13443549.final.stat
# 50731638 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13443553.final.stat
# 50592991 + 0 in total (QC-passed reads + QC-failed reads)
# SRR13443554.final.stat    reads不够
# 47570952 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362271.final.stat    reads不够
# 43903760 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362272.final.stat
# 97205176 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362275.final.stat
# 67693822 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362276.final.stat
# 78843819 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362281.final.stat
# 62889953 + 0 in total (QC-passed reads + QC-failed reads)
# SRR14362282.final.stat
# 61430470 + 0 in total (QC-passed reads + QC-failed reads)
# SRR3595211.final.stat    reads不够
# 10203178 + 0 in total (QC-passed reads + QC-failed reads)
# SRR3595212.final.stat    reads不够
# 14113255 + 0 in total (QC-passed reads + QC-failed reads)
# SRR3595213.final.stat    reads不够
# 5237407 + 0 in total (QC-passed reads + QC-failed reads)
# SRR3595214.final.stat    reads不够
# 7013056 + 0 in total (QC-passed reads + QC-failed reads)
```
* 比对率 The alignment rate, or percentage of mapped reads, should be greater than 95%, though values >80% may be acceptable. √
* IDR value Replicate concordance is measured by calculating IDR values (Irreproducible Discovery Rate). The experiment passes if both rescue and self consistency ratios are less than 2.
* Various peak files must meet certain requirements. Please visit the section on output files under the pipeline overview for more information on peak files.  
** The number of peaks within a replicated peak file should be >150,000, though values >100,000 may be acceptable.  
```bash
#  差不多
cd /mnt/xuruizhi/ATAC_brain/mouse/peaks
wc -l *.narrowPeak
#    126857 SRR11179780_peaks.narrowPeak
#    111102 SRR11179781_peaks.narrowPeak
#    174276 SRR13049359_peaks.narrowPeak
#    245020 SRR13049362_peaks.narrowPeak
#    193402 SRR13049363_peaks.narrowPeak
#    172142 SRR13049364_peaks.narrowPeak
#     95584 SRR13443549_peaks.narrowPeak
#    181380 SRR13443553_peaks.narrowPeak
#     80903 SRR13443554_peaks.narrowPeak
#    197143 SRR14362271_peaks.narrowPeak
#    230309 SRR14362272_peaks.narrowPeak
#    206593 SRR14362275_peaks.narrowPeak
#    203466 SRR14362276_peaks.narrowPeak
#    199754 SRR14362281_peaks.narrowPeak
#    169213 SRR14362282_peaks.narrowPeak
#    101655 SRR3595211_peaks.narrowPeak
#    117018 SRR3595212_peaks.narrowPeak
#     58525 SRR3595213_peaks.narrowPeak
#     80972 SRR3595214_peaks.narrowPeak
```

** The number of peaks within an IDR peak file should be >70,000, though values >50,000 may be acceptable.  
** A nucleosome free region (NFR) must be present.  
** A mononucleosome peak must be present in the fragment length distribution. These are reads that span a single nucleosome, so they are longer than 147 bp but shorter than 147*2 bp. Good ATAC-seq datasets have reads that span nucleosomes (which allows for calling nucleosome positions in addition to open regions of chromatin).  


* The fraction of reads in called peak regions (FRiP score) should be >0.3, though values greater than 0.2 are acceptable. TSS enrichment remains in place as a key signal to noise measure.

* Transcription start site (TSS) enrichment values are dependent on the reference files used; cutoff values for high quality data are listed in the table below.  


| Annotation used              | Value  | Resulting Data Status  |
|------------------------------|--------|------------------------|
| mm10 Refseq TSS annotation   | < 10   | Concerning             |
| mm10 Refseq TSS annotation   | 10--15 | Acceptable             |
| mm10 Refseq TSS annotation   | > 15   | Ideal                  |

1. fragment length distribution 一般只对于双端测序有意义
```bash
# 在Linux中画图  
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/frag_length
cd /mnt/xuruizhi/ATAC_brain/mouse/final

cat MOUSE.list | parallel -k -j 6 "
  echo {} 
  java -jar /mnt/d/biosoft/picard/picard.jar CollectInsertSizeMetrics \
  -I {}.final.bam \
  -O ../frag_length/{}.insert_size_metrics.txt \
  -H ../frag_length/{}.insert_size_histogram.pdf"
```


2.  FRiP 影响不大
FRiP（Fraction of reads in peaks，Fraction of all mapped reads that fall into the called peak regions）表示的是位于peak区域的reads的比例，FRiP score是一个比值，其分子是位于peak区域的reads总数，分母是比对到参考基因组上的reads总数。  

* 推荐使用未去重、只比对完后的bam文件，后续使用diffbind来看

* 输出结果错误
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/FRiP
cd /mnt/xuruizhi/ATAC_brain/mouse/Tn5_shift
cp ../final/MOUSE.list ./
cat MOUSE.list | parallel -k -j 6 "
  wc -l {}.Tn5.bed | awk '{print $1}' >> ../FRiP/bed_totalReads.txt
  bedtools intersect -wa -a {}.Tn5.bed -b ../peaks/{}_peaks.narrowPeak \
  | wc -l | awk '{print $1}' >> ../FRiP/bed_peakReads.txt

  paste MOUSE.list ../FRiP/bed_totalReads.txt ../FRiP/bed_peakReads.txt > ../FRiP/bed_FRiP.txt
  cat ../FRiP/bed_FRiP.txt | awk '{print $1,$2,$3,$2/$3*100"%"}' > ../FRiP/bed_FRiP.txt
"
```

3. 统计平均长度与长度分布情况
* bash统计平均值
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/length
cd /mnt/xuruizhi/ATAC_brain/mouse/peaks

for i in *_peaks.narrowPeak
do
  echo $i
  cat $i | awk '{print $3-$2}' | awk '{sum += $1} END {print avg = sum/NR}'
done
# SRR11179780_peaks.narrowPeak
# 470.198
# SRR11179781_peaks.narrowPeak
# 392.47
# SRR13049359_peaks.narrowPeak
# 363.649
# SRR13049362_peaks.narrowPeak
# 425.259
# SRR13049363_peaks.narrowPeak
# 411.096
# SRR13049364_peaks.narrowPeak
# 391.281
# SRR13443549_peaks.narrowPeak
# 454.323
# SRR13443553_peaks.narrowPeak
# 427.244
# SRR13443554_peaks.narrowPeak
# 480.342
# SRR14362271_peaks.narrowPeak
# 412.807
# SRR14362272_peaks.narrowPeak
# 462.077
# SRR14362275_peaks.narrowPeak
# 457.008
# SRR14362276_peaks.narrowPeak
# 453.776
# SRR14362281_peaks.narrowPeak
# 480.939
# SRR14362282_peaks.narrowPeak
# 439.371
# SRR3595211_peaks.narrowPeak
# 339.392
# SRR3595212_peaks.narrowPeak
# 380.324
# SRR3595213_peaks.narrowPeak
# 286.26
# SRR3595214_peaks.narrowPeak
# 288.26
```
* R中统计与画图
！ 将R储存在 D:/ATAC_brain文件夹中
```r

> summary(a)
       V1       
 Min.   :  167  
 1st Qu.:  696  
 Median :  995  
 Mean   : 1110  
 3rd Qu.: 1384  
 Max.   :11927

```
```r
# HIPP长度分布情况
# 画图
getwd()    #[1] "D:/brain/brain/R_analyse"
# 以HIPP为例
a<-read.table('../IDR_final/mouse/HIPP_pool_merge_peak_length.txt')
dim(a)
png('hist.png')
hist(abs(as.numeric(a[,1])),breaks=500,xlab = "Fragment length(bp)",ylab = "Frequency",main = "HIPP peak length")

# 其他样本类似
b <-read.table('../IDR_final/mouse/PFC_pool_merge_peak_length.txt')
hist(abs(as.numeric(b[,1])),breaks=500,xlab = "Fragment length(bp)",ylab = "Frequency",main = "PFC peak length")
c <-read.table('../IDR_final/mouse/cortex_pool_merge_peak_length.txt')
hist(abs(as.numeric(c[,1])),breaks=500,xlab = "Fragment length(bp)",ylab = "Frequency",main = "cortex peak length")
```
4. TSS enrichment在下一节展示，很重要的判断质量的标准

# 8. Visualization    
1.  filterbam2Bw    
```bash 
mkdir -p  /mnt/xuruizhi/ATAC_brain/mouse/bw
cd /mnt/xuruizhi/ATAC_brain/mouse/final
conda activate py3.8
ls *.bam | while read id; 
do 
  bamCoverage -p 6  -b $id \
  -o ../bw/${id%%.*}.bw \
  --binSize 20 \
  --smoothLength 60 \
  --normalizeUsing RPKM \
  --centerReads 
  >> ../bw/bamCoverage.log 2>&1
done
```

2. TSS enrichment 

下载TSS注释文件：the BED file which contains the coordinates for all genes [下载地址](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=6884883_WoMR8YyIAAVII92Rr1Am3Kd0jr5H&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr12%3A56703576-56703740&hgta_outputType=primaryTable&hgta_outFileName=)   
[参数选择](https://www.jianshu.com/p/d6cb795af22a)   

将`mm10.reseq.bed`保存在 /mnt/d/ATAC/TSS 文件夹内。 
```bash
# 在py3.8环境下运行
conda activate py3.8
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/TSS
cp /mnt/d/ATAC/TSS/mm10.refseq.bed /mnt/xuruizhi/ATAC_brain/mouse/TSS

cd /mnt/xuruizhi/ATAC_brain/mouse/bw
ls *.bw | while read id; 
do 
  computeMatrix reference-point --referencePoint TSS -p 6 \
    -b 1000  -a 1000 \
    -R ../TSS/mm10.refseq.bed \
    -S $id \
    --skipZeros \
    -o ../TSS/${id%%.*}_matrix.gz \
    --outFileSortedRegions ../TSS/${id%%.*}_regions.bed
done


# profile plot
cd ../TSS/
cat ../final/MOUSE.list | while read id
do 
  plotProfile -m ${id}_matrix.gz \
    -out ${id}_profile.png \
    --perGroup \
    --colors green \
    --plotTitle "" \
    --refPointLabel "TSS" \
    -T "${id} read density" \
    -z ""
done

# heatmap
cat ../final/MOUSE.list | while read id
do 
  plotHeatmap -m ${id}_matrix.gz \
  -out ${id}_heatmap.png \
  --colorMap RdBu \
  --whatToShow 'heatmap and colorbar' \
  --zMin -8 --zMax 8  
done


# heatmap and profile plot
cat ../final/MOUSE.list | while read id
do 
  plotHeatmap -m ${id}_matrix.gz \
    -out ${id}_all.png \
    --colorMap RdBu \
    --zMin -12 --zMax 12
done


# 画 `gene body` 区，使用 `scale-regions`  
cd ../bw
mkdir -p ../genebody
# create a matrix 
ls *.bw | while read id; 
do
computeMatrix scale-regions -p 6 \
    -b 10000  -a 10000 \
    -R ../TSS/mm10.refseq.bed \
    -S ${id} \
    --skipZeros \
    -o ../genebody/${id%%.*}_matrix.gz 
done

ls ../genebody/*.gz | while read id
do
  plotHeatmap -m ${id} -out ${id%%.*}_heatmap.png 
done
```

# 9. 寻找rep间consensus peak
1. 对样本进行整理
* 见"mouse数据质量情况.xlsx文件"
```bash
# Sort peak by -log10(p-value)
mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/IDR
cd /mnt/xuruizhi/ATAC_brain/mouse/peaks
cat ../final/MOUSE.list | parallel -j 6 -k "
sort -k8,8nr {}_peaks.narrowPeak > ../IDR/{}.narrowPeak 
"

cd ../IDR

vim idr.sh
#!/bin/bash
# 获取输入参数
input_files=("$@")  # 从命令行参数获取输入文件列表
output_dir="${input_files[-1]}"  # 最后一个参数作为输出目录
unset 'input_files[${#input_files[@]}-1]'  # 移除输出目录参数

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"

# 函数：运行IDR分析
run_idr_analysis() {
    local input1="$1"
    local input2="$2"
    local output_name="$3"
    
    idr --samples "$input1" "$input2" \
        --input-file-type narrowPeak \
        --rank p.value \
        --soft-idr-threshold 0.05 \
        --use-best-multisummit-IDR \
        --output-file "${output_dir}/${output_name}_0.05.txt" \
        --log-output-file "${output_dir}/${output_name}_0.05.log" \
        --plot
    
    echo "IDR analysis completed for samples: $input1, $input2"
    echo "Output files: ${output_dir}/${output_name}_0.05.txt, ${output_dir}/${output_name}_0.05.log"
}

# 运行IDR分析
for ((i=0; i<${#input_files[@]}; i++)); do
    for ((j=i+1; j<${#input_files[@]}; j++)); do
        file1="${input_files[i]}"
        file2="${input_files[j]}"
        file1_name=$(basename "$file1" | sed 's/\..*//')  # 提取文件名，去除文件扩展名
        file2_name=$(basename "$file2" | sed 's/\..*//')  # 提取文件名，去除文件扩展名
        output_name="${file1_name}_${file2_name}"
        
        run_idr_analysis "$file1" "$file2" "$output_name"
    done
done
```
2. 代码
① HIPP: SRR11179780 + SRR11179781
```bash
cd /mnt/xuruizhi/ATAC_brain/mouse/IDR
chmod +x idr.sh 
./idr.sh SRR11179780.narrowPeak SRR11179781.narrowPeak .
mv SRR11179780_SRR11179781.txt HIPP.txt # 66547

awk '{if($5 >= 540) print $0}' HIPP.txt > HIPP_common0.05.txt
cut -f 1,2,3 HIPP_common0.05.txt > HIPP_common0.05.bed
cat HIPP_common0.05.bed | tsv-summarize -g 1 --count
cat HIPP_common0.05.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "chrX_GL456233_random" | grep -v "chr4_GL456216_random" > HIPP_pool.bed
wc -l *.bed
#  14761 HIPP_common0.05.bed
#  14723 HIPP_pool.bed
sort -k1,1 -k2,2n HIPP_pool.bed > HIPP_pool_sort.bed
bedtools merge -i HIPP_pool_sort.bed -d 50 > HIPP_pool_merge.bed
wc -l  HIPP_pool_merge.bed
#  14720 HIPP_pool_merge.bed
```

② cortex: SRR13049359 + SRR13049362
```bash
cd /mnt/xuruizhi/ATAC_brain/mouse/IDR 
./idr.sh SRR13049359.narrowPeak SRR13049362.narrowPeak .
mv SRR13049359_SRR13049362.txt cortex.txt 
# 140856

awk '{if($5 >= 540) print $0}' cortex.txt > cortex_common0.05.txt
cut -f 1,2,3 cortex_common0.05.txt > cortex_common0.05.bed
cat cortex_common0.05.bed | tsv-summarize -g 1 --count
cat cortex_common0.05.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "chrX_GL456233_random" | grep -v "chr4_GL456216_random" | grep -v "chr4_JH584295_random" > cortex_pool.bed
wc -l cortex*.bed
#   41485 cortex_common0.05.bed
#   41470 cortex_pool.bed
sort -k1,1 -k2,2n cortex_pool.bed > cortex_pool_sort.bed
bedtools merge -i cortex_pool_sort.bed -d 50 > cortex_pool_merge.bed
wc -l  cortex_pool_merge.bed
#  41455 
```
③ STR: SRR13049363 + SRR13049364
```bash
cd /mnt/xuruizhi/ATAC_brain/mouse/IDR 
./idr.sh SRR13049363.narrowPeak SRR13049364.narrowPeak .
mv SRR13049363_SRR13049364.txt STR.txt 
# 128132

awk '{if($5 >= 540) print $0}' STR.txt > STR_common0.05.txt
cut -f 1,2,3 STR_common0.05.txt > STR_common0.05.bed
cat STR_common0.05.bed | tsv-summarize -g 1 --count
cat STR_common0.05.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "chrX_GL456233_random" | grep -v "chr4_GL456216_random" | grep -v "chr4_JH584295_random" > STR_pool.bed
wc -l STR*.bed
#   53187 STR_common0.05.bed
#   53166 STR_pool.bed
sort -k1,1 -k2,2n STR_pool.bed > STR_pool_sort.bed
bedtools merge -i STR_pool_sort.bed -d 50 > STR_pool_merge.bed
wc -l  STR_pool_merge.bed
#  53130 
```
④ PFC：SRR14362271 + SRR14362272 + SRR14362275 + SRR14362276 + SRR14362281 + SRR14362282

```bash
cd /mnt/xuruizhi/ATAC_brain/mouse/IDR 
./idr.sh SRR14362271.narrowPeak SRR14362272.narrowPeak SRR14362275.narrowPeak SRR14362276.narrowPeak SRR14362281.narrowPeak SRR14362282.narrowPeak .
mv SRR14362271_SRR14362272.txt PFC1.txt 
mv SRR14362271_SRR14362275.txt PFC2.txt
mv SRR14362271_SRR14362276.txt PFC3.txt 
mv SRR14362271_SRR14362281.txt PFC4.txt
mv SRR14362271_SRR14362282.txt PFC5.txt 
mv SRR14362272_SRR14362275.txt PFC6.txt
mv SRR14362272_SRR14362276.txt PFC7.txt 
mv SRR14362272_SRR14362281.txt PFC8.txt
mv SRR14362272_SRR14362282.txt PFC9.txt 
mv SRR14362275_SRR14362276.txt PFC10.txt
mv SRR14362275_SRR14362281.txt PFC11.txt 
mv SRR14362275_SRR14362282.txt PFC12.txt
mv SRR14362276_SRR14362281.txt PFC13.txt 
mv SRR14362276_SRR14362282.txt PFC14.txt
mv SRR14362281_SRR14362282.txt PFC15.txt 


for i in PFC*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common0.05.txt"
done

for i in PFC*_common0.05.txt; do
  cut -f 1,2,3 "$i" > "${i%%.*}.bed"
done
cat PFC*_common0.bed | tsv-summarize -g 1 --count
cat PFC*_common0.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "chrX_GL456233_random" | grep -v "chr4_GL456216_random" | grep -v "chr4_JH584295_random" > PFC_pool.bed
wc -l PFC*.bed
#    54800 PFC10_common0.bed
#    56601 PFC11_common0.bed
#    48761 PFC12_common0.bed
#    59358 PFC13_common0.bed
#    50197 PFC14_common0.bed
#    48530 PFC15_common0.bed
#    54767 PFC1_common0.bed
#    51428 PFC2_common0.bed
#    50577 PFC3_common0.bed
#    50965 PFC4_common0.bed
#    45747 PFC5_common0.bed
#    60667 PFC6_common0.bed
#    64515 PFC7_common0.bed
#    66529 PFC8_common0.bed
#    53607 PFC9_common0.bed
#   816562 PFC_pool.bed

sort -k1,1 -k2,2n PFC_pool.bed > PFC_pool_sort.bed
bedtools merge -i PFC_pool_sort.bed -d 50 > PFC_pool_merge.bed
wc -l  PFC_pool_merge.bed
#  91720 
```

⑤ DG：SRR3595211 + SRR3595212 + SRR3595213 + SRR3595214

```bash
cd /mnt/xuruizhi/ATAC_brain/mouse/IDR 
./idr.sh SRR3595211.narrowPeak SRR3595212.narrowPeak SRR3595213.narrowPeak SRR3595214.narrowPeak .
mv SRR3595211_SRR3595212.txt DG1.txt 
mv SRR3595211_SRR3595213.txt DG2.txt
mv SRR3595211_SRR3595214.txt DG3.txt 
mv SRR3595212_SRR3595213.txt DG4.txt
mv SRR3595212_SRR3595214.txt DG5.txt 
mv SRR3595213_SRR3595214.txt DG6.txt


for i in DG*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common0.05.txt"
done

for i in DG*_common0.05.txt; do
  cut -f 1,2,3 "$i" > "${i%%.*}.bed"
done
cat DG*_common0.bed | tsv-summarize -g 1 --count
cat DG*_common0.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "chrX_GL456233_random" | grep -v "chr4_GL456216_random"  > DG_pool.bed
wc -l DG*.bed
#   17444 DG1_common0.bed
#    8474 DG2_common0.bed
#   10532 DG3_common0.bed
#    8766 DG4_common0.bed
#   11193 DG5_common0.bed
#    8682 DG6_common0.bed
#   65051 DG_pool.bed

sort -k1,1 -k2,2n DG_pool.bed > DG_pool_sort.bed
bedtools merge -i DG_pool_sort.bed -d 50 > DG_pool_merge.bed
wc -l  DG_pool_merge.bed
#  21918

































