- [0. 安装下载](#0-安装下载)
- [1. 建目录](#1-建目录)
- [2. 下载数据](#2-下载数据)
- [3. 比对前质控](#3-比对前质控)
- [4. 比对](#4-比对)
- [5. 比对后过滤](#5-post-alignment-processing)
- [6. Tn5转换](#6-shift-reads)
- [7. call peak](#7-call-peaks)
- [8. 可视化](#8-visualization)
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
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam

vim mouse_common.sh
#!/usr/bin bash

samtools view -h -f 2 -q 30 ../rmdup/{}.rmdup.bam | grep -v  chrM | samtools sort -@ 96 -O bam  -o ../filter/{}.filter.bam
samtools index -@ 96 ../filter/{}.filter.bam
samtools flagstat -@ 96 ../filter/{}.filter.bam > ../filter/{}.filter.stat
```


3. 大批量处理
```bash
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam

# 因为picard提交任务报错，直接跑
picard MarkDuplicates -I ./SRR3595211.sort.bam -O ../rmdup/SRR3595211.rmdup.bam  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE ../rmdup/SRR3595211.log

cp MOUSE.list ../rmdup/
cat MOUSE.list | while read id
do
  sed "s/{}/${id}/g" mouse_common.sh > ${id}_filter.sh
  bsub -q mpi -n 48 -o ../filter "
  bash ${id}_filter.sh >> ../filter/filter.log 2>&1"
done
```


4. Blacklist filtering

```bash

cat name_new.list  | while read id; do sed "s/{}/${id}/g" filter.sh > ${id}_filter.sh; done

cat name_new.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/filter_new/mouse "bash ${id}_filter.sh"
done
# Job <8604250-263> is submitted to queue <mpi>.
```





























批量提交
```bash
rmdup_dir=/scratch/wangq/xrz/ATAC_brain/mouse/rmdup
cd /scratch/wangq/xrz/ATAC_brain/mouse/sort_bam

cat name.list  | while read id; do sed "s/{}/${id}/g" rmdup.sh > ${id}_rmdup.sh; done
cat name.list | while read id
do
  bsub -q largemem -n 48 -J rmdup -o $rmdup_dir "bash ${id}.sh"
done
```