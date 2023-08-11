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



vim 2.list 先不处理，因为此脑区不是很重要
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
cat HIPP.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ${id}_qc_trim.sh
done
cat HIPP.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_qc_trim.sh >> ./HIPP.list.log 2>&1"

cd ../fastqc
multiqc .

cd  ../fastqc_again
multiqc .
# 因为质控不合格，去除203&204
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


cat HIPP.list  | while read id
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
  picard MarkDuplicates -I /scratch/wangq/xrz/ATAC_brain/human/sort_bam/{}.sort.bam -O /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.rmdup.bam -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE /scratch/wangq/xrz/ATAC_brain/human/rmdup/{}.log" :::: HIPP.list
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
cat ../align/HIPP.list  | while read id
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

vim human.sh
# local
# rm blklist
bedtools intersect -wa -a {}.filter.bam -b ../blklist/hg38.blacklist.bed | \
wc -l  > ../blklist/{}.intersect.list
bedtools intersect -v -a {}.filter.bam -b ../blklist/hg38.blacklist.bed > ../final/{}.final.bam
samtools index -@ 6 ../final/{}.final.bam
samtools flagstat -@ 6 ../final/{}.final.bam > ../final/{}.final.stat


cat ../sequence/1.list  | while read id
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
```bash
cd /mnt/xuruizhi/ATAC_brain/human/final
# 3.list
samtools merge -@ 4 -o HAB_rep1.bam SRR21163263.final.bam SRR21163264.final.bam
samtools merge -@ 4 -o HAB_rep2.bam SRR21163339.final.bam SRR21163340.final.bam
samtools merge -@ 4 -o PSC_rep1.bam SRR21163311.final.bam SRR21163312.final.bam
samtools merge -@ 4 -o PSC_rep2.bam SRR21163357.final.bam SRR21163358.final.bam
samtools merge -@ 4 -o PSC_rep3.bam SRR21163371.final.bam SRR21163372.final.bam
```

# 6. 对于合并和未合并都Call peaks 

```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/bed
mkdir -p /mnt/xuruizhi/ATAC_brain/human/Tn5_shift
mkdir -p /mnt/xuruizhi/ATAC_brain/human/peaks

cd /mnt/xuruizhi/ATAC_brain/human/final
vim 1.list
PSM
VLPFC
CERE

vim 2.list
PAC
PVC

vim 3.list
HAB_rep1
HAB_rep2
PSC_rep1
PSC_rep2
PSC_rep3

vim 4.list
OFC
NACC
CN

vim 5.list
AMY
DLPFC
HIPP


cat 3.list | while read id
do
  samtools flagstat -@ 6 ${id}.bam > ${id}.stat
done



vim human.sh
# shift bed
# cd /mnt/xuruizhi/ATAC_brain/human/final
bedtools bamtobed -i {}.bam > ../bed/{}.bed
cat ../bed/{}.bed | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/{}.Tn5.bed
macs2 callpeak  -g mm --shift -75 --extsize 150 --nomodel \
--nolambda --keep-dup all -n {} -t ../Tn5_shift/{}.Tn5.bed --outdir ../peaks/


cat 3.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ../peaks/${id}_callpeaks.sh
done
cat 3.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash ../peaks/{}_callpeaks.sh >> ../peaks/peaks.log 2>&1"
```





































# 12. 对应RNA-seq数据  

```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human
```


## 12.1 下载数据

```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/

cd /mnt/xuruizhi/RNA_brain/human
vim 1.list
# PSM
SRR21161730
SRR21161731
SRR21161738
SRR21161739
SRR21161881
SRR21161882
SRR21161914
SRR21161915
# VLPFC
SRR21161734
SRR21161735
SRR21161750
SRR21161751
SRR21161759
SRR21161760
SRR21161765
SRR21161766
SRR21161909
SRR21161910
SRR21161931
SRR21161932
# CERE
SRR21161742
SRR21161743
SRR21161767
SRR21161768
SRR21161780
SRR21161781
SRR21161961
SRR21161962


vim 4.list
# OFC
SRR21163196
SRR21163197
SRR21161771
SRR21161772
SRR21161782
SRR21161783
SRR21161887
SRR21161888
SRR21161941
SRR21161942
# NACC
SRR21161776
SRR21161777
SRR21161811
SRR21161812
SRR21161937
SRR21161938
# CN
SRR21161784
SRR21161785
SRR21161823
SRR21161824
SRR21161831
SRR21161832
SRR21161920
SRR21161921
SRR21161963
SRR21161964

vim 5.list
# AMY
SRR21161763
SRR21161764
SRR21161790
SRR21161791
SRR21161795
SRR21161796
SRR21161803
SRR21161804
SRR21161905
SRR21161906
SRR21161929
SRR21161930
# DLPFC
SRR21161801
SRR21161802
SRR21161821
SRR21161822
SRR21161829
SRR21161830
SRR21161912
SRR21161913
SRR21161943
SRR21161944
SRR21161978
SRR21161979
# HIPP
SRR21161825
SRR21161826
SRR21161843
SRR21161844
SRR21161853
SRR21161854
SRR21161893
SRR21161894
SRR21161923
SRR21161924



vim 2.list 先不处理，因为此脑区不是很重要
SRR21161716
SRR21161717
SRR21161736
SRR21161737
SRR21161744
SRR21161745
SRR21161753
SRR21161754
SRR21161918
SRR21161919
SRR21161957
SRR21161958
SRR21161718
SRR21161719
SRR21161724
SRR21161725
SRR21161732
SRR21161733
SRR21161873
SRR21161874
SRR21161955
SRR21161956


vim 3.list
SRR21161837
SRR21161838
SRR21161846
SRR21161847
SRR21161933
SRR21161934
SRR21161746
SRR21161747
SRR21161773
SRR21161774
SRR21161899
SRR21161900
SRR21161951
SRR21161952
SRR21161968
SRR21161969

mkdir -p ./sra
cd ./sra
cp /mnt/d/perl/perl_scripts/download_srr.pl ./
cp ../*.list ./ 
# 批量下载
cat 1.list | parallel -k -j 6 "
  echo {} >> ./download.log
  perl download_srr.pl --output-dir . --srr {} >> ./download.log 2>&1
"
cat download.log | grep " downloaded successfully"
```

```bash
# 在本地转换格式
mkdir -p /mnt/xuruizhi/RNA_brain/human/sequence
cd /mnt/xuruizhi/RNA_brain/human/sequence
cp ../human.list ./

fqdir=/mnt/xuruizhi/RNA_brain/human/sequence
cd $fqdir
vim sra2fq.sh
#!/usr/bin bash
#sra2fq
fastq-dump --gzip --split-3 -O $fqdir ../sra/{}/{}.sra

cat human.list | while read id
do
  sed "s/{}/${id}/g" sra2fq.sh > ${id}_sra2fq.sh
done

cat human.list | parallel --no-run-if-empty --linebuffer -k -j 8 " 
  bash {}_sra2fq.sh"

```
2. 基因组数据



## 12.2 比对前质控

* 本地循环
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/fastqc
mkdir -p /mnt/xuruizhi/RNA_brain/human/trim
mkdir -p /mnt/xuruizhi/RNA_brain/human/fastqc_again

# 双端
ls *_1.fastq.gz | sed 's/_1.fastq.gz//g'  >pair.list
cd /mnt/xuruizhi/RNA_brain/human/sequence
# 循环脚本
vim human_pair.sh
#!/usr/bin bash
# This script is for pari-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/human/sequence /mnt/xuruizhi/ATAC_brain/human/sra/{}/{}.sra

# fastqc
 fastqc -t 6 -o ../fastqc {}_1.fastq.gz
 fastqc -t 6 -o ../fastqc {}_2.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o ../trim  {}_1.fastq.gz  {}_2.fastq.gz

# fatsqc_again
fastqc -t 6 -o ../fastqc_again ../trim/{}_1_val_1.fq.gz
fastqc -t 6 -o ../fastqc_again ../trim/{}_2_val_2.fq.gz




cat pair.list | while read id
do
  sed "s/{}/${id}/g" human_pair.sh > ${id}_pair_trim.sh
done
cat pair.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_pair_trim.sh >> ../trim/trim_fastqc.log 2>&1"


# 单端
cd /mnt/xuruizhi/RNA_brain/human/sequence
cat >single.list <<EOF
SRR13443447
SRR13443448
SRR13443449
EOF

# 循环脚本
vim human_single.sh
#!/usr/bin bash
# This script is for single-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/human/sequence /mnt/xuruizhi/ATAC_brain/human/sra/{}/{}.sra

# fastqc
fastqc -t 6 -o ../fastqc {}.fastq.gz

# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 -o ../trim {}.fastq.gz

# fatsqc_again
fastqc -t 6 -o ../fastqc_again ../trim/{}_trimmed.fq.gz




cat single.list | while read id
do
  sed "s/{}/${id}/g" human_single.sh > ${id}_single_trim.sh
  bash ${id}_single_trim.sh >> ../trim/trim_fastqc.log 2>&1
done


multiqc .

# 因为CG含量质控不合格，删除了
# SRR14494965 2质量不合格，但是1还可以，保留
SRR14494966 
SRR14494968 
SRR14494970 
SRR14494971 
```
```bash
# 传到超算
# sequence后来的补充内容没有上传到超算
rsync -av /mnt/xuruizhi/ATAC_brain/human/trim \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -av /mnt/xuruizhi/ATAC_brain/human/fastqc \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
rsync -av /mnt/xuruizhi/ATAC_brain/human/fastqc_again \
wangq@202.119.37.251:/scratch/wangq/xrz/ATAC_brain/human/
```

## 12.3 比对

1. 参考基因组  
```bash
# 小鼠 mm10 
mkdir -p /mnt/xuruizhi/RNA_brain/human/genome
mkdir -p /mnt/d/RNA_brain/human/genome
cd /mnt/d/RNA_brain/human/genome
wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
tar xvzf mm10_genome.tar.gz
rm mm10_genome.tar.gz
cp -r /mnt/d/RNA_brain/human/genome/* /mnt/xuruizhi/RNA_brain/human/genome

# 通过beyond compare传输超算
mkdir -p /scratch/wangq/xrz/RNA_brain/human/


cd /mnt/xuruizhi/RNA_brain/human/trim
vim pair.list
SRR14494965 
SRR14494974 
SRR3595255 
SRR3595256 
SRR3595258 
SRR11179785 
SRR11179786 
SRR11179787


vim single.list
SRR13443447 
SRR13443448 
SRR13443449 
```

2. 比对
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/align
# 循环 
cd /mnt/xuruizhi/RNA_brain/human/trim

# 双端
hisat2 -p 6 -t -x ../genome/mm10/genome -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz -S ../align/{}.sam 2>../align/{}.log 2>&1

# 单端
hisat2 -p 6 -t -x ../genome/mm10/genome -U {}_trimmed.fq.gz -S ../align/{}.sam 2>../align/{}.log 2>&1
```

3. sort_transfertobam_index  
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/sort_bam
# 双端与单端一致
samtools sort -@ 6 ../align/{}.sam > ../sort_bam/{}.sort.bam
samtools index -@ 6 ../sort_bam/{}.sort.bam
samtools flagstat  -@ 6 ../sort_bam/{}.sort.bam > ../sort_bam/{}.raw.stat
```
4. 大批量处理
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/align
mkdir -p /mnt/xuruizhi/RNA_brain/human/sort_bam
cd /mnt/xuruizhi/RNA_brain/human/trim

vim human_pair.sh
#!/usr/bin bash
# This script is for pair-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/human/sequence /mnt/xuruizhi/ATAC_brain/human/sra/{}/{}.sra

# fastqc
#  fastqc -t 6 -o ../fastqc {}_1.fastq.gz
#  fastqc -t 6 -o ../fastqc {}_2.fastq.gz

# # trim
# trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o ../trim  {}_1.fastq.gz  {}_2.fastq.gz

# # fatsqc_again
# fastqc -t 6 -o ../fastqc_again ../trim/{}_1_val_1.fq.gz
# fastqc -t 6 -o ../fastqc_again ../trim/{}_2_val_2.fq.gz

# align and sort
hisat2 -p 6 -t -x ../genome/mm10/genome -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz -S ../align/{}.sam 
samtools sort -@ 6 ../align/{}.sam > ../sort_bam/{}.sort.bam
samtools index -@ 6 ../sort_bam/{}.sort.bam
samtools flagstat  -@ 6 ../sort_bam/{}.sort.bam > ../sort_bam/{}.raw.stat


vim human_single.sh
#!/usr/bin bash
# This script is for single-end sequence.

# sra2fq.sh
# fastq-dump --gzip --split-3 -O /mnt/xuruizhi/ATAC_brain/human/sequence /mnt/xuruizhi/ATAC_brain/human/sra/{}/{}.sra

# fastqc
# fastqc -t 6 -o ../fastqc {}.fastq.gz

# # trim
# trim_galore --phred33 --length 35 -e 0.1 --stringency 3 -o ../trim {}.fastq.gz

# # fatsqc_again
# fastqc -t 6 -o ../fastqc_again ../trim/{}_trimmed.fq.gz

# align
hisat2 -p 6 -t -x ../genome/mm10/genome -U {}_trimmed.fq.gz -S ../align/{}.sam 
samtools sort -@ 6 ../align/{}.sam > ../sort_bam/{}.sort.bam
samtools index -@ 6 ../sort_bam/{}.sort.bam
samtools flagstat  -@ 6 ../sort_bam/{}.sort.bam > ../sort_bam/{}.raw.stat
```
```bash
cd /mnt/xuruizhi/RNA_brain/human/trim

# 双端 
cat pair.list  | while read id
do 
  sed "s/{}/${id}/g" human_pair.sh > ${id}_pair_align.sh
  bash  ${id}_pair_align.sh >> ../align/align_pair.log 2>&1
done
# 单端
cat single.list  | while read id
do 
  sed "s/{}/${id}/g" human_single.sh > ${id}_single_align.sh
  bash ${id}_single_align.sh >> ../align/align_single.log 2>&1
done
```

## 12.4 表达量统计
1. 下载注释文件
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/annotation
cd /mnt/xuruizhi/RNA_brain/human/annotation
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ensGene.gtf.gz
gzip -dc mm10.ensGene.gtf.gz > mm10.ensGene.gtf
```

2. 统计
```bash
cd /mnt/xuruizhi/RNA_brain/human/sort_bam
mkdir -p /mnt/xuruizhi/RNA_brain/human/HTseq

parallel -j 8 "
    htseq-count -s no -f bam {1}.sort.bam ../annotation/mm10.ensGene.gtf > ../HTseq/{1}.count  2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

mkdir -p /mnt/d/RNA_brain/human/HTseq
cp /mnt/xuruizhi/RNA_brain/human/HTseq/* /mnt/d/RNA_brain/human/HTseq
```
```r
rm(list=ls())
setwd("D:/RNA_brain/human/HTseq")

# 得到文件样本编号
files <- list.files(".", "*.count")
f_lists <- list()
for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[prefix]] = i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  data[[count]] <- a
}

# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]],by="gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```
## 12.5 差异表达分析

1. 数据筛选
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/Deseq2
mkdir -p /mnt/d/RNA_brain/human/Deseq2

# 与RNA-seq对应：HIPP,DG,PFC,SEN
vim human.list 
SRR14494965 PFC.
SRR14494974 PFC.
SRR3595255 DG.
SRR3595256 DG.
SRR3595258 DG.
SRR11179785 HIPP
SRR11179786 HIPP.
SRR11179787 HIPP.
SRR13443447 PUT. 
SRR13443448 MOB.
SRR13443449 SEN.

# 因为考虑到后续Deseq2只能对有重复的样本进行分析，因此还是先使用有重复的三个脑区
cd /mnt/d/RNA_brain/human/Deseq2
vim human.list 
SRR14494965
SRR14494974
SRR3595255
SRR3595256
SRR3595258
SRR11179785
SRR11179786
SRR11179787

while read -r i
do
  cp ../HTseq/${i}.count ./
done < human.list
```
```r
rm(list=ls())
setwd("D:/RNA_brain/human/Deseq2")

# 得到文件样本编号
files <- list.files(".", "*.count")
f_lists <- list()
for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[prefix]] = i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  data[[count]] <- a
}

# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]],by="gene_id")
}
write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```
2. 差异分析，以DG为control
* coldata
```bash
cd /mnt/d/RNA_brain/human/Deseq2 
vim coldata.csv
"ids","state","condition","treatment"
"SRR14494965","WT","PFC","treatment"
"SRR14494974","WT","PFC","treatment"
"SRR3595255","WT","DG","treatment"
"SRR3595256","WT","DG","treatment"
"SRR3595258","WT","DG","treatment"
"SRR11179785","WT","HIPP","treatment"
"SRR11179786","WT","HIPP","treatment"
"SRR11179787","WT","HIPP","treatment"
```
```r
BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# 导入countdata文件
dataframe <- read.csv("merge.csv", header=TRUE, row.names = 1)
countdata <- dataframe[-(1:5),]
countdata <- countdata[rowSums(countdata) > 0,]
head(countdata)
# 导入coltdata文件
coldata <- read.table("coldata.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ condition)
dds
# 归一化
rld <- rlog(dds, blind=FALSE)

# PCA分析 
# intgroup分组
pcaData <- plotPCA(rld, intgroup=c("condition"),returnData = T) 
pcaData <- pcaData[order(pcaData$condition,decreasing=F),]
table(pcaData$condition)
# PCA1
plot(pcaData[,1:2],pch = 19,col= c(rep("red",3),rep("green",3),rep("blue",2)))+ text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=0.75, font = 1)+legend(-30,-5,inset = .02,pt.cex= 1.5,legend = c("DG","HIPP","PFC"), col = c( "red","green","blue"),pch = 19, cex=0.75,bty="n")
# PCA2
plotPCA(rld, intgroup="condition") + ylim(-30, 30)+text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=0.5, font = 1)

# 聚类热图
library("RColorBrewer")
# 得到数据对象中基因的计数的转化值
gene_data_transform <- assay(rld)
# 使用dist方法求样本之间的距离
sampleDists <- dist(t(gene_data_transform))
# 转化为矩阵用于后续pheatmap()方法的输入
sampleDistMatrix <- as.matrix(sampleDists)
# 将矩阵的名称进行修改
rownames(sampleDistMatrix) <- rld$condition
colnames(sampleDistMatrix) <- rld$condition
# 设置色盘
colors <- colorRampPalette(rev(brewer.pal(8, "Blues")) )(255)
# 绘制热图与聚类
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


dds$condition <- factor(as.vector(dds$condition), levels = c("DG","HIPP","PFC")) 
dds$condition
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"            "condition_HIPP_vs_DG"
# [3] "condition_PFC_vs_DG"
```



* HIPP vs DG 
```r
result <- results(dds, name="condition_HIPP_vs_DG", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4287, 16%
# LFC < 0 (down)     : 5070, 19%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 6192, 23%
# (mean count < 2)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 11049  9357
write.csv(result_order, file="HIPP_vs_DG.csv", quote = F)

# 火山图
voldata <- read.csv("HIPP_vs_DG.csv", header=TRUE, row.names = 1)
voldata$plog <- (-log10(voldata$padj))
voldata$compare <- ifelse(abs(voldata$log2FoldChange) >1 & voldata$padj < 0.05, ifelse(voldata$log2FoldChange >1, "up","down"), "no")
ggplot(voldata, aes(x=log2FoldChange, y=plog,color= compare)) +
  geom_point(alpha=.5) +
  theme(panel.grid.major = element_blank(),
        axis.ticks.x = element_line(size=1),
        axis.text.x = element_text(angle=30, hjust=1, vjust=1),
        axis.title.x=element_text(face="italic", colour="darkred", size=14), # 字体
        axis.line = element_line(color="black"),
        plot.title = element_text(colour="red", size=8, face="bold")) +
  ylab("-log10(padj)") + 
  xlab("log2FC") +
  ggtitle("differencial genes") +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = -1:1) + 
  scale_color_manual(values = c("blue", "grey", "red"))
# MA图  
png("MA_HIPP_DG.png")
plotMA(result_order,ylim=c(-12,12))
dev.off()

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
# 查看数据框的大小
dim(diff_gene)   #3585    6
write.csv(diff_gene, file="diff_HIPP_DG.csv", quote = F)
```




* PFC_vs_DG
```r
result <- results(dds, name="condition_PFC_vs_DG", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7127, 27%
# LFC < 0 (down)     : 6635, 25%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 6192, 23%
# (mean count < 2)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  6644 13762

write.csv(result_order, file="PFC_vs_DG.csv", quote = F)
png("MA_PFC_DG.png")
plotMA(result_order,ylim=c(-12,12))
dev.off()

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #10617    6
write.csv(diff_gene, file="diff_PFC_DG.csv", quote = F)
```


2. 差异分析，以HIPP为control
```r
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# 其他和前文一致
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ condition)
dds$condition <- factor(as.vector(dds$condition), levels = c("HIPP","DG","PFC")) 
dds$condition
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"             "condition_DG_vs_HIPP"  "condition_PFC_vs_HIPP"
```



* DG_vs_HIPP
```r
result <- results(dds, name="condition_DG_vs_HIPP", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order) # 和前文HIPP vs DG刚好反过来
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5070, 19%
# LFC < 0 (down)     : 4287, 16%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 6192, 23%
# (mean count < 2)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 11049  9357
write.csv(result_order, file="DG_vs_HIPP.csv", quote = F)

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #3585    6
write.csv(diff_gene, file="diff_DG_HIPP.csv", quote = F)
```



* PFC_vs_HIPP
```r
result <- results(dds, name="condition_PFC_vs_HIPP", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 6782, 25%
# LFC < 0 (down)     : 6938, 26%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 5676, 21%
# (mean count < 1)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  7202 13720 

write.csv(result_order, file="PFC_vs_HIPP.csv", quote = F)
png("MA_PFC_HIPP.png")
plotMA(result_order,ylim=c(-12,12))
dev.off()

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #10400    6
write.csv(diff_gene, file="diff_PFC_HIPP.csv", quote = F)
```

3. 差异分析，以PFC为control
```r
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# 其他和前文一致
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ condition)
dds$condition <- factor(as.vector(dds$condition), levels = c("PFC","HIPP","DG")) 
dds$condition
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"             "condition_HIPP_vs_PFC" "condition_DG_vs_PFC" 
```



* HIPP_vs_PFC
```r
result <- results(dds, name="condition_HIPP_vs_PFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order) 
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 6938, 26%
# LFC < 0 (down)     : 6782, 25%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 5676, 21%
# (mean count < 1)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  7202 13720 
write.csv(result_order, file="HIPP_vs_PFC.csv", quote = F)

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #10400    6
write.csv(diff_gene, file="diff_HIPP_PFC.csv", quote = F)
```



* DG_vs_PFC
```r
result <- results(dds, name="condition_DG_vs_PFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 26613 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 6635, 25%
# LFC < 0 (down)     : 7127, 27%
# outliers [1]       : 15, 0.056%
# low counts [2]     : 6192, 23%
# (mean count < 2)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  6644 13762 

write.csv(result_order, file="DG_vs_PFC.csv", quote = F)

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 10617    6
write.csv(diff_gene, file="diff_DG_PFC.csv", quote = F)
```


4. 脑区差异基因

找到某一脑区对另外两个脑区来说，都差异表达的基因。

```bash
cd /mnt/d/RNA_brain/human/Deseq2 
# diff_HIPP_DG.csv diff_HIPP_PFC.csv

# 区分up/down
for i in *.csv
do
sed 's/,/\t/g' "$i" | tail -n +2  > "${i%.csv}.txt"
done
dos2unix *.txt

for i in diff*.txt 
do
  echo " ==> $i <== " 
  tsv-filter --is-numeric 3 --gt 3:0 $i > ${i%%.*}_up.txt
  tsv-filter --is-numeric 3 --lt 3:0 $i > ${i%%.*}_down.txt
done

wc -l *_up.txt
   1114 diff_DG_HIPP_up.txt
   5077 diff_DG_PFC_up.txt
   2471 diff_HIPP_DG_up.txt
   5348 diff_HIPP_PFC_up.txt
   5540 diff_PFC_DG_up.txt
   5052 diff_PFC_HIPP_up.txt

wc -l *_down.txt
   2471 diff_DG_HIPP_down.txt
   5540 diff_DG_PFC_down.txt
   1114 diff_HIPP_DG_down.txt
   5052 diff_HIPP_PFC_down.txt
   5077 diff_PFC_DG_down.txt
   5348 diff_PFC_HIPP_down.txt
```

* HIPP对于其他两个脑区的差异基因
```bash
# up
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_HIPP_DG_up.txt diff_HIPP_PFC_up.txt > diff_HIPP_up.txt # 932
# down
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_HIPP_DG_down.txt diff_HIPP_PFC_down.txt > diff_HIPP_down.txt # 397

cat diff_HIPP_up.txt diff_HIPP_down.txt > diff_HIPP.txt
```
```r
library(biomaRt)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)


  setwd("D:/RNA_brain/human/Deseq2")
  region <- c("HIPP")
  data <- read.table(paste0("diff_",region,".txt"), header=FALSE)
  
  ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Mm.eg.db")
    return(a)
  }
  region_ensembl_id_transform <- ensembl_id_transform(data$V1)
  write.csv(ensembl_id_transform(data$V1), file =  paste0("diff_",region,"_ensemblID.tsv"))

  mart <- useDataset("mmusculus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
  region_biomart_ensembl_id_transform <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
    filters = "ensembl_gene_id",
    values = data$V1,
    mart = mart
  )
  write.csv(region_biomart_ensembl_id_transform, file =  paste0("diff_",region,"_biomartID.tsv"))

  # GO analysis and barplot
  region_biomart <- enrichGO(
    gene = region_biomart_ensembl_id_transform$entrezgene_id, 
    keyType = "ENTREZID",
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  pdf(file = paste0(region, "_biomart.pdf"))
  barplot(region_biomart, showCategory = 40, font.size = 6, title = paste("The GO BP enrichment analysis", sep = ""))
  dev.off()

  region_transform <- enrichGO(
    gene = region_ensembl_id_transform$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  pdf(file = paste0(region, "_transform.pdf"))
  barplot(region_transform, showCategory = 40, font.size = 6, title = paste("The GO BP enrichment analysis", sep = ""))
  dev.off()

  region_kegg <- enrichKEGG(
    gene = region_ensembl_id_transform$ENTREZID,
    organism = 'mmu',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(file = paste0(region, "_kegg.pdf"),width = 80, height = 120)
  barplot(region_kegg, showCategory = 20, font.size = 120,title = "KEGG Pathway Enrichment Analysis")
  dev.off()
```

* DG 对于其他两个脑区的差异基因
```bash
# up
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_DG_HIPP_up.txt diff_DG_PFC_up.txt > diff_DG_up.txt # 597
# down
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_DG_HIPP_down.txt diff_DG_PFC_down.txt > diff_DG_down.txt # 1363

cat diff_DG_up.txt diff_DG_down.txt > diff_DG.txt
```
```r
  region <- c("DG")
  data <- read.table(paste0("diff_",region,".txt"), header=FALSE)
```

* PFC 对于其他两个脑区的差异基因
```bash
# up
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_PFC_DG_up.txt diff_PFC_HIPP_up.txt > diff_PFC_up.txt # 4070
# down
awk 'NR==FNR {a[$1]=1; next} a[$1]' diff_PFC_DG_down.txt diff_PFC_HIPP_down.txt > diff_PFC_down.txt # 4154

cat diff_PFC_up.txt diff_PFC_down.txt > diff_PFC.txt
```
```r
  region <- c("PFC")
  data <- read.table(paste0("diff_",region,".txt"), header=FALSE)
```
```bash
cp /mnt/d/RNA_brain/human/Deseq2/* /mnt/xuruizhi/RNA_brain/human/Deseq2
```
