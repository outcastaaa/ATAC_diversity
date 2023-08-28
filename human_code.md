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

# 5. 合并neuron和non-neuron 还未全部完成

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
# 1.list
samtools merge -@ 4 -o PMC_rep1.bam SRR21163180.final.bam SRR21163181.final.bam
samtools merge -@ 4 -o PMC_rep2.bam SRR21163186.final.bam SRR21163187.final.bam
samtools merge -@ 4 -o PMC_rep3.bam SRR21163293.final.bam SRR21163294.final.bam
samtools merge -@ 4 -o VLPFC_rep1.bam SRR21163184.final.bam SRR21163185.final.bam
samtools merge -@ 4 -o VLPFC_rep2.bam SRR21163207.final.bam SRR21163208.final.bam
samtools merge -@ 4 -o VLPFC_rep3.bam SRR21163320.final.bam SRR21163321.final.bam
samtools merge -@ 4 -o VLPFC_rep4.bam SRR21163337.final.bam SRR21163338.final.bam
samtools merge -@ 4 -o CRBLM_rep1.bam SRR21163190.final.bam SRR21163191.final.bam
samtools merge -@ 4 -o CRBLM_rep2.bam SRR21163209.final.bam SRR21163210.final.bam
samtools merge -@ 4 -o CRBLM_rep3.bam SRR21163216.final.bam SRR21163217.final.bam
samtools merge -@ 4 -o CRBLM_rep4.bam SRR21163365.final.bam SRR21163366.final.bam
samtools merge -@ 4 -o HIPP_rep1.bam SRR21163267.final.bam SRR21163268.final.bam
samtools merge -@ 4 -o HIPP_rep2.bam SRR21163304.final.bam SRR21163305.final.bam
ls *rep*.bam | parallel --no-run-if-empty --linebuffer -k -j 2 '
samtools index -@ 4 {}
samtools flagstat -@ 4 {} > {}.stat'
```

# 6. 对于合并和未合并都Call peaks 

```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/bed
mkdir -p /mnt/xuruizhi/ATAC_brain/human/Tn5_shift
mkdir -p /mnt/xuruizhi/ATAC_brain/human/peaks

cd /mnt/xuruizhi/ATAC_brain/human/final
# 1_all.list 3_all.list HIPP_all.list对应详细SRR名称
vim 1.list
PMC_rep1
PMC_rep2
PMC_rep3
VLPFC_rep1
VLPFC_rep2
VLPFC_rep3
VLPFC_rep4
CRBLM_rep1
CRBLM_rep2
CRBLM_rep3
CRBLM_rep4

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

vim HIPP.list
HIPP_rep1
HIPP_rep2

cat 1.list | while read id
do
  samtools flagstat -@ 6 ${id}.bam > ${id}.stat
done



vim human.sh
# shift bed
# cd /mnt/xuruizhi/ATAC_brain/human/final
bedtools bamtobed -i {}.final.bam | awk -v OFS="\t" '{
    if ($6 == "+") {
        print $1, $2+4, $3+4;
    } else if ($6 == "-") {
        print $1, $2-5, $3-5;
    }
}' > ../Tn5_shift/{}.Tn5.bed
### 一定记得换物种要换-g参数
macs2 callpeak  -g hs --shift -75 --extsize 150 --nomodel \
--nolambda --keep-dup all -n {} -t ../Tn5_shift/{}.Tn5.bed --outdir ../peaks/


cat HIPP_all.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ../peaks/${id}_callpeaks.sh
done
cat HIPP_all.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash ../peaks/{}_callpeaks.sh >> ../peaks/peaks.log 2>&1"

# 因为前面做的都出错，因此直接全部重新call peak
cd ../Tn5_shift
ls *.Tn5.bed | parallel --no-run-if-empty --linebuffer -k -j 6 " 
macs2 callpeak  -g hs --shift -75 --extsize 150 --nomodel --nolambda --keep-dup all -n {} -t ./{} --outdir ../peaks/ >> ../peaks/peaks_new.log 2>&1 "
```

## 后续内容先只处理1.list（PMC,VLPFC,CRBLM）+ 3.list(HAB,PSC)+ HIPP.lsit

# 7. Quality check 
1. fragment length distribution —— 此步骤不受peak影响
```bash
# 在Linux中画图  
mkdir -p /mnt/xuruizhi/ATAC_brain/human/frag_length
cd /mnt/xuruizhi/ATAC_brain/human/final

cat 3_all.list | parallel -k -j 6 "
  echo {} 
  java -jar /mnt/d/biosoft/picard/picard.jar CollectInsertSizeMetrics \
  -I {}.final.bam \
  -O ../frag_length/{}.insert_size_metrics.txt \
  -H ../frag_length/{}.insert_size_histogram.pdf"

cat 3.list | parallel -k -j 6 "
  echo {} 
  java -jar /mnt/d/biosoft/picard/picard.jar CollectInsertSizeMetrics \
  -I {}.bam \
  -O ../frag_length/{}.insert_size_metrics.txt \
  -H ../frag_length/{}.insert_size_histogram.pdf"
```

3. 统计平均长度与长度分布情况
* bash统计平均值
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/length
cd /mnt/xuruizhi/ATAC_brain/human/peaks

for i in *_peaks.narrowPeak
do
  echo $i
  cat $i | awk '{print $3-$2}' | awk '{sum += $1} END {print avg = sum/NR}' > ../length/${i%%.*}_ave.txt
  awk '{print $3- $2}' $i > ${i%%.*}_length.txt
done

rsync -avP /mnt/xuruizhi/ATAC_brain/human/peaks /mnt/d/ATAC_brain/human
```
* R中统计与画图
！ 将R储存在 D:/ATAC_brain/human文件夹中
```r
setwd("D:/ATAC_brain/human/peaks")
file_list <- list.files(pattern = "\\d+_peaks_length\\.txt")
summary_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE)
  summary_list[[file]] <- summary(data)
}

output_file <- "summary_output.txt" 

output_lines <- character(length(file_list))
for (i in seq_along(file_list)) {
  output_lines[i] <- paste("Summary for", file_list[i], ":\n", capture.output(print(summary_list[[file_list[i]]])), collapse = "\n")
}

writeLines(output_lines, con = output_file)


# 画peak length 直方图
summary_list <- list()
file_list <- list.files(path = "./", pattern = "\\d+_peaks_length\\.txt", full.names = TRUE)

for (file in file_list) {
  data <- read.table(file)
  dim_data <- dim(data)
  summary_list[[file]] <- dim_data

  # 生成直方图并保存为png文件
  png(paste0(file, "_hist.png"))
  # pdf(paste0(file, "_hist.pdf"))
  hist(abs(as.numeric(data[, 1])), breaks = 500, xlab = "peak length (bp)", ylab = "Frequency", main = file)
  dev.off()
}

for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}
```
4. TSS enrichment在下一节展示，很重要的判断质量的标准


# 8. Visualization    
1.  filterbam2Bw  —— 此步骤不受peak影响
```bash 
mkdir -p  /mnt/xuruizhi/ATAC_brain/human/bw
cd /mnt/xuruizhi/ATAC_brain/human/final
conda activate py3.8
cat 3.list | while read id; 
do 
  bamCoverage -p 6  -b $id.bam \
  -o ../bw/${id}.bw \
  --binSize 20 \
  --smoothLength 60 \
  --normalizeUsing RPKM \
  --centerReads 
  >> ../bw/bamCoverage1.log 2>&1
done

cat 3_all.list | parallel --no-run-if-empty --linebuffer -k -j 4 '
  bamCoverage -p 2  -b {}.final.bam \
  -o ../bw/{}.bw \
  --binSize 20 \
  --smoothLength 60 \
  --normalizeUsing RPKM \
  --centerReads 
  >> ../bw/bamCoverage1all.log 2>&1 '
```

2. TSS enrichment   —— 此步骤不受peak影响

下载TSS注释文件： [下载地址](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1682652860_KcLBsnsbOiOTzM1hWST65rw7RhuK&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=bed&hgta_outFileName=hg38%C3%83%C6%92%C3%82%C2%A2%C3%83%C2%A2%C3%A2%E2%82%AC%C5%A1%C3%82%C2%AC%C3%83%C2%A2%C3%A2%E2%80%9A%C2%AC%C3%82%C2%9D%C3%83%C6%92%C3%82%C2%A2%C3%83%C2%A2%C3%A2%E2%82%AC%C5%A1%C3%82%C2%AC%C3%83%C2%A2%C3%A2%E2%80%9A%C2%AC%C3%82%C2%9Ducsc.bed)   

```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/TSS
cd /mnt/d/ATAC_brain/human/TSS
gzip -dc hg38_ucsc.bed.gz > hg38.TSS.bed
cp ./hg38.TSS.bed  /mnt/xuruizhi/ATAC_brain/human/TSS


cd /mnt/xuruizhi/ATAC_brain/human/bw
cat ../final/3_all.list | while read id; 
do 
  computeMatrix reference-point --referencePoint TSS -p 6 \
    -b 1000  -a 1000 \
    -R ../TSS/hg38.tss.bed \
    -S ${id}.bw \
    --skipZeros \
    -o ../TSS/${id}_matrix.gz \
    --outFileSortedRegions ../TSS/${id}_regions.bed
done


# profile plot
cd ../TSS/
cat ../final/1_all.list | while read id
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
cat ../final/1_all.list | while read id
do 
  plotHeatmap -m ${id}_matrix.gz \
  -out ${id}_heatmap.png \
  --colorMap RdBu \
  --whatToShow 'heatmap and colorbar' \
  --zMin -8 --zMax 8  
done


# heatmap and profile plot
cat ../final/1_all.list | while read id
do 
  plotHeatmap -m ${id}_matrix.gz \
    -out ${id}_all.png \
    --colorMap RdBu \
    --zMin -12 --zMax 12
done


# 画 `gene body` 区，使用 `scale-regions`  
cd ../bw
mkdir -p ../genebody
ls *.bw | while read id; 
do
computeMatrix scale-regions -p 6 \
    -b 10000  -a 10000 \
    -R ../TSS/hg38.tss.bed  \
    -S ${id} \
    --skipZeros \
    -o ../genebody/${id%%.*}_matrix.gz 
done

ls ../genebody/*.gz | while read id
do
  plotHeatmap -m ${id} -out ${id%%.*}_heatmap.png 
done
```

# 9. 寻找rep间consensus peak —— IDR
1. 对单独样本进行整理
* 见"human数据质量情况.xlsx文件"  
  
此处IDR分析为：单独神经元；单独非神经元；合并后bulk；分别进行IDR

```bash
# Sort peak by -log10(p-value)
mkdir -p /mnt/xuruizhi/ATAC_brain/human/IDR
cd /mnt/xuruizhi/ATAC_brain/human/peaks
cat ../final/1.list | parallel -j 6 -k "
sort -k8,8nr {}_peaks.narrowPeak > ../IDR/{}.narrowPeak 
"

cd ../IDR
cp /mnt/xuruizhi/ATAC_brain/mouse/IDR/idr.sh ./
```

2. PMC
① 神经元
```bash
# SRR21163181,SRR21163187,SRR21163294
sudo chmod +x idr.sh
conda activate py3.8
./idr.sh SRR21163181.narrowPeak SRR21163187.narrowPeak SRR21163294.narrowPeak .
mv SRR21163181_SRR21163187.txt PMC_neu1.txt 
mv SRR21163181_SRR21163294.txt PMC_neu2.txt 
mv SRR21163187_SRR21163294.txt PMC_neu3.txt 

for i in PMC_neu*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in PMC_neu*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat PMC_neu*_common.bed | tsv-summarize -g 1 --count
cat PMC_neu*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"  > PMC_neu_pool.bed
sort -k1,1 -k2,2n PMC_neu_pool.bed | bedtools merge -i stdin -d 50 > PMC_neu_pool_merge.bed 

wc -l PMC_neu*.bed
  # 11728 PMC_neu1_common.bed
  #  6858 PMC_neu2_common.bed
  #  5856 PMC_neu3_common.bed
  # 24398 PMC_neu_pool.bed
  # 14971 PMC_neu_pool_merge.bed
```
② 非神经元
```bash
# SRR21163180,SRR21163186,SRR21163293

./idr.sh SRR21163180.narrowPeak SRR21163186.narrowPeak SRR21163293.narrowPeak .
mv SRR21163180_SRR21163186.txt PMC_non1.txt 
mv SRR21163180_SRR21163293.txt PMC_non2.txt 
mv SRR21163186_SRR21163293.txt PMC_non3.txt 

for i in PMC_non*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in PMC_non*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat PMC_non*_common.bed | tsv-summarize -g 1 --count
cat PMC_non*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"  > PMC_non_pool.bed
sort -k1,1 -k2,2n PMC_non_pool.bed | bedtools merge -i stdin -d 50 > PMC_non_pool_merge.bed 
wc -l  PMC_non*.bed
  # 15946 PMC_non1_common.bed
  # 13837 PMC_non2_common.bed
  # 15940 PMC_non3_common.bed
  # 45464 PMC_non_pool.bed
  # 23898 PMC_non_pool_merge.bed
```
③ 合并
```bash
# PMC_rep1.narrowPeak,PMC_rep2.narrowPeak,PMC_rep3.narrowPeak
./idr.sh PMC_rep1.narrowPeak PMC_rep2.narrowPeak PMC_rep3.narrowPeak .
mv PMC_rep1_PMC_rep2.txt PMC_all1.txt 
mv PMC_rep1_PMC_rep3.txt PMC_all2.txt 
mv PMC_rep2_PMC_rep3.txt PMC_all3.txt 

for i in PMC_all*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in PMC_all*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat PMC_all*_common.bed | tsv-summarize -g 1 --count
cat PMC_all*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"  > PMC_all_pool.bed
sort -k1,1 -k2,2n PMC_all_pool.bed | bedtools merge -i stdin -d 50 > PMC_all_pool_merge.bed 
wc -l PMC_all*.bed
  # 22589 PMC_all1_common.bed
  # 19224 PMC_all2_common.bed
  # 18608 PMC_all3_common.bed
  # 60195 PMC_all_pool.bed
  # 32487 PMC_all_pool_merge.bed
```

3. VLPFC

① 神经元
```bash
# SRR21163185,SRR21163208,SRR21163321,SRR21163338
sudo chmod +x idr.sh
conda activate py3.8
./idr.sh SRR21163185.narrowPeak SRR21163208.narrowPeak SRR21163321.narrowPeak SRR21163338.narrowPeak .
mv SRR21163185_SRR21163208.txt VLPFC_neu1.txt 
mv SRR21163185_SRR21163321.txt VLPFC_neu2.txt 
mv SRR21163185_SRR21163338.txt VLPFC_neu3.txt
mv SRR21163208_SRR21163321.txt VLPFC_neu4.txt
mv SRR21163208_SRR21163338.txt VLPFC_neu5.txt 
mv SRR21163321_SRR21163338.txt VLPFC_neu6.txt 

for i in VLPFC_neu*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in VLPFC_neu*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat VLPFC_neu*_common.bed | tsv-summarize -g 1 --count
cat VLPFC_neu*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"   > VLPFC_neu_pool.bed
sort -k1,1 -k2,2n VLPFC_neu_pool.bed | bedtools merge -i stdin -d 50 > VLPFC_neu_pool_merge.bed 
wc -l VLPFC_neu*.bed
  #  5649 VLPFC_neu1_common.bed
  #  4481 VLPFC_neu2_common.bed
  #  5092 VLPFC_neu3_common.bed
  #  7915 VLPFC_neu4_common.bed
  #  8603 VLPFC_neu5_common.bed
  #  9778 VLPFC_neu6_common.bed
  # 41322 VLPFC_neu_pool.bed
  # 16208 VLPFC_neu_pool_merge.bed

```
② 非神经元
```bash
# SRR21163184 SRR21163207  SRR21163320  SRR21163337

./idr.sh SRR21163184.narrowPeak SRR21163207.narrowPeak SRR21163320.narrowPeak SRR21163337.narrowPeak .
mv SRR21163184_SRR21163207.txt VLPFC_non1.txt 
mv SRR21163184_SRR21163320.txt VLPFC_non2.txt 
mv SRR21163184_SRR21163337.txt VLPFC_non3.txt 
mv SRR21163207_SRR21163320.txt VLPFC_non4.txt 
mv SRR21163207_SRR21163337.txt VLPFC_non5.txt 
mv SRR21163320_SRR21163337.txt VLPFC_non6.txt  


for i in VLPFC_non*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in VLPFC_non*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat VLPFC_non*_common.bed | tsv-summarize -g 1 --count
cat VLPFC_non*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"   > VLPFC_non_pool.bed
sort -k1,1 -k2,2n VLPFC_non_pool.bed | bedtools merge -i stdin -d 50 > VLPFC_non_pool_merge.bed 
wc -l VLPFC_non*.bed
  # 11863 VLPFC_non1_common.bed
  #  9754 VLPFC_non2_common.bed
  # 12114 VLPFC_non3_common.bed
  #  8783 VLPFC_non4_common.bed
  # 11400 VLPFC_non5_common.bed
  # 11981 VLPFC_non6_common.bed
  # 65667 VLPFC_non_pool.bed
  # 21944 VLPFC_non_pool_merge.bed
```
③ 合并
```bash
# VLPFC_rep1.narrowPeak VLPFC_rep2.narrowPeak VLPFC_rep3.narrowPeak VLPFC_rep4.narrowPeak
./idr.sh VLPFC_rep1.narrowPeak VLPFC_rep2.narrowPeak VLPFC_rep3.narrowPeak VLPFC_rep4.narrowPeak .
mv VLPFC_rep1_VLPFC_rep2.txt VLPFC_all1.txt 
mv VLPFC_rep1_VLPFC_rep3.txt VLPFC_all2.txt 
mv VLPFC_rep1_VLPFC_rep4.txt VLPFC_all3.txt 
mv VLPFC_rep2_VLPFC_rep3.txt VLPFC_all4.txt 
mv VLPFC_rep3_VLPFC_rep4.txt VLPFC_all5.txt 
mv VLPFC_rep2_VLPFC_rep4.txt VLPFC_all6.txt 


for i in VLPFC_all*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in VLPFC_all*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat VLPFC_all*_common.bed | tsv-summarize -g 1 --count
cat VLPFC_all*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"  > VLPFC_all_pool.bed
cat VLPFC_all_pool.bed | sort | tsv-summarize -g 1 --count 
sort -k1,1 -k2,2n VLPFC_all_pool.bed | bedtools merge -i stdin -d 50 > VLPFC_all_pool_merge.bed 
wc -l VLPFC_all*.bed
  # 13062 VLPFC_all1_common.bed
  # 13963 VLPFC_all2_common.bed
  # 17161 VLPFC_all3_common.bed
  # 17207 VLPFC_all4_common.bed
  # 20876 VLPFC_all5_common.bed
  # 18209 VLPFC_all6_common.bed
  # 99993 VLPFC_all_pool.bed
  # 34274 VLPFC_all_pool_merge.bed
```

4. CRBLM 同上

！发现SRR21163210和SRR21163209和其他样本重复性很差，因此删除

① 神经元
```bash
# SRR21163191.narrowPeak SRR21163210.narrowPeak SRR21163217.narrowPeak SRR21163366.narrowPeak
./idr.sh SRR21163191.narrowPeak SRR21163210.narrowPeak SRR21163217.narrowPeak SRR21163366.narrowPeak .
mv SRR21163191_SRR21163210.txt CRBLM_neu1.txt 
mv SRR21163191_SRR21163217.txt CRBLM_neu2.txt 
mv SRR21163191_SRR21163366.txt CRBLM_neu3.txt
mv SRR21163210_SRR21163217.txt CRBLM_neu4.txt 
mv SRR21163210_SRR21163366.txt CRBLM_neu5.txt 
mv SRR21163217_SRR21163366.txt CRBLM_neu6.txt


for i in CRBLM_neu*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in CRBLM_neu*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
cat CRBLM_neu*_common.bed | tsv-summarize -g 1 --count
cat CRBLM_neu*_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"   > CRBLM_neu_pool.bed
sort -k1,1 -k2,2n CRBLM_neu_pool.bed | bedtools merge -i stdin -d 50 > CRBLM_neu_pool_merge.bed 

wc -l CRBLM_neu*.bed
  #  1828 CRBLM_neu1_common.bed
  # 25872 CRBLM_neu2_common.bed
  # 17825 CRBLM_neu3_common.bed
  #  1926 CRBLM_neu4_common.bed
  #  3622 CRBLM_neu5_common.bed
  # 16828 CRBLM_neu6_common.bed
  # 67752 CRBLM_neu_pool.bed
  # 33225 CRBLM_neu_pool_merge.bed

# 删除210，保留CRBLM_neu2_common.bed，CRBLM_neu3_common.bed，CRBLM_neu6_common.bed
cat CRBLM_neu2_common.bed CRBLM_neu3_common.bed CRBLM_neu6_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"   > CRBLM_neu_pool.bed
sort -k1,1 -k2,2n CRBLM_neu_pool.bed | bedtools merge -i stdin -d 50 > CRBLM_neu_pool_merge.bed 
wc -l CRBLM_neu*.bed
  # 25872 CRBLM_neu2_common.bed
  # 17825 CRBLM_neu3_common.bed
  # 16828 CRBLM_neu6_common.bed
  # 60376 CRBLM_neu_pool.bed
  # 33130 CRBLM_neu_pool_merge.bed
```
② 非神经元
```bash
# SRR21163190 SRR21163209 SRR21163216 SRR21163365
./idr.sh SRR21163190.narrowPeak SRR21163209.narrowPeak SRR21163216.narrowPeak SRR21163365.narrowPeak .
mv SRR21163190_SRR21163209.txt CRBLM_non1.txt 
mv SRR21163190_SRR21163216.txt CRBLM_non2.txt 
mv SRR21163190_SRR21163365.txt CRBLM_non3.txt
mv SRR21163209_SRR21163216.txt CRBLM_non4.txt 
mv SRR21163209_SRR21163365.txt CRBLM_non5.txt 
mv SRR21163216_SRR21163365.txt CRBLM_non6.txt


for i in CRBLM_non*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in CRBLM_non*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done
wc -l CRBLM_non*.bed
  #  2142 CRBLM_non1_common.bed
  # 41473 CRBLM_non2_common.bed
  # 34733 CRBLM_non3_common.bed
  #  2324 CRBLM_non4_common.bed
  #  3447 CRBLM_non5_common.bed
  # 29920 CRBLM_non6_common.bed
  # 79431 CRBLM_non_pool.bed
  # 39409 CRBLM_non_pool_merge.bed

# 删除209，保留CRBLM_non2_common.bed，CRBLM_non3_common.bed，CRBLM_non6_common.bed
cat  CRBLM_non2_common.bed CRBLM_non3_common.bed CRBLM_non6_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"   > CRBLM_non_pool.bed
sort -k1,1 -k2,2n CRBLM_non_pool.bed | bedtools merge -i stdin -d 50 > CRBLM_non_pool_merge.bed 
wc -l CRBLM_non*.bed
#   41473 CRBLM_non2_common.bed
#   34733 CRBLM_non3_common.bed
#   29920 CRBLM_non6_common.bed
#  106006 CRBLM_non_pool.bed
#   52879 CRBLM_non_pool_merge.bed
```
③ 合并
```bash
# CRBLM_rep1.narrowPeak CRBLM_rep2.narrowPeak CRBLM_rep3.narrowPeak CRBLM_rep4.narrowPeak

./idr.sh CRBLM_rep1.narrowPeak CRBLM_rep3.narrowPeak CRBLM_rep4.narrowPeak .
mv CRBLM_rep1_CRBLM_rep3.txt CRBLM_all2.txt 
mv CRBLM_rep1_CRBLM_rep4.txt CRBLM_all3.txt
mv CRBLM_rep3_CRBLM_rep4.txt CRBLM_all6.txt


for i in CRBLM_all*.txt; do
    awk '{if($5 >= 540) print $0}' "$i" > "${i%%.*}_common.txt"
done
for i in CRBLM_all*_common.txt; do
  cut -f 1,2,3 "$i" > ${i%%.*}.bed
done


# 删除CRBLM_rep2.narrowPeak
cat CRBLM_all2_common.bed CRBLM_all3_common.bed CRBLM_all6_common.bed | grep -v "chrUn_*" | grep -v "chrY" | grep -v "random"  > CRBLM_all_pool.bed
cat CRBLM_all_pool.bed | sort | tsv-summarize -g 1 --count 
sort -k1,1 -k2,2n CRBLM_all_pool.bed | bedtools merge -i stdin -d 50 > CRBLM_all_pool_merge.bed 
wc -l CRBLM_all*.bed
#   56132 CRBLM_all2_common.bed
#   38131 CRBLM_all3_common.bed
#   36760 CRBLM_all6_common.bed
#  130779 CRBLM_all_pool.bed
#   68379 CRBLM_all_pool_merge.bed
#  330181 total
```

5. 长度统计
```bash
cd /mnt/xuruizhi/ATAC_brain/human/IDR
mkdir -p /mnt/xuruizhi/ATAC_brain/human/peak_length
for i in *_pool_merge.bed
do
  echo $i
  awk '{print $3- $2}' $i > ../peak_length/${i%%.*}_length.txt
done
mkdir -p /mnt/d/ATAC_brain/human/peak_length
cp /mnt/xuruizhi/ATAC_brain/human/peak_length/*_pool_merge_length.txt  /mnt/d/ATAC_brain/human/peak_length/
```
```r
setwd("D:/ATAC_brain/human/peak_length")
file_list <- list.files(path = "./", pattern = ".*_pool_merge_length\\.txt$", full.names = TRUE)
summary_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE)
  summary_list[[file]] <- summary(data)
}

output_file <- "summary_output.txt" 

output_lines <- character(length(file_list))
for (i in seq_along(file_list)) {
  output_lines[i] <- paste("Summary for", file_list[i], ":\n", capture.output(print(summary_list[[file_list[i]]])), collapse = "\n")
}

writeLines(output_lines, con = output_file)


# 画直方图
summary_list <- list()
file_list <- list.files(path = "./", pattern = ".*_pool_merge_length\\.txt$", full.names = TRUE)
for (file in file_list) {
  data <- read.table(file)
  dim_data <- dim(data)
  summary_list[[file]] <- dim_data

  png(paste0(file, "_hist.png"))
  hist(abs(as.numeric(data[, 1])), breaks = 500, xlab = "Fragment length (bp)", ylab = "Frequency", main = file)
  dev.off()
}
for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}
```
6. 富集分析
```bash
mkdir -p /mnt/d/ATAC_brain/human/GO
cp /mnt/xuruizhi/ATAC_brain/human/IDR/*_pool_merge.bed  /mnt/d/ATAC_brain/human/GO/
```
```r
BiocManager::install("org.Hg.eg.db", force = TRUE)
library(biomaRt)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
region <- c("CRBLM_all")
region <- c("CRBLM_neu")
region <- c("CRBLM_non")
region <- c("PMC_all")
region <- c("PMC_neu")
region <- c("PMC_non")
region <- c("VLPFC_all")
region <- c("VLPFC_neu")
region <- c("VLPFC_non")

  setwd("D:/ATAC_brain/human/GO")
  region_peak <- readPeakFile(paste0("D:/ATAC_brain/human/GO/", region, "_pool_merge.bed"), sep = "")

  png(paste0(region, "_covplot.png"))
  covplot(region_peak)
  dev.off()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000)
  tagMatrix <- getTagMatrix(region_peak, windows = promoter)
  png(paste0(region, "_promoter.png"))
  tagHeatmap(tagMatrix, xlim = c(-1000, 1000), color = "red")
  dev.off()


  png(paste0(region, "_avg_prof_plot.png"))
  plotAvgProf(
    tagMatrix,
    xlim = c(-1000, 1000),
    xlab = "Genomic Region (5'->3')",
    ylab = "Peak Frequency"
  )
  dev.off()

  region_peakAnnolist <- annotatePeak(
    region_peak,
    tssRegion = c(-1000, 1000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
  write.table(
    as.data.frame(region_peakAnnolist),
    file = paste0(region, "_allpeak.annotation.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  png(paste0(region, "_plotAnnoPie.png"))
  plotAnnoPie(region_peakAnnolist)
  dev.off()

  region_peakAnno <- as.data.frame(region_peakAnnolist)

  ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
    return(a)
  }
  region_ensembl_id_transform <- ensembl_id_transform(region_peakAnno$ENSEMBL)
  write.csv(ensembl_id_transform(region_peakAnno$ENSEMBL), file = paste0(region, "_allpeak_geneID.tsv"), quote = FALSE)

  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
  region_biomart_ensembl_id_transform <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
    filters = "ensembl_gene_id",
    values = region_peakAnno$ENSEMBL,
    mart = mart
  )
  write.csv(region_biomart_ensembl_id_transform, file = paste0(region, "_allpeak_biomart_geneID.tsv"), quote = FALSE)

  # GO analysis and barplot
  region_biomart <- enrichGO(
    gene = region_biomart_ensembl_id_transform$entrezgene_id, 
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
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
    OrgDb = org.Hs.eg.db,
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
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(file = paste0(region, "_kegg.pdf"))
  barplot(region_kegg, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
  dev.off()
```

# 10. 每个脑区独有peak
## 10.1 三个脑区：PMC,CRBLM,VLPFC

1. total_diff
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/diff_peak1
cp /mnt/xuruizhi/ATAC_brain/human/IDR/*_pool_merge.bed  /mnt/xuruizhi/ATAC_brain/human/diff_peak1
cd /mnt/xuruizhi/ATAC_brain/human/diff_peak1
cp /mnt/xuruizhi/ATAC_brain/mouse/diff_peak1/totaldiff_peaks.sh ./

vim file_neu.list
CRBLM_neu_pool_merge.bed
PMC_neu_pool_merge.bed
VLPFC_neu_pool_merge.bed

vim file_non.list
CRBLM_non_pool_merge.bed
PMC_non_pool_merge.bed
VLPFC_non_pool_merge.bed

vim file_all.list
CRBLM_all_pool_merge.bed
PMC_all_pool_merge.bed
VLPFC_all_pool_merge.bed

for i in file_*.list
do
  echo ${i} 
  cat ${i} | bash totaldiff_peaks.sh
done
```
2. 数目统计

```bash
cd /mnt/xuruizhi/ATAC_brain/human/diff_peak1
find . -type f -name "*_pool_merge_totaldiff.bed" -exec sh -c 'mv "$0" "${0%_pool_merge_totaldiff.bed}_totaldiff.bed"' {} \;

wc -l *  
  # 68379 CRBLM_all_pool_merge.bed
  # 38040 CRBLM_all_totaldiff.bed
  # 32487 PMC_all_pool_merge.bed
  #  3766 PMC_all_totaldiff.bed
  # 34274 VLPFC_all_pool_merge.bed
  #  3149 VLPFC_all_totaldiff.bed

  # 33130 CRBLM_neu_pool_merge.bed
  # 22631 CRBLM_neu_totaldiff.bed
  # 14971 PMC_neu_pool_merge.bed
  #  2690 PMC_neu_totaldiff.bed
  # 16208 VLPFC_neu_pool_merge.bed
  #  3418 VLPFC_neu_totaldiff.bed

  # 52879 CRBLM_non_pool_merge.bed
  # 28673 CRBLM_non_totaldiff.bed
  # 23898 PMC_non_pool_merge.bed
  #  1174 PMC_non_totaldiff.bed
  # 21944 VLPFC_non_pool_merge.bed
  #   767 VLPFC_non_totaldiff.bed
```

3. 长度计算

```bash
mkdir -p /mnt/d/ATAC_brain/human/GO_totaldiff
cd /mnt/xuruizhi/ATAC_brain/human/diff_peak1
ls *_totaldiff.bed | while read id
do
  awk '{print ($3-$2)}' ${id} > ${id%%.*}_length.txt
done
cp /mnt/xuruizhi/ATAC_brain/human/diff_peak1/*_length.txt /mnt/d/ATAC_brain/human/GO_totaldiff
cp /mnt/xuruizhi/ATAC_brain/human/diff_peak1/*_totaldiff.bed /mnt/d/ATAC_brain/human/GO_totaldiff
```

```r
setwd("D:/ATAC_brain/human/GO_totaldiff")
file_list <- list.files(path = "./", pattern = ".*_totaldiff_length\\.txt$", full.names = TRUE)
# 其他同上代码
```

4. 富集分析
```bash
cp /mnt/xuruizhi/ATAC_brain/human/diff_peak1/*_totaldiff.bed /mnt/d/ATAC_brain/human/GO_totaldiff/
```
```r
region <- c("CRBLM_all")
region <- c("CRBLM_neu")
region <- c("CRBLM_non")
region <- c("PMC_all")
region <- c("PMC_neu")
region <- c("PMC_non")
region <- c("VLPFC_all")
region <- c("VLPFC_neu")
region <- c("VLPFC_non")

setwd("D:/ATAC_brain/human/GO_totaldiff/")
region_peak <- readPeakFile(paste0("D:/ATAC_brain/human/GO_totaldiff/", region, "_totaldiff.bed"), sep = "")
```

# 11. 每个脑区共有peak
```bash
# 建目录
mkdir -p /mnt/xuruizhi/ATAC_brain/human/common_0.5

cd /mnt/xuruizhi/ATAC_brain/human/IDR
cp ./*_pool_merge.bed /mnt/xuruizhi/ATAC_brain/human/common_0.5/
```

```bash
cd /mnt/xuruizhi/ATAC_brain/human/common_0.5/
vim common_peak_0.5.sh
#!/bin/bash
echo "请输入文件名列表，以空格或换行符分隔："
read -r file_list
IFS=$' \n' read -ra file_array <<< "$file_list"
count=1
for ((i=0; i<${#file_array[@]}; i++)); do
    for ((j=i+1; j<${#file_array[@]}; j++)); do
        sample_a="${file_array[i]}"
        sample_b="${file_array[j]}"
        output_file="${count}$(basename "${sample_a%.*}_${sample_b%.*}_0.5.bed")"
        ((count++))
        bedtools intersect -a "$sample_a" -b "$sample_b" -sorted -f 0.5 -r > "$output_file"
    done
done

bash common_peak_0.5.sh
# CRBLM_all_pool_merge.bed PMC_all_pool_merge.bed VLPFC_all_pool_merge.bed
# CRBLM_neu_pool_merge.bed PMC_neu_pool_merge.bed VLPFC_neu_pool_merge.bed
# CRBLM_non_pool_merge.bed PMC_non_pool_merge.bed VLPFC_non_pool_merge.bed


wc -l *
  #  20848 1CRBLM_all_pool_merge_PMC_all_pool_merge_0.5.bed
  #   7251 1CRBLM_neu_pool_merge_PMC_neu_pool_merge_0.5.bed
  #  17542 1CRBLM_non_pool_merge_PMC_non_pool_merge_0.5.bed
  #  21560 2CRBLM_all_pool_merge_VLPFC_all_pool_merge_0.5.bed
  #   6808 2CRBLM_neu_pool_merge_VLPFC_neu_pool_merge_0.5.bed
  #  16199 2CRBLM_non_pool_merge_VLPFC_non_pool_merge_0.5.bed
  #  22062 3PMC_all_pool_merge_VLPFC_all_pool_merge_0.5.bed
  #   9457 3PMC_neu_pool_merge_VLPFC_neu_pool_merge_0.5.bed
  #  16371 3PMC_non_pool_merge_VLPFC_non_pool_merge_0.5.bed
  #  68379 CRBLM_all_pool_merge.bed
  #  33130 CRBLM_neu_pool_merge.bed
  #  52879 CRBLM_non_pool_merge.bed
  #  32487 PMC_all_pool_merge.bed
  #  14971 PMC_neu_pool_merge.bed
  #  23898 PMC_non_pool_merge.bed
  #  34274 VLPFC_all_pool_merge.bed
  #  16208 VLPFC_neu_pool_merge.bed
  #  21944 VLPFC_non_pool_merge.bed


vim common_peak_nextstep.sh
#!/bin/bash
echo "请输入文件名列表，以空格或换行符分隔："
read -r file_list
IFS=$' \n' read -ra file_array <<< "$file_list"
first_file="${file_array[0]}"
for ((i=1; i<${#file_array[@]}; i++)); do
    output_file="cross$((i)).bed"
    input_file="${file_array[i]}"
    bedtools intersect -a "$first_file" -b "$input_file" > "$output_file"
    first_file="$output_file"
done


./common_peak_nextstep.sh 
# 1CRBLM_all_pool_merge_PMC_all_pool_merge_0.5.bed 2CRBLM_all_pool_merge_VLPFC_all_pool_merge_0.5.bed 3PMC_all_pool_merge_VLPFC_all_pool_merge_0.5.bed
# 结果 15685 cross2.bed

# 1CRBLM_neu_pool_merge_PMC_neu_pool_merge_0.5.bed 2CRBLM_neu_pool_merge_VLPFC_neu_pool_merge_0.5.bed 3PMC_neu_pool_merge_VLPFC_neu_pool_merge_0.5.bed
# 结果 5068 cross2.bed

# 1CRBLM_non_pool_merge_PMC_non_pool_merge_0.5.bed 2CRBLM_non_pool_merge_VLPFC_non_pool_merge_0.5.bed 3PMC_non_pool_merge_VLPFC_non_pool_merge_0.5.bed
# 12620 cross2.bed
```
* 回比到原有peak集
```bash
ls *non_pool_merge.bed | while read id 
do
bedtools intersect -wa -u -a ${id} -b cross2.bed -sorted > ${id%%.*}_commonpeak.bed
done
# 12620  *non_pool_merge_commonpeak.bed 
# 都是一一对应的关系
```
* 统计交集长度
```bash
mv cross2.bed all_common.bed
mv cross2.bed neu_common.bed
mv cross2.bed non_common.bed
```

```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/common_peak_length
cd /mnt/xuruizhi/ATAC_brain/human/common_0.5
for i in *_common.bed
do
  echo $i
  awk '{print $3- $2}' $i > ../common_peak_length/${i%%.*}_length.txt
done
mkdir -p /mnt/d/ATAC_brain/human/GO_common
cp /mnt/xuruizhi/ATAC_brain/human/common_peak_length/*  /mnt/d/ATAC_brain/human/GO_common
cp /mnt/xuruizhi/ATAC_brain/human/common_0.5/*_common.bed  /mnt/d/ATAC_brain/human/GO_common
```
```r
setwd("D:/ATAC_brain/human/GO_common")
file_list <- list.files(path = "./", pattern = ".*_common_length\\.txt$", full.names = TRUE)
```
* 富集分析

```r
setwd("D:/ATAC_brain/human/GO_common")
region <- c("all")
region <- c("neu")
region <- c("non")

region_peak <- readPeakFile(paste0("D:/ATAC_brain/human/GO_common/", region, "_common.bed"), sep = "")

  png(paste0(region, "_covplot.png"))
  covplot(region_peak)
  dev.off()

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  promoter <- getPromoters(TxDb = txdb, upstream = 1000, downstream = 1000)
  tagMatrix <- getTagMatrix(region_peak, windows = promoter)
  png(paste0(region, "_promoter.png"))
  tagHeatmap(tagMatrix, xlim = c(-1000, 1000), color = "red")
  dev.off()


  png(paste0(region, "_avg_prof_plot.png"))
  plotAvgProf(
    tagMatrix,
    xlim = c(-1000, 1000),
    xlab = "Genomic Region (5'->3')",
    ylab = "Peak Frequency"
  )
  dev.off()

  region_peakAnnolist <- annotatePeak(
    region_peak,
    tssRegion = c(-1000, 1000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
  write.table(
    as.data.frame(region_peakAnnolist),
    file = paste0(region, "_allpeak.annotation.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  png(paste0(region, "_plotAnnoPie.png"))
  plotAnnoPie(region_peakAnnolist)
  dev.off()

  region_peakAnno <- as.data.frame(region_peakAnnolist)

  ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
    return(a)
  }
  region_ensembl_id_transform <- ensembl_id_transform(region_peakAnno$ENSEMBL)
  write.csv(ensembl_id_transform(region_peakAnno$ENSEMBL), file = paste0(region, "_allpeak_geneID.tsv"), quote = FALSE)

  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
  region_biomart_ensembl_id_transform <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "description"),
    filters = "ensembl_gene_id",
    values = region_peakAnno$ENSEMBL,
    mart = mart
  )
  write.csv(region_biomart_ensembl_id_transform, file = paste0(region, "_allpeak_biomart_geneID.tsv"), quote = FALSE)

  # GO analysis and barplot
  region_biomart <- enrichGO(
    gene = region_biomart_ensembl_id_transform$entrezgene_id, 
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
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
    OrgDb = org.Hs.eg.db,
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
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(file = paste0(region, "_kegg.pdf"))
  barplot(region_kegg, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
  dev.off()
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
SRR21161738.
SRR21161739
SRR21161881
SRR21161882..
SRR21161914.
SRR21161915
# VLPFC
SRR21161734.
SRR21161735
SRR21161750
SRR21161751
SRR21161759.
SRR21161760
SRR21161765
SRR21161766
SRR21161909.
SRR21161910
SRR21161931
SRR21161932
# CERE
SRR21161742.
SRR21161743.
SRR21161767
SRR21161768..
SRR21161780..
SRR21161781.
SRR21161961
SRR21161962.


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


## 12.2 比对前质控
```bash
# 在本地转换格式
mkdir -p /mnt/xuruizhi/RNA_brain/human/sequence
mkdir -p /mnt/xuruizhi/RNA_brain/human/fastqc
mkdir -p /mnt/xuruizhi/RNA_brain/human/trim
mkdir -p /mnt/xuruizhi/RNA_brain/human/fastqc_again

cd /mnt/xuruizhi/RNA_brain/human/sequence

vim human.sh
#!/usr/bin bash
#sra2fq
fastq-dump --gzip --split-3 -O /mnt/xuruizhi/RNA_brain/human/sequence /mnt/xuruizhi/RNA_brain/human/sra/{}/{}.sra
# fastqc
fastqc -o ../fastqc {}_1.fastq.gz
fastqc -o ../fastqc {}_2.fastq.gz
# trim
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o ../trim  {}_1.fastq.gz  {}_2.fastq.gz
# fatsqc_again
fastqc -o ../fastqc_again ../trim/{}_1_val_1.fq.gz
fastqc -o ../fastqc_again ../trim/{}_2_val_2.fq.gz



cat ../sra/1.list | while read id
do
  sed "s/{}/${id}/g" human.sh > ${id}_prealign.sh
done
cat ../sra/1.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  bash {}_prealign.sh >> ../trim/trim_fastqc.log 2>&1"


cd ../fastqc
multiqc .
cd ../fastqc_again
multiqc .

# 因为CG含量质控不合格，对文件5'端进行5bp删除
cd /scratch/wangq/xrz/RNA_brain/human/trim
mkdir -p ../cut

cat 1.list | parallel --no-run-if-empty --linebuffer -k -j 6 " 
  cutadapt --cut -5 {}_1_val_1.fq.gz  -o ../cut/{}_1_cut_1.fq.gz >> ../cut/cut_5bp.log 2>&1"
fastqc -t 6 -o ../fastqc2 *.fq.gz


# cutadapt双端3’删除
cd /scratch/wangq/xrz/RNA_brain/human/trim
mkdir -p ../cut2
bsub -q mpi -n 24 -o ../cut2 '
cat 1.list | parallel --no-run-if-empty --linebuffer -k -j 20 " 
  cutadapt -u 5  -o ../cut2/{}_1_cut_1.fq.gz -p ../cut2/{}_2_cut_2.fq.gz {}_1_val_1.fq.gz {}_2_val_2.fq.gz  >> ../cut2/cut_5bp.log 2>&1"'

# 处理完发现，处理后结果不如处理前，还是保留使用原结果
```

## 12.3 比对

1. 参考基因组  
```bash
# 人类 hg38 
mkdir -p /mnt/xuruizhi/RNA_brain/human/genome
cd /mnt/xuruizhi/RNA_brain/human/genome
wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz

# 传输超算
cd /scratch/wangq/xrz/RNA_brain/human/genome
tar xvzf hg38_genome.tar.gz
rm hg38_genome.tar.gz
```

2. 比对
```bash
mkdir -p /scratch/wangq/xrz/RNA_brain/human/sort_bam
mkdir -p /scratch/wangq/xrz/RNA_brain/human/align
cp /scratch/wangq/xrz/RNA_brain/human/*.list ../trim
cd /scratch/wangq/xrz/RNA_brain/human/trim

vim human.sh
#!/usr/bin bash
hisat2 -p 20 -t -x ../genome/hg38/genome -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz -S ../align/{}.sam 2>../align/{}.log 
samtools sort -@ 20 ../align/{}.sam > ../sort_bam/{}.sort.bam
samtools index -@ 8 ../sort_bam/{}.sort.bam
samtools flagstat  -@ 8 ../sort_bam/{}.sort.bam > ../sort_bam/{}.raw.stat



cat 1.list  | while read id
do 
  sed "s/{}/${id}/g" human.sh > ${id}_align.sh
  bsub -q mpi -n 24  "
  bash  ${id}_align.sh >> ./align.log 2>&1"
done
```

## 12.4 表达量统计
1. 下载注释文件
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/annotation
cd /mnt/xuruizhi/RNA_brain/human/annotation
# UCSC注释文件比较混乱，改为ensembl，和小鼠保持一致
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
gzip -dc hg38.ensGene.gtf.gz > hg38.gtf
```

2. 统计
```bash
cd /mnt/xuruizhi/RNA_brain/human/sort_bam
mkdir -p /mnt/xuruizhi/RNA_brain/human/HTseq

parallel -j 8 "
    htseq-count -s no -r pos -f bam {1}.sort.bam ../annotation/hg38.gtf > ../HTseq/{1}.count  2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

#下载HTSeq
cd /mnt/d/biosoft
wget https://files.pythonhosted.org/packages/de/13/cbeb753eb29d31bb252bc18443e79a3739bf5bd3dbc99e65fcf13db81be1/HTSeq-2.0.2.tar.gz
tar -zxvf HTSeq-2.0.2.tar.gz
rsync -avP /mnt/d/biosoft/HTSeq-2.0.2 wangq@202.119.37.251:~/xuruizhi/biosoft
# 进入超算
cd ~/xuruizhi/biosoft/HTSeq-2.0.2
python setup.py install –user  #安装
# 因为超算版本太低还是算了


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

cd /mnt/d/RNA_brain/human/Deseq2
vim 1_neu.list
SRR21161731
SRR21161739
SRR21161882
SRR21161915
SRR21161735
SRR21161751
SRR21161760
SRR21161766
SRR21161910
SRR21161932
SRR21161743
SRR21161768
SRR21161781
SRR21161962

vim 1_non.list
SRR21161730
SRR21161738
SRR21161881
SRR21161914
SRR21161734
SRR21161750
SRR21161759
SRR21161765
SRR21161909
SRR21161931
SRR21161742
SRR21161767
SRR21161780
SRR21161961


while read -r i
do
  cp ../HTseq/${i}.count ./
done < 1_non.list
```
```r
rm(list=ls())
setwd("D:/RNA_brain/human/Deseq2")

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

data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]],by="gene_id")
}
write.csv(data_merge, "non.csv", quote = FALSE, row.names = FALSE)
# neu同理
```
2. 初步分析
* coldata
```bash
cd /mnt/d/RNA_brain/human/Deseq2 
vim coldata_neu.csv
"ids","state","condition","treatment"
"SRR21161731","WT","neuron","PSM"
"SRR21161739","WT","neuron","PSM"
"SRR21161882","WT","neuron","PSM"
"SRR21161915","WT","neuron","PSM"
"SRR21161735","WT","neuron","VLPFC"
"SRR21161751","WT","neuron","VLPFC"
"SRR21161760","WT","neuron","VLPFC"
"SRR21161766","WT","neuron","VLPFC"
"SRR21161910","WT","neuron","VLPFC"
"SRR21161932","WT","neuron","VLPFC"
"SRR21161743","WT","neuron","CRBLM"
"SRR21161768","WT","neuron","CRBLM"
"SRR21161781","WT","neuron","CRBLM"
"SRR21161962","WT","neuron","CRBLM"

vim coldata_non.csv
"ids","state","condition","treatment"
"SRR21161730","WT","non-neuron","PSM"
"SRR21161738","WT","non-neuron","PSM"
"SRR21161881","WT","non-neuron","PSM"
"SRR21161914","WT","non-neuron","PSM"
"SRR21161734","WT","non-neuron","VLPFC"
"SRR21161750","WT","non-neuron","VLPFC"
"SRR21161759","WT","non-neuron","VLPFC"
"SRR21161765","WT","non-neuron","VLPFC"
"SRR21161909","WT","non-neuron","VLPFC"
"SRR21161931","WT","non-neuron","VLPFC"
"SRR21161742","WT","non-neuron","CRBLM"
"SRR21161767","WT","non-neuron","CRBLM"
"SRR21161780","WT","non-neuron","CRBLM"
"SRR21161961","WT","non-neuron","CRBLM"
```
```r
BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# neu
dataframe <- read.csv("neu.csv", header=TRUE, row.names = 1)
countdata <- dataframe[-(1:5),]
countdata <- countdata[rowSums(countdata) > 0,]
head(countdata)
# 导入coltdata文件
coldata <- read.table("coldata_neu.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds



# PCA分析 
# 归一化
rld <- rlog(dds, blind=FALSE)
# intgroup分组
pcaData <- plotPCA(rld, intgroup=c("treatment"),returnData = T) 
pcaData <- pcaData[order(pcaData$treatment,decreasing=F),]
table(pcaData$treatment)
# PCA1
plot(pcaData[, 1:2], pch = 19, col = c(rep("red", 4), rep("green", 4), rep("blue", 6)))
text(pcaData[, 1], pcaData[, 2], row.names(pcaData), cex = 0.75, font = 1)
legend(0, 0, inset = 0.02, pt.cex = 1.5, legend = c("CRBLM", "PSM", "VLPFC"), col = c("red", "green", "blue"), pch = 19, cex = 0.75, bty = "n")
# PCA2
plotPCA(rld, intgroup="treatment") + ylim(-30, 30)+text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=0.5, font = 1)
# 聚类热图
library("RColorBrewer")
gene_data_transform <- assay(rld)
sampleDists <- dist(t(gene_data_transform))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$treatment
colnames(sampleDistMatrix) <- rld$treatment
colors <- colorRampPalette(rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
         
3. 差异分析，以CRBLM为control
```r
dds$treatment <- factor(as.vector(dds$treatment), levels = c("CRBLM","PSM","VLPFC")) 
dds$treatment
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"   "treatment_PSM_vs_CRBLM"   "treatment_VLPFC_vs_CRBLM"
```



* PSM_vs_CRBLM
```r
result <- results(dds, name="treatment_PSM_vs_CRBLM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 52215 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9032, 17%
# LFC < 0 (down)     : 4912, 9.4%
# outliers [1]       : 140, 0.27%
# low counts [2]     : 13160, 25%
# (mean count < 1)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 24971 13944 
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 11268    6
write.csv(diff_gene, file="neu_diff_PSM_vs_CRBLM.csv", quote = F)
```


* VLPFC_vs_CRBLM
```r
result <- results(dds, name="treatment_VLPFC_vs_CRBLM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 52215 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9603, 18%
# LFC < 0 (down)     : 5202, 10%
# outliers [1]       : 140, 0.27%
# low counts [2]     : 13160, 25%
# (mean count < 1)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  24110 14805
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 11652   6
write.csv(diff_gene, file="neu_diff_VLPFC_vs_CRBLM.csv", quote = F)
```


4. 差异分析，以PSM为control
```r
# 其他和前文一致
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds$treatment <- factor(as.vector(dds$treatment), levels = c("PSM","CRBLM","VLPFC")) 
dds$treatment

dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"              "treatment_CRBLM_vs_PSM"
# [3] "treatment_VLPFC_vs_PSM"
```



* CRBLM_vs_PSM
```r
result <- results(dds, name="treatment_CRBLM_vs_PSM", pAdjustMethod = "fdr", alpha = 0.05)
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
# 24971 13944

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #11268    6
write.csv(diff_gene, file="neu_diff_CRBLM_vs_PSM.csv", quote = F)
```


* VLPFC_vs_PSM
```r
result <- results(dds, name="treatment_VLPFC_vs_PSM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 52215 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0019%
# outliers [1]       : 140, 0.27%
# low counts [2]     : 0, 0%
# (mean count < 0)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  52074     1

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #1    6
write.csv(diff_gene, file="neu_diff_VLPFC_vs_PSM.csv", quote = F)
```

5. 差异分析，以VLPFC为control
```r
# 其他和前文一致
dds$treatment <- factor(as.vector(dds$treatment), levels = c("VLPFC","PSM","CRBLM")) 
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"                "treatment_PSM_vs_VLPFC"  
# [3] "treatment_CRBLM_vs_VLPFC"
```



* PSM_vs_VLPFC
```r
result <- results(dds, name="treatment_PSM_vs_VLPFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order) 
# out of 52215 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0019%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 140, 0.27%
# low counts [2]     : 0, 0%
# (mean count < 0)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 52074     1 

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #1   6
write.csv(diff_gene, file="neu_diff_PSM_vs_VLPFC.csv", quote = F)
```


* CRBLM_vs_VLPFC
```r
result <- results(dds, name="treatment_CRBLM_vs_VLPFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 52215 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5202, 10%
# LFC < 0 (down)     : 9603, 18%
# outliers [1]       : 140, 0.27%
# low counts [2]     : 13160, 25%
# (mean count < 1)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 24110 14805

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #11652    6
write.csv(diff_gene, file="neu_diff_CRBLM_vs_VLPFC.csv", quote = F)
```


4. 脑区差异基因

找到某一脑区对另外两个脑区来说，都差异表达的基因。

```bash
cd /mnt/d/RNA_brain/human/Deseq2 
# diff_HIPP_DG.csv diff_HIPP_PFC.csv

# 区分up/down
for i in neu_diff_*.csv
do
sed 's/,/\t/g' "$i" | tail -n +2  > "${i%.csv}.txt"
done
dos2unix *.txt

for i in neu_diff*.txt 
do
  echo " ==> $i <== " 
  tsv-filter --is-numeric 3 --gt 3:0 $i > ${i%%.*}_up.txt
  tsv-filter --is-numeric 3 --lt 3:0 $i > ${i%%.*}_down.txt
done

wc -l *_up.txt
  #  3169 neu_diff_CRBLM_vs_PSM_up.txt
  #  3417 neu_diff_CRBLM_vs_VLPFC_up.txt
  #  8099 neu_diff_PSM_vs_CRBLM_up.txt
  #     1 neu_diff_PSM_vs_VLPFC_up.txt
  #  8235 neu_diff_VLPFC_vs_CRBLM_up.txt
  #     0 neu_diff_VLPFC_vs_PSM_up.txt

wc -l *_down.txt
  #  8099 neu_diff_CRBLM_vs_PSM_down.txt
  #  8235 neu_diff_CRBLM_vs_VLPFC_down.txt
  #  3169 neu_diff_PSM_vs_CRBLM_down.txt
  #     0 neu_diff_PSM_vs_VLPFC_down.txt
  #  3417 neu_diff_VLPFC_vs_CRBLM_down.txt
  #     1 neu_diff_VLPFC_vs_PSM_down.txt
```

* CRBLM对于其他两个脑区的差异基因
```bash
# up
awk 'NR==FNR {a[$1]=1; next} a[$1]' neu_diff_CRBLM_vs_PSM_up.txt neu_diff_CRBLM_vs_VLPFC_up.txt > neu_diff_CRBLM_up.txt #2500
# down
awk 'NR==FNR {a[$1]=1; next} a[$1]' neu_diff_CRBLM_vs_PSM_down.txt neu_diff_CRBLM_vs_VLPFC_down.txt > neu_diff_CRBLM_down.txt # 7187

cat neu_diff_CRBLM_up.txt neu_diff_CRBLM_down.txt > neu_diff_CRBLM.txt
```
```r
library(biomaRt)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd("D:/RNA_brain/human/Deseq2")
region <- c("neu_diff_CRBLM_up")
region <- c("neu_diff_CRBLM_down")
region <- c("neu_diff_CRBLM")


data <- read.table(paste0(region,".txt"), header=FALSE)
  
  ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
    return(a)

  }
  region_ensembl_id_transform <- ensembl_id_transform(data$V1)
  write.csv(ensembl_id_transform(data$V1), file =  paste0(region,"_ensemblID.tsv"))

  region_transform <- enrichGO(
    gene = data$V1, 
    keyType = "ENSEMBL",
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
# 两种都可以
  region_transform <- enrichGO(
    gene = region_ensembl_id_transform$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
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
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(file = paste0(region, "_kegg.pdf"),width = 80, height = 120)
  barplot(region_kegg, showCategory = 20, font.size = 120,title = "KEGG Pathway Enrichment Analysis")
  dev.off()
```

* PSM和VLPFC没有找到对于其他两个脑区共有的差异基因
```r
# neu_diff_PSM_vs_CRBLM_up.txt
region <- c("neu_diff_VLPFC_vs_CRBLM_down")
data <- read.table(paste0(region,".txt"), header=FALSE)
```








































```r
BiocManager::install("DESeq2")
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

# neu
dataframe <- read.csv("non.csv", header=TRUE, row.names = 1)
countdata <- dataframe[-(1:5),]
countdata <- countdata[rowSums(countdata) > 0,]
head(countdata)
# 导入coltdata文件
coldata <- read.table("coldata_non.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds



# PCA分析 
# 归一化
rld <- rlog(dds, blind=FALSE)
# intgroup分组
pcaData <- plotPCA(rld, intgroup=c("treatment"),returnData = T) 
pcaData <- pcaData[order(pcaData$treatment,decreasing=F),]
table(pcaData$treatment)
# PCA1
plot(pcaData[, 1:2], pch = 19, col = c(rep("red", 4), rep("green", 4), rep("blue", 6)))
text(pcaData[, 1], pcaData[, 2], row.names(pcaData), cex = 0.75, font = 1)
legend(0, 0, inset = 0.02, pt.cex = 1.5, legend = c("CRBLM", "PSM", "VLPFC"), col = c("red", "green", "blue"), pch = 19, cex = 0.75, bty = "n")
# PCA2
plotPCA(rld, intgroup="treatment") + ylim(-40, 40)+text(pcaData[,1],pcaData[,2],row.names(pcaData),cex=0.5, font = 1)
# 聚类热图
library("RColorBrewer")
gene_data_transform <- assay(rld)
sampleDists <- dist(t(gene_data_transform))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$treatment
colnames(sampleDistMatrix) <- rld$treatment
colors <- colorRampPalette(rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
         
3. 差异分析，以CRBLM为control
```r
dds$treatment <- factor(as.vector(dds$treatment), levels = c("CRBLM","PSM","VLPFC")) 
dds$treatment
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"   "treatment_PSM_vs_CRBLM"   "treatment_VLPFC_vs_CRBLM"
```



* PSM_vs_CRBLM
```r
result <- results(dds, name="treatment_PSM_vs_CRBLM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 514, 1%
# LFC < 0 (down)     : 975, 1.9%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 17954, 35%
# (mean count < 3)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 31882  1489
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 1337    6
write.csv(diff_gene, file="non_diff_PSM_vs_CRBLM.csv", quote = F)
```


* VLPFC_vs_CRBLM
```r
result <- results(dds, name="treatment_VLPFC_vs_CRBLM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 707, 1.4%
# LFC < 0 (down)     : 925, 1.8%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 17954, 35%
# (mean count < 3

table(result_order$padj<0.05)
# FALSE  TRUE 
#  31739  1632 
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 11652   6
write.csv(diff_gene, file="non_diff_VLPFC_vs_CRBLM.csv", quote = F)
```


4. 差异分析，以PSM为control
```r
# 其他和前文一致
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds$treatment <- factor(as.vector(dds$treatment), levels = c("PSM","CRBLM","VLPFC")) 
dds$treatment

dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"              "treatment_CRBLM_vs_PSM"
# [3] "treatment_VLPFC_vs_PSM"
```



* CRBLM_vs_PSM
```r
result <- results(dds, name="treatment_CRBLM_vs_PSM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order) # 和前文HIPP vs DG刚好反过来
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 975, 1.9%
# LFC < 0 (down)     : 514, 1%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 17954, 35%
# (mean count < 3)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 31882  1489

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #1337    6
write.csv(diff_gene, file="non_diff_CRBLM_vs_PSM.csv", quote = F)
```


* VLPFC_vs_PSM
```r
result <- results(dds, name="treatment_VLPFC_vs_PSM", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0019%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 0, 0%
# (mean count < 0)

table(result_order$padj<0.05)
# FALSE  TRUE 
#  51324     1

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #1    6
write.csv(diff_gene, file="non_diff_VLPFC_vs_PSM.csv", quote = F)
```

5. 差异分析，以VLPFC为control
```r
# 其他和前文一致
dds$treatment <- factor(as.vector(dds$treatment), levels = c("VLPFC","PSM","CRBLM")) 
dds <- DESeq(dds)
resultsNames(dds) 
# [1] "Intercept"                "treatment_PSM_vs_VLPFC"  
# [3] "treatment_CRBLM_vs_VLPFC"
```



* PSM_vs_VLPFC
```r
result <- results(dds, name="treatment_PSM_vs_VLPFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order) 
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.0019%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 0, 0%
# (mean count < 0)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 51324     1 

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   #1   6
write.csv(diff_gene, file="non_diff_PSM_vs_VLPFC.csv", quote = F)
```


* CRBLM_vs_VLPFC
```r
result <- results(dds, name="treatment_CRBLM_vs_VLPFC", pAdjustMethod = "fdr", alpha = 0.05)
result_order <- result[order(result$pvalue),]
summary(result_order)
# out of 51657 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 925, 1.8%
# LFC < 0 (down)     : 707, 1.4%
# outliers [1]       : 332, 0.64%
# low counts [2]     : 17954, 35%
# (mean count < 3)

table(result_order$padj<0.05)
# FALSE  TRUE 
# 31739  1632

diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene)   # 1459    6
write.csv(diff_gene, file="non_diff_CRBLM_vs_VLPFC.csv", quote = F)
```


4. 脑区差异基因

找到某一脑区对另外两个脑区来说，都差异表达的基因。

```bash
cd /mnt/d/RNA_brain/human/Deseq2 

# 区分up/down
for i in non_diff_*.csv
do
sed 's/,/\t/g' "$i" | tail -n +2  > "${i%.csv}.txt"
done
dos2unix *.txt

for i in non_diff*.txt 
do
  echo " ==> $i <== " 
  tsv-filter --is-numeric 3 --gt 3:0 $i > ${i%%.*}_up.txt
  tsv-filter --is-numeric 3 --lt 3:0 $i > ${i%%.*}_down.txt
done

wc -l non*_up.txt
  #  898 non_diff_CRBLM_vs_PSM_up.txt
  #  825 non_diff_CRBLM_vs_VLPFC_up.txt
  #  439 non_diff_PSM_vs_CRBLM_up.txt
  #    0 non_diff_PSM_vs_VLPFC_up.txt
  #  634 non_diff_VLPFC_vs_CRBLM_up.txt
  #    1 non_diff_VLPFC_vs_PSM_up.txt

wc -l non*_down.txt
  #  439 non_diff_CRBLM_vs_PSM_down.txt
  #  634 non_diff_CRBLM_vs_VLPFC_down.txt
  #  898 non_diff_PSM_vs_CRBLM_down.txt
  #    1 non_diff_PSM_vs_VLPFC_down.txt
  #  825 non_diff_VLPFC_vs_CRBLM_down.txt
  #    0 non_diff_VLPFC_vs_PSM_down.txt
```

* CRBLM对于其他两个脑区的差异基因
```bash
# up
awk 'NR==FNR {a[$1]=1; next} a[$1]' non_diff_CRBLM_vs_PSM_up.txt non_diff_CRBLM_vs_VLPFC_up.txt > non_diff_CRBLM_up.txt #589
# down
awk 'NR==FNR {a[$1]=1; next} a[$1]' non_diff_CRBLM_vs_PSM_down.txt non_diff_CRBLM_vs_VLPFC_down.txt > non_diff_CRBLM_down.txt # 293

cat non_diff_CRBLM_up.txt non_diff_CRBLM_down.txt > non_diff_CRBLM.txt
```
```r
library(biomaRt)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd("D:/RNA_brain/human/Deseq2")
region <- c("non_diff_CRBLM_up")
region <- c("non_diff_CRBLM_down")
region <- c("non_diff_CRBLM")


data <- read.table(paste0(region,".txt"), header=FALSE)
  
  ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
    return(a)

  }
  region_ensembl_id_transform <- ensembl_id_transform(data$V1)
  write.csv(ensembl_id_transform(data$V1), file =  paste0(region,"_ensemblID.tsv"))

  region_transform <- enrichGO(
    gene = data$V1, 
    keyType = "ENSEMBL",
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
# 两种都可以
  region_transform <- enrichGO(
    gene = region_ensembl_id_transform$ENTREZID, 
    keyType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
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
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(file = paste0(region, "_kegg.pdf"),width = 80, height = 120)
  barplot(region_kegg, showCategory = 20, font.size = 120,title = "KEGG Pathway Enrichment Analysis")
  dev.off()
```

* PSMY和VLPFC没有找到对于其他两个脑区共有的差异基因
```r
# neu_diff_PSM_vs_CRBLM_up.txt
region <- c("neu_diff_VLPFC_vs_CRBLM_down")
data <- read.table(paste0(region,".txt"), header=FALSE)
```


