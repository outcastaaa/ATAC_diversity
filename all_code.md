# 1. 建目录
```bash
cd /mnt/xuruizhi # 挂载到NAS
mkdir -p /mnt/xuruizhi/brain/sequence/mouse
mkdir -p /mnt/xuruizhi/brain/sequence/rat
mkdir -p /mnt/xuruizhi/brain/sequence/macaque
mkdir -p /mnt/xuruizhi/brain/sequence/human
```

# 2. 下载数据

```bash
# mouse 21个
cd /mnt/xuruizhi/brain/sequence
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
SRR14362275
SRR14362276
SRR14362281
SRR14362282
SRR14614715
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443553
SRR13443549
SRR13443554
EOF

prefetch SRR11179779  SRR13049359  SRR13049362  SRR14362275  SRR14362282  SRR3595212 \
SRR11179780  SRR13049360  SRR13049363  SRR14362276  SRR14614715  SRR3595213 \
SRR11179781  SRR13049361  SRR13049364  SRR14362281  SRR3595211   SRR3595214 \
SRR13443553 SRR13443549 SRR13443554
cd ~/data/sra
cp *.sra /mnt/xuruizhi/brain/sequence/mouse
rm *.sra



# rat 12个
cd /mnt/xuruizhi/brain/sequence
cat >RAT.list <<EOF
SRR9050823
SRR9050824
SRR9050825
SRR9050826
SRR9050827
SRR9050828
SRR9050829
SRR9050830
SRR14136034
SRR14136035
SRR14136036
SRR14136037
EOF

prefetch SRR9050823 SRR9050824 SRR9050825 SRR9050826 \
SRR9050827 SRR9050828 SRR9050829 SRR9050830 \
SRR14136034 SRR14136035 SRR14136036 SRR14136037
cd ~/data/sra
cp *.sra /mnt/xuruizhi/brain/sequence/rat 
rm *.sra

# macaque只有peak，无原始文件


# human
cd /mnt/xuruizhi/brain/sequence
cat >HUMAN.list <<EOF
SRR5367709
SRR5367714
SRR5367715
SRR5367718
SRR5367720
SRR5367725
SRR5367726
SRR5367729
SRR5367732
SRR5367733
SRR5367737
SRR5367740
SRR5367742
SRR5367744
SRR5367745
SRR5367749
SRR5367751
SRR5367754
SRR5367756
SRR5367757
SRR5367762
SRR5367763
SRR5367766
SRR5367768
SRR5367769
SRR5367774
SRR5367775
SRR5367778
SRR5367780
SRR5367785
SRR5367789
SRR5367791
SRR5367795
SRR5367798
SRR5367800
SRR5367801
SRR5367802
SRR5367807
SRR5367808
SRR5367811
SRR5367813
SRR5367814
SRR5367815
SRR5367820
SRR5367821
SRR5367824
EOF


cd /mnt/xuruizhi/brain/sequence
cat HUMAN.list | while read id;
do 
  echo $id
  prefetch $id
  cp /home/xuruizhi/data/sra/$id.sra /mnt/xuruizhi/brain/sequence/human
  1>log_human.txt 2>&1
done
```

# 在超算中转换格式
```bash
# 代码示例
#批量转换，echo是打印命令，while循环的意义是生成脚本
cd ~/xuruizhi/brain
mkdir -p ~/xuruizhi/brain/fastq
fqdir= ~/xuruizhi/brain/fastq/物种
ls *.sra | while read id
do
 echo "fastq-dump --gzip --split-3 -O ${fqdir} ${id}"
done >sra2fq.sh
# 提交后台运行命令，脚本文件后缀为.sh，日志文件后缀为.log，运行脚本的命令为sh
sh sra2fq.sh > sra2fq.log 
#查看输出的fastq的gz压缩文件，用zless命令
zless -S SRRxxx.fastq.gz



# 超算跑的有问题,在本地
mkdir -p /mnt/xuruizhi/brain/fastq
cd /mnt/xuruizhi/brain/fastq
cat >species.list <<EOF
mouse
rat
human
EOF
cat species.list | while read line
do 
  mkdir -p /mnt/xuruizhi/brain/fastq/$line
  cd /mnt/xuruizhi/brain/sequence/$line
  fqdir=/mnt/xuruizhi/brain/fastq/$line
  ls *.sra | while read id
  do
    echo "fastq-dump --gzip --split-3 -O ${fqdir} ${id}"
  done >${fqdir}/sra2fq.sh
  nohup sh sra2fq.sh > sra2fq.log & 
done 



# 超算中处理
rsync -av /mnt/xuruizhi/brain/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/ # copy的时候目录没有设置好
rsync -av /mnt/d/biosoft \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/biosoft
## 正式代码
# 三个物种批处理
mkdir -p ~/xuruizhi/brain/barin/fastq
cd ~/xuruizhi/brain/barin/fastq
cat >species.list <<EOF
mouse
rat
human
EOF
cat species.list | while read line
do 
  mkdir -p ~/xuruizhi/brain/brain/fastq/$line
  cd ~/xuruizhi/brain/brain/sequence/$line
  fqdir=~/xuruizhi/brain/brain/fastq/$line
  ls *.sra | while read id
  do
    echo "/share/home/wangq/bin/fastq-dump.3.0.0 --gzip --split-3 -O ${fqdir} ${id}"
  done >${fqdir}/sra2fq.sh
done  



# mouse
# bsub -q mpi -n 24 -J sra2fq_mouse -o ~/xuruizhi/brain/brain/fastq/mouse "sh sra2fq.sh > sra2fq.log"
cat  MOUSE.list | perl -ne ' \
chomp;
print "$_" . ".sra\n";
'
cd ~/xuruizhi/brain/brain/sequence
bsub -q mpi -n 24 -J sra2fq_mouse -o ~/xuruizhi/brain/brain/fastq/mouse " cat MOUSE.list |  parallel --no-run-if-empty --linebuffer -k -j 10 '/share/home/wangq/bin/fastq-dump.3.0.0 --gzip --split-3 -O ~/xuruizhi/brain/brain/fastq/mouse ~/xuruizhi/brain/brain/sequence/{}'"




# rat
cd ~/xuruizhi/brain/brain/fastq/rat
bsub -q mpi -n 24 -J sra2fq_rat -o ~/xuruizhi/brain/brain/fastq/rat "sh sra2fq.sh > sra2fq.log"
# human
cd ~/xuruizhi/brain/brain/fastq/human
bsub -q mpi -n 24 -J sra2fq_human -o ~/xuruizhi/brain/brain/fastq/human "sh sra2fq.sh > sra2fq.log"


# 大循环超算，未尝试
cd ~/xuruizhi/brain/fastq
cat >species.list <<EOF
mouse
rat
human
EOF
cat species.list | while read line
do 
  mkdir -p ~/xuruizhi/brain/fastq/$line
  cd ~/xuruizhi/brain/sequence/$line
  fqdir=~/xuruizhi/brain/fastq/$line
  ls *.sra | while read id
  do
    echo "fastq-dump --gzip --split-3 -O ${fqdir} ${id}"
  done >sra2fq.sh
  bsub -q mpi -n 24 -J sra2fq -o ~/xuruizhi/brain/fastq/$line "sh sra2fq.sh > sra2fq.log"
done
```
# 3. 


# 4. Pre-alinment

## 4.1 quality control checking  

1. 目的：whether the sequencing quality is qualified or not    

2. 使用软件：`FastQC`,  `multiqc`  
FastQC可用于可视化测序数据中的`碱基质量评分`、`GC含量`、序列长度分布、序列重复水平、k-mer的过度表达，及`引物、接头的污染`。  
3. 代码
```bash
# 新建目录  
mkdir /mnt/d/ATAC/fastqc

# ！注意！一定在存储fastq.gz的文件夹路径下执行下面的命令
cd ~/data/sra
fastqc -t 4 -o /mnt/d/ATAC/fastqc/ *.gz

cd /mnt/d/ATAC/fastqc/
multiqc .
```
4. 结果解读：分析看[该文章](https://github.com/outcastaaa/bioinformatics-learning/blob/main/RNA-seq/Tools/fastqc.md)  

可看到，四个样本拆分得到的8个测序数据：   

`sequence quality`都很高，前期因为使用默认参数进行base calling, 所以这部分碱基的质量一般有些较低也不影响，下一步trim掉即可；  
`Per base sequence content`前期波动，后期都很平稳。使用转座酶片段化的文库在读取开始位置都有内在偏差。这种偏差与绝对序列无关，而是在读数的 5' 末端提供了许多不同 K-mer 的富集。 虽然这是一个真正的技术偏差，但它不可以通过trimming来纠正的，并且在大多数情况下不会对下游分析产生不利影响。  
`Per sequence GC content`和上面原因一样，建库偏差导致GC含量和预测不同。
`Sequence Duplication Levels`更specific的子集富集或低复杂性污染物在图的右侧产生尖峰，蓝线中存在红线中消失，可以被去重。  
`Adapter Content`报错是因为含有ATAC-seq转座酶序列，后续去接头即可。  






## 4.2 pre-alinment QC
1. 目的：adapters and low quality reads trimming
2. 使用软件：`Trim Galore`，因为在最后一部卡住，选择使用`cutadapt + trimmomatic`分步修剪  [参考1](https://www.jianshu.com/p/4ee2f4d2292f)  

* Trim Galore可以自动检测接头序列，质控和去除接头两个步骤一起,适用于多种组学去接头.  
3. 代码：  

```bash
mkdir -p /mnt/d/ATAC/trim/
# 构建循环
cd ~/data/sra
ls ./*_1.fastq.gz > ./1
ls ./*_2.fastq.gz > ./2

paste 1 1 2 >config.raw

# 执行trim代码，有时候会卡住，要有耐心
cat config.raw | while read id;
do echo $id 
 arr=($id)
 fq1=${arr[1]}
 fq2=${arr[2]}
 sample=${arr[0]}

trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o /mnt/d/ATAC/trim/  $fq1 $fq2 &
done

#-q 质量；--length 去除长度小于35的reads；-e 允许的最大误差；--paired 双端测序；-o 输出目录；后接 fastq_file1和file2

# 再次质控
fastqc -t 4 -o /mnt/d/ATAC/fastqc_again/ /mnt/d/ATAC/trim2/*.gz
cd /mnt/d/ATAC/fastqc_again/
multiqc .
```
* 或者采用cutadapt + trimmomatic
```bash
# 去接头
mkdir -p /mnt/d/ATAC/cutadapt/

# 构建循环
cd ~/data/sra
ls *_1.fastq.gz > 1
ls *_2.fastq.gz > 2

paste 1 1 2 >config.raw

# 执行cutadapt代码
cat config.raw | while read id;
do echo $id 
 arr=($id)
 fq1=${arr[1]}
 fq2=${arr[2]}
 sample=${arr[0]}

cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -j 6\
    --minimum-length 30 --overlap 3 -e 0.1 --trim-n \
    -o /mnt/d/ATAC/cutadapt/$fq1 -p /mnt/d/ATAC/cutadapt/$fq2 $fq1 $fq2
done
    # --minimum-length 如果剔除接头后read长度低于30，这条read将会被丢弃
    # --overlap        如果两端的序列与接头有4个碱基的匹配将会被剔除
    # --trim-n         剔除两端的N
    # -a 去除3端引物序列
    # -e 容错率，默认为0.1
	# --discard-untrimmed 去除没有adapter的reads，不要乱用


# 去低质量reads
mkdir -p /mnt/d/ATAC/trimmomatic/paired
mkdir -p /mnt/d/ATAC/trimmomatic/unpaired
cat config.raw | while read id;
do echo $id 
 arr=($id)
 fq1=${arr[1]}
 fq2=${arr[2]}
 sample=${arr[0]}

 # Trimmomatic-0.38记得更改路径
    java -jar /mnt/d/biosoft/Trimmomatic-0.38/Trimmomatic-0.38.jar \
    PE -threads 4 -phred33 /mnt/d/ATAC/cutadapt/$fq1 /mnt/d/ATAC/cutadapt/$fq2 \
	/mnt/d/ATAC/trimmomatic/paired/$fq1 /mnt/d/ATAC/trimmomatic/unpaired/$fq1 \
	/mnt/d/ATAC/trimmomatic/paired/$fq2 /mnt/d/ATAC/trimmomatic/unpaired/$fq2 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 
done

  # LEADING:20，从序列的开头开始去掉质量值小于 20 的碱基
  # TRAILING:20，从序列的末尾开始去掉质量值小于 20 的碱基
  # SLIDINGWINDOW:5:15，从 5' 端开始以 5bp 的窗口计算碱基平均质量，如果此平均值低于 15，则从这个位置截断read
  # MINLEN:30， 如果 reads 长度小于 30bp 则扔掉整条 read。


# 再次质控
fastqc -t 4 -o /mnt/d/ATAC/fastqc_again/again  /mnt/d/ATAC/trimmomatic/paired/*.gz
cd /mnt/d/ATAC/fastqc_again/again
multiqc . 
```
 
4. 结果：  
* 储存在[trim_galore](https://github.com/outcastaaa/ATAC/tree/main/trim_galore)和[cutadapt+Trimmomatic](https://github.com/outcastaaa/ATAC/tree/main/cutadapt%2BTrimmomatic)  

* 再次质控结果  
其他没有合格的板块或不影响下游分析，或可以通过后续步骤解决，这一步成功去除了`Adapter Content`版块的接头序列。





























# 参考文章  


1. 恒河猴初级运动皮层，预印本

ATAC-seq filtering: Duplicate reads, mitochondrial reads, and multi-mapping reads were filtered out as a part of the ENCODE ATAC-seq pipeline. We made the following read filtering changes from the pipeline defaults: "atac.multimapping" : 0.
ATAC-seq peak calling: Peaks were called on each sample individually using MACS2 within the ENCODE ATAC-seq pipeline workflow. We made the following changes from pipeline defaults: "atac.cap_num_peak" : 300000 and "atac.smooth_win" : 150.
ATAC-seq peak filtering: Reproducible peaks among replicates were determined within the ENCODE ATAC-seq pipeline workflow using IDR with the following changes from pipeline defaults: "atac.enable_idr" : true and "atac.idr_thresh" : 0.1.
Other: More information on all data processing and conclusions can be found in the publication.
Genome_build: mm10, rn6


2. Nat Neurosci 2017 Mar	神经元活动改变了成年大脑中的染色质可及性景观

The transposome assay was performed as previously described17. Around 50,000 nuclei from each dentate gyrus were used. qPCR was used to estimate of the number of additional cycles needed to generate products at 25% saturation. Typically, four to seven additional PCR cycles were added to the initial set of five cycles. The library was purified on AMPure XP beads (Beckman A63881) and analyzed on an Agilent Bioanalyzer, and 50-bp paired-end sequencing was performed using an Illumina HiSeq 2500 platform according to standard operating procedures. Sequencing reads were mapped to a mouse genome assembly (mm9) from the UCSC genome browser (http://genome.ucsc.edu/) using Bowtie2 (ref. 46). Duplicate reads were marked and removed by PICARD Tools (http://broadinstitute.github.io/picard/). Open chromatin peaks were analyzed using MACS (ref. 47) software. Differential peaks between different conditions were generated by diffReps (ref. 48). Significant overlap between two sets of genomic regions was tested using GAT (ref. 49). HOMER25 was used to annotate those peaks for motif discovery analysis. IGV (ref. 50) and ngsplot (ref. 51) were used for visualization of raw intensities. Statistical analyses were performed using an in-house R script unless otherwise specified.  

3. Mol Psychiatry 2022 Nov	揭示酒精诱导的抗焦虑分解过程中的表观基因组和转录组相互作用

使用BWA MEM（Burrows-Wheeler Aligner）（v6.0.0）将读取与Rattus norvegicus Rnor_7.15基因组进行比对，并调整比对以考虑转座子结合：+链比对+4bp，-链比对-5bp.使用Macs2（v2.1.1）执行峰值调用。在分数阈值>5处过滤峰。在差异分析之前，使bedtools（v2.27.1）合并将峰合并到并集中，并使用featureCounts（v1.5.0）根据读取对齐将峰丰度量化为原始计数。执行归一化并读取表示为 CPM（每百万计数），此外，使用 edgeR 计算 TMM 归一化和差异统计。我们使用错误发现率（FDR）校正调整了多个测试的p值。  

4. Elife 2021 Dec	5-羟甲基胞嘧啶介导的活性去甲基化是哺乳动物神经元分化和功能所必需的

ATAC-Seq 数据：读取使用参数“--严格 3 --fastqc --paired”的trim_galore进行处理。使用带有参数“-X 10 --no-mixed --no-discordant”的 bowtie2（版本 2.1.0）将修剪后的读取映射到 mm2000。使用samtools删除了重复项。使用deepTools bamCompare模块将读取规范化为RPKM，忽略chrX，chrY，chrM并过滤读取，最低映射质量为30。Metagene和热图配置文件是使用deepTools模块computeMatrix，plotProfile和plotHeatmap生成的。100 nt以下的唯一片段用于调用带有参数“--nomodel -q 2.0 --call-summits”的macs01峰值。宽峰也用macs2以及参数“--nomodel -f BAM --keep-dup all --broad -g mm -B -q 0.01”调用。对峰进行滤波，以去除 chrY 和 chrM 峰以及与 ENCODE 中的 mm10 黑名单重叠的峰.  

5. Nucleic Acids Res 2020 Sep	增强子RNA预测增强子-基因调控链接，对神经元系统中的增强子功能至关重要

Low-quality bases (Phred <20) and Nextera adapters (5’-CTGTCTCTTATA-3’) were identified and trimmed using FastQC and TrimGalore (v0.4.5).FASTQ files containing trimmed sequences were then aligned to Rn6 Ensembl genome assembly (v95) to generate binary alignment map (BAM) files with Bowtie2 (v.2.3.4.2) with custom options: `--local --very-sensitive-local`, and `-X 3000`。Post-trimming QC with FastQC。Sequences ordered by genomic position using Samtools (v.1.9)。BAM files for each brain region were merged to generate three metasamples that were used for downstream data analysis. Peaks for each region were called using macs2 (Zhang et al., 2008) callpeak with the options - -qvalue 0.00001 - -gsize 2729862089 - -format BAMPE, ignoring duplicates.Within each brain region, peaks closer than 1000bps were merged with bedtools and peaks less than 146bp, the specific length of DNA wrapped around a single nucleosome, were removed in R. Finally, peaks from each brain region were merged with bedtools to create an additive peak set containing 192,830 peaks.Genome_build: rn6。Supplementary_files_format_and_content: BED file cotaining ATAC-seq peaks (output of MACS2). Data are in BED format - columns are chromosome, start, end, and score for how many cell types (of 3) contained the ATAC peak.  

6. 
