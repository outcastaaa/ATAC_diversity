- [1. 建目录](#1-建目录)
- [2. 下载数据](#2-下载数据)
- [3. 比对前质控](#3-比对前质控)
- [4. 比对](#4-比对)
- [5. 比对后过滤](#5-post-alignment-processing)
- [6. Tn5转换](#6-shift-reads)
- [7. call peak](#7-call-peaks)
- [8. 可视化](#8-visualization)
- [9. diffbind](#9-使用diffbind做主成分分析)


# 1. 建目录
```bash
cd /mnt/xuruizhi # 挂载到NAS
mkdir -p /mnt/xuruizhi/brain/sequence/mouse
mkdir -p /mnt/xuruizhi/brain/sequence/rat
mkdir -p /mnt/xuruizhi/brain/sequence/macaque
mkdir -p /mnt/xuruizhi/brain/sequence/human
```

# 2. 下载数据
1. 测序数据
```bash
# mouse 21个，后三个是单端测序
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
SRR13443549
SRR13443553
SRR13443554
EOF

prefetch SRR11179779  SRR13049359  SRR13049362  SRR14362275  SRR14362282  SRR3595212 \
SRR11179780  SRR13049360  SRR13049363  SRR14362276  SRR14614715  SRR3595213 \
SRR11179781  SRR13049361  SRR13049364  SRR14362281  SRR3595211   SRR3595214 \
SRR13443553 SRR13443549 SRR13443554
cd ~/data/sra
cp *.sra /mnt/xuruizhi/brain/sequence/mouse
rm *.sra



# rat 12个,全是双端测序
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


# human 本数据库提供的为bam文件，下载了暂且不使用
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

# new_human

```
2. 基因组数据
```bash
hg19：
    $ cd /your/path/of/reference/
    $ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    $ tar zvfx chromFa.tar.gz
    $ cat *.fa > hg19.fa
    $ rm chr*.fa
建立bwa索引：
    $ bwa index -a bwtsw  hg19.fa
    # 产生.bwt .pac .ann .amb .sa五个新文件
    # -a：两种构建index算法，bwtsw以及is，bwtsw适用大于10MB的参考基因组，比如人，is适用于小于2GB的数据库，是默认的算法，速度较快，需要较大的内存，
    # -p：输出数据库的前缀，默认与输入文件名一致，这里我没有加这个参数，直接输出到当前目录
建立bowtie2索引：
    $ bowtie2-build hg19.fa hg19.fa
    #  bowtie2-build命令在安装bowtie2的目录下找到
    # 第一个hg19.fa代表输入的参考序列
    # 第二个hg19.fa代表输出的索引文件前缀
    # 产生六个.bt2新文件
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
  bsub -q mpi -n 24 -J sra2fq -o ~/xuruizhi/brain/fastq/$line "bash sra2fq.sh > sra2fq.log"
done
```
# 3. 比对前质控

```bash
# 由于相邻转座的最小间隔为 38 bp，通常 38 bp 以下的片段直接删除
rsync -av /mnt/xuruizhi/brain/fastq/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/brain/fastq/

# 新建目录  
mkdir -p ~/xuruizhi/brain/brain/fastqc/mouse

# ！注意！一定在存储fastq.gz的文件夹路径下执行下面的命令
cd ~/xuruizhi/brain/brain/fastq/mouse/
bsub -q mpi -n 24 -J fastqc -o ~/xuruizhi/brain/brain/fastqc/mouse " \
fastqc -t 4 -o ~/xuruizhi/brain/brain/fastqc/mouse/ *.gz"
# job <8592824>
```

* 超算达咩  

```bash
# 超算下载安装包，发现搞不来在本地跑吧
cd /mnt/d/biosoft
wget -c https://files.pythonhosted.org/packages/a3/30/4a889a6916d7480c153774777e634b89865f95cb02f2c3209762c7ef984b/cutadapt-4.1.tar.gz
tar -zxvf cutadapt-4.1.tar.gz
cd cutadapt-4.1
python setup.py install --user

mkdir -p ~/xuruizhi/brain/brain/trim/mouse

# 双端测序
# 构建循环
cd ~/xuruizhi/brain/brain/fastq/mouse/
ls ./*_1.fastq.gz > ./1
ls ./*_2.fastq.gz > ./2
paste 1 1 2 >config.raw
# 执行trim代码，有时候会卡住，要有耐心
bsub -q mpi -n 24 -J fastqc -o ~/xuruizhi/brain/brain/trim/mouse " \
cat config.raw | while read id;
do echo $id 
 arr=($id)
 fq1=${arr[1]}
 fq2=${arr[2]}
 sample=${arr[0]}

~/xuruizhi/biosoft/biosoft/TrimGalore-0.6.6/trim_galore \
--phred33 --length 35 -e 0.1 --stringency 3 --paired -o ~/xuruizhi/brain/brain/trim/mouse  $fq1 $fq2  --path_to_cutadapt ~/xuruizhi/biosoft/cutadapt-4.1
done "

#-q 质量；--length 去除长度小于35的reads；-e 允许的最大误差；--paired 双端测序；-o 输出目录；后接 fastq_file1和file2

# 再次质控
mkdir -p ~/xuruizhi/brain/brain/fastqc_again/mouse/
fastqc -t 4 -o ~/xuruizhi/brain/brain/fastqc_again/mouse/ ~/xuruizhi/brain/brain/trim/mouse/*.gz
cd  ~/xuruizhi/brain/brain/fastqc_again/mouse/
multiqc .
```
* 本地质控  
```bash
# 新建目录  
mkdir -p /mnt/xuruizhi/brain/fastqc/
# 通过beyond compare传输到/mnt/d/brain/brain/fastqc/mouse
cd /mnt/d/brain/brain/fastqc/mouse
multiqc .
cp /mnt/d/brain/brain/fastqc/mouse/* /mnt/xuruizhi/brain/fastqc/mouse/


mkdir -p /mnt/xuruizhi/brain/trim/mouse
# 双端测序
# 构建循环
cd /mnt/xuruizhi/brain/fastq/mouse/
ls *_1.fastq.gz > 1
sed 's/_1.fastq.gz//g' ./1 >name.list
# 执行trim代码，有时候会卡住，要有耐心
cat name.list | parallel --no-run-if-empty --linebuffer -k -j 8 ' \
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired \
-o /mnt/xuruizhi/brain/trim/mouse  ./{}_1.fastq.gz ./{}_2.fastq.gz'
#-q 质量；--length 去除长度小于35的reads；-e 允许的最大误差；--paired 双端测序；-o 输出目录；后接 fastq_file1和file2

# 单端测序
cd /mnt/xuruizhi/brain/fastq/mouse
cat >single.list <<EOF
SRR13443549.fastq.gz
SRR13443553.fastq.gz
SRR13443554.fastq.gz
EOF
cat single.list | while read id
do 
  trim_galore --phred33 --length 35 -e 0.1 --stringency 3 \
  -o /mnt/xuruizhi/brain/trim/mouse ${id}
done

# 再次质控
mkdir -p /mnt/xuruizhi/brain/fastqc_again/mouse/
fastqc -t 6 -o /mnt/xuruizhi/brain/fastqc_again/mouse/ /mnt/xuruizhi/brain/trim/mouse/*.gz
cd  /mnt/xuruizhi/brain/fastqc_again/mouse/
multiqc .

# 传到超算
mkdir -p /share/home/wangq/xuruizhi/brain/brain/trim/
rsync -av /mnt/xuruizhi/brain/trim/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/brain/trim/
mkdir -p /share/home/wangq/xuruizhi/brain/brain/fastqc_again/mouse/
rsync -av /mnt/xuruizhi/brain/fastqc_again/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/brain/fastqc_again/
```
# 4. 比对
## 4.1 比对
1. 参考基因组  
```bash
# 小鼠 mm10 的 bowtie2 的 index 已经建立过，在超算中处理
mkdir -p /mnt/xuruizhi/brain/genome/mouse/
cp /mnt/d/atac/genome/* /mnt/xuruizhi/brain/genome/mouse/
# 进入超算
mkdir -p /share/home/wangq/xuruizhi/brain/brain/genome/mouse/
rsync -av /mnt/d/atac/genome/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/brain/genome/mouse/
```
2. 比对，在超算
```bash
mkdir -p ~/xuruizhi/brain/brain/alignment/mouse
# 循环 
cd ~/xuruizhi/brain/brain/fastq/mouse
ls *_1.fastq.gz > 1
sed 's/_1.fastq.gz//g' ./1 >name.list
cp ~/xuruizhi/brain/brain/fastq/mouse/name.list ~/xuruizhi/brain/brain/trim/mouse/pair.list

cd ~/xuruizhi/brain/brain/trim/mouse/
ls *_1_val_1.fq.gz > ./1
ls *_2_val_2.fq.gz > ./2
paste pair.list 1 2  > name.list
cat name.list
# SRR11179779
# SRR11179780
# SRR11179781
# SRR13049359
# SRR13049360
# SRR13049361
# SRR13049362
# SRR13049363
# SRR13049364
# SRR14362275
# SRR14362276
# SRR14362281
# SRR14362282
# SRR14614715
# SRR3595211
# SRR3595212
# SRR3595213
# SRR3595214

cd ~/xuruizhi/brain/brain/trim/mouse/
bsub -q mpi -n 24 -J align -o ~/xuruizhi/brain/brain/alignment/mouse " 
cat ~/xuruizhi/brain/brain/trim/mouse/pair.list | \
parallel --no-run-if-empty --linebuffer -k -j 4 ' bowtie2  -p 6  \
-x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz \
2> ~/xuruizhi/brain/brain/alignment/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment/mouse/{}.sam'"
# Job <8595215>
# better alignment results are frequently achieved with --very-sensitive
# use -X 2000 to allow larger fragment size (default is 500)

# 单端测序
cd ~/xuruizhi/brain/brain/trim/mouse/
bowtie2_index=~/xuruizhi/brain/brain/genome/mouse
align_dir=~/xuruizhi/brain/brain/alignment/mouse

ls *_trimmed.fq.gz > single
sed 's/_trimmed.fq.gz//g' single >single.list

bsub -q mpi -n 24 -J align_single -o ~/xuruizhi/brain/brain/alignment/mouse " 
cat ~/xuruizhi/brain/brain/trim/mouse/single.list | \
parallel --no-run-if-empty --linebuffer -k -j 4 ' bowtie2  -p 1  \
-x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -U {}_trimmed.fq.gz \
2> ~/xuruizhi/brain/brain/alignment/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment/mouse/{}.sam '"
# Job <8595222>




# 超级大批量处理，修雷教的
cd ~/xuruizhi/brain/brain/trim/mouse/ 
mkdir -p ~/xuruizhi/brain/brain/alignment_new/mouse
# 写双端批量
cat 111.sh
#!/usr/bin bash
bowtie2  -p 48 -x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz \
2> ~/xuruizhi/brain/brain/alignment_new/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment_new/mouse/{}.sam

cat pair.list  | while read id; do sed "s/{}/${id}/g" 111.sh > ${id}.sh; done
# Job <8601500> is submitted to queue <mpi>.
# Job <8601501> is submitted to queue <mpi>.
# Job <8601502> is submitted to queue <mpi>.
# Job <8601503> is submitted to queue <mpi>.
# Job <8601504> is submitted to queue <mpi>.
# Job <8601505> is submitted to queue <mpi>.
# Job <8601506> is submitted to queue <mpi>.
# Job <8601507> is submitted to queue <mpi>.
# Job <8601508> is submitted to queue <mpi>.
# Job <8601509> is submitted to queue <mpi>.
# Job <8601510> is submitted to queue <mpi>.
# Job <8601511> is submitted to queue <mpi>.
# Job <8601512> is submitted to queue <mpi>.
# Job <8601513> is submitted to queue <mpi>.
# Job <8601514> is submitted to queue <mpi>.
# Job <8601515> is submitted to queue <mpi>.
# Job <8601516> is submitted to queue <mpi>.


# 写双端批量
cat 222.sh
#!/usr/bin bash
bowtie2  -p 48 -x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -U {}_trimmed.fq.gz \
2> ~/xuruizhi/brain/brain/alignment_new/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment_new/mouse/{}.sam

cat single.list  | while read id; do sed "s/{}/${id}/g" 222.sh > ${id}.sh; done

cat single.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/alignment_new/mouse/ "bash ${id}.sh"
done
# Job <8602391> is submitted to queue <mpi>.
# Job <8602392> is submitted to queue <mpi>.
# Job <8602393> is submitted to queue <mpi>.
```
* SRR14614715.sam太大了，比对已经完成，但是后续步骤先不进行

# 除了SRR14614715.sam重新比对,在~/xuruizhi/brain/brain/alignment_new/mouse/

## 4.2 sort_transfertobam_index  
```bash
mkdir -p ~/xuruizhi/brain/brain/sort_bam/mouse
cd ~/xuruizhi/brain/brain/alignment/mouse/

cat >name.list <<EOF
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
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF


bamdir=~/xuruizhi/brain/brain/sort_bam/mouse
# sam to bam
bsub -q mpi -n 24 -J sort_index -o $bamdir " 
  cat name.list | parallel --no-run-if-empty --linebuffer -k -j 4 '
  samtools sort -@ 6 {}.sam > $bamdir/{}.sort.bam
  samtools index -@ 6 $bamdir/{}.sort.bam
  samtools flagstat  -@ 6 $bamdir/{}.sort.bam > $bamdir/{}.raw.stat'"
# job <8600760>

# samtools index为已经基于坐标排序后bam或者cram的文件创建索引，默认在当前文件夹产生*.bai的index文件
# raw.stat记录匹配后原始文件情况
```
## 大批量转化
```bash
mkdir -p h
cd ~/xuruizhi/brain/brain/alignment_new/mouse/
cat >name.list <<EOF
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
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF


cat samtobam.sh
#!/usr/bin bash
samtools sort -@ 48 {}.sam > ~/xuruizhi/brain/brain/sort_bam_new/mouse/{}.sort.bam
samtools index -@ 48 ~/xuruizhi/brain/brain/sort_bam_new/mouse/{}.sort.bam
samtools flagstat  -@ 48 ~/xuruizhi/brain/brain/sort_bam_new/mouse/{}.sort.bam > ~/xuruizhi/brain/brain/sort_bam_new/mouse/{}.raw.stat

cat name.list  | while read id; do sed "s/{}/${id}/g" samtobam.sh > ${id}.sh; done

cat name.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/sort_bam_new/mouse/ "bash ${id}.sh"
done
```

# 5. Post-alignment processing 
1. 目的：  

去除没有匹配到的、匹配得分较低的、重复的reads(如果两条reads具有相同的长度而且比对到了基因组的同一位置，那么就认为这样的reads是由PCR扩增而来)；去除线粒体中染色质可及区域及ENCODE blacklisted regions。    

## 5.1 remove PCR-duplicate reads
目的：去除因为PCR偏好性导致的reads重复扩增  

```bash
mkdir -p ~/xuruizhi/brain/brain/rmdup_new/mouse
cd ~/xuruizhi/brain/brain/sort_bam_new/mouse

cat >name.list <<EOF
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
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
EOF



rmdup_dir=~/xuruizhi/brain/brain/rmdup/mouse
cd ~/xuruizhi/brain/brain/sort_bam/mouse
bsub -q mpi -n 24 -J rmdup -o $rmdup_dir " 
cat name.list | parallel -k -j 8 '
picard MarkDuplicates -I {}.sort.bam \
	 -O $rmdup_dir/{}.rmdup.bam \
	 -REMOVE_DUPLICATES true \
   -VALIDATION_STRINGENCY LENIENT \
   -METRICS_FILE $rmdup_dir/{}.log
samtools index -@ 3 $rmdup_dir/{}.rmdup.bam
samtools flagstat -@ 3 $rmdup_dir/{}.rmdup.bam > $rmdup_dir/{}.rmdup.stat'"
# job <8600145>

picard MarkDuplicates -I  SRR13049360.sort.bam \
	 -O $rmdup_dir/SRR13049360.rmdup.bam \
	 -REMOVE_DUPLICATES true \
   -VALIDATION_STRINGENCY LENIENT \
   -METRICS_FILE $rmdup_dir/SRR13049360.log
samtools index -@ 3 $rmdup_dir/SRR13049360.rmdup.bam
samtools flagstat -@ 3 $rmdup_dir/SRR13049360.rmdup.bam > $rmdup_dir/SRR13049360.rmdup.stat
#--VALIDATION_STRINGENCY <验证严格性>此程序读取的所有 SAM 文件的验证严格性。
#将严格性设置为 SILENT 可以提高处理 BAM 文件时的性能，其中可变长度数据（读取、质量、标签）不需要解码。
#默认值：严格。 可能的值：{STRICT、LENIENT、SILENT}
```
## 5.1 大批量去重
```bash
mkdir -p ~/xuruizhi/brain/brain/rmdup_new/mouse
cd ~/xuruizhi/brain/brain/sort_bam_new/mouse

cat >name.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049360.
SRR13049361
SRR13049362..
SRR13049363..
SRR13049364..
SRR14362275.
SRR14362276
SRR14362281
SRR14362282
SRR3595211..
SRR3595212..
SRR3595213
SRR3595214..
SRR13443549
SRR13443553
SRR13443554
EOF



cd ~/xuruizhi/brain/brain/sort_bam_new/mouse
cat rmdup.sh
#!/usr/bin bash
parallel -k -j 24 'picard MarkDuplicates -I {}.sort.bam -O ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam \
-REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT \
-METRICS_FILE ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.log'
cat name.list  | while read id; do sed "s/{}/${id}/g" rmdup.sh > ${id}.sh; done

cat name.list | while read id
do
  bsub -q mpi -n 24 -o ~/xuruizhi/brain/brain/rmdup_new/mouse "bash ${id}.sh"
done
# 或者
rmdup_dir=~/xuruizhi/brain/brain/rmdup_new/mouse
cd ~/xuruizhi/brain/brain/sort_bam_new/mouse
bsub -q largemem -n 24 -J rmdup -o $rmdup_dir " 
cat name.list | parallel -k -j 12 '
picard MarkDuplicates -I {}.sort.bam \
	 -O $rmdup_dir/{}.rmdup.bam \
	 -REMOVE_DUPLICATES true \
   -VALIDATION_STRINGENCY LENIENT \
   -METRICS_FILE $rmdup_dir/{}.log'"


# picard老是出问题，先跑成功的

cd ~/xuruizhi/brain/brain/rmdup_new/mouse
cat >name_new.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049361
SRR14362276
SRR14362281
SRR14362282
SRR3595213
SRR13443549
SRR13443553
SRR13443554
EOF

cat >index.sh <<EOF
#!/usr/bin bash

samtools index -@ 48 ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam
samtools flagstat -@ 48 ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam > ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.stat
EOF

cat name_new.list  | while read id; do sed "s/{}/${id}/g" index.sh > ${id}.sh; done

cat name_new.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/rmdup_new/mouse "bash ${id}.sh"
done
# Job <8603828-39> is submitted to queue <mpi>.
```


## 5.2 remove bad quality reads
* 目的：保留都比对到同一个染色体的paired reads（proper paired），同时质量较高的reads (mapping quality>=30) 

```bash
samtools view -f 2 -q 30 -o test.filter.bam test.rmdup.bam
# -f Retain properly paired reads -f 2
# -q 取mapping质量大于30的reads
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804) 看情况取舍
```

## 5.3 remove chrM reads
* 目的：去除比对到线粒体上的reads，这一步一定要做，线粒体上长度小，极大概率覆盖很多reads，造成虚假peak。由于mtDNA读段的百分比是文库质量的指标，我们通常在比对后删除线粒体读段。  

* 统计chrM reads，使用没有去除PCR重复，先不做了
```bash
mkdir -p ~/xuruizhi/brain/brain/stat/mouse
cp ~/xuruizhi/brain/brain/sort_bam/mouse/*.raw.stat ~/xuruizhi/brain/brain/stat/mouse/
cp ~/xuruizhi/brain/brain/sort_bam/mouse/name.list ~/xuruizhi/brain/brain/stat/mouse/

cd ~/xuruizhi/brain/brain/sort_bam/mouse

bsub -q mpi -n 24 -J stat -o ~/xuruizhi/brain/brain/stat/mouse " 
cat name.list | parallel --no-run-if-empty --linebuffer -k -j 12 '
  samtools idxstats {}.sort.bam | grep 'chrM' | cut -f 3  
  samtools idxstats {}.sort.bam | awk '{SUM += $3} END {print SUM}''"

# 第一列是染色体名称，第二列是序列长度，第三列是mapped reads数，第四列是unmapped reads数
```


* 将上一步和这一步结合起来
```bash
mkdir -p ~/xuruizhi/brain/brain/filter_new/mouse
# ~/xuruizhi/brain/brain/filter/mouse
cp ~/xuruizhi/brain/brain/rmdup_new/mouse/name_new.list ~/xuruizhi/brain/brain/filter_new/mouse/

cd ~/xuruizhi/brain/brain/rmdup_new/mouse
cat >filter.sh <<EOF
#!/usr/bin bash

samtools view -h -f 2 -F 1804 -q 30 {}.rmdup.bam | grep -v  chrM | samtools sort -@ 6 -O bam  -o ~/xuruizhi/brain/brain/filter_new/mouse/{}.filter.bam
samtools index -@ 48 ~/xuruizhi/brain/brain/filter_new/mouse/{}.filter.bam
samtools flagstat -@ 48 ~/xuruizhi/brain/brain/filter_new/mouse/{}.filter.bam > ~/xuruizhi/brain/brain/filter_new/mouse/{}.filter.stat
EOF
# -F 1804包含read unmapped， mate unmapped，not primary alignment，read fails platform/vendor quality checks，  read is PCR or optical duplicate适合双端测序
cat name_new.list  | while read id; do sed "s/{}/${id}/g" filter.sh > ${id}_filter.sh; done

cat name_new.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/filter_new/mouse "bash ${id}_filter.sh"
done
# Job <8604250-263> is submitted to queue <mpi>.
```

## 5.4 Blacklist filtering

1. 目的：去除ENCODE blacklisted 区域，通过blacklist的过滤，可以进一步降低peak calling的假阳性。    

本流程采用的方法是：在peak calling之前去除，比对后的reads 去除PCR重复等后单独去除 blacklist region，再 call peak.  
   

```bash
# 前面filter删除的太狠了,SRR13443549 SRR13443553 SRR13443554
# 在本地跑
# 下载对应物种的 blacklist.bed文件
cd /mnt/d/ATAC/blklist
wget https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gzip -dc mm10.blacklist.bed.gz > mm10.blacklist.bed
rm *.gz
wc -l  mm10.blacklist.bed #164

mkdir -p /mnt/xuruizhi/brain/blklist/mouse
cp /mnt/d/ATAC/blklist/mm10.blacklist.bed /mnt/xuruizhi/brain/blklist/mouse

# 通过beyond compare，将~/xuruizhi/brain/brain/filter_new/mouse传输到/mnt/d/brain/brain/filter_new/mouse
mkdir -p /mnt/xuruizhi/brain/filter_new/mouse
cp /mnt/d/brain/brain/filter_new/mouse/* /mnt/xuruizhi/brain/filter_new/mouse


cd /mnt/xuruizhi/brain/filter_new/mouse
cat >name_new.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049361
SRR14362276
SRR14362281
SRR14362282
SRR3595213
SRR13443549
SRR13443553
SRR13443554
EOF

cat name_new.list | while read id;
do 
  echo $id 
  echo "${id}.filter.bam"

  # 取交集看bam文件和blacklist有多少重合部分
  bedtools intersect -wa -a ${id}.filter.bam \
  -b /mnt/xuruizhi/brain/blklist/mouse/mm10.blacklist.bed | \
  wc -l  > /mnt/xuruizhi/brain/blklist/mouse/${id}.intersect.list
done

mkdir -p /mnt/xuruizhi/brain/final/mouse
cd /mnt/xuruizhi/brain/filter_new/mouse
cat name_new.list | while read id;
do 
  echo "$id"
  # 凡是bam中含有blacklist都删除
  bedtools intersect -v -a ${id}.filter.bam -b ../../blklist/mouse/mm10.blacklist.bed > ../../final/mouse/${id}.final.bam
  samtools index -@ 7 ../../final/mouse/${id}.final.bam
  samtools flagstat -@ 7 ../../final/mouse/${id}.final.bam > ../../final/mouse/${id}.final.stat
done
```
6. 结果解读：  
```bash
# 原比对文件数据，以SRR13049361为例
171031266 + 0 in total (QC-passed reads + QC-failed reads)
171031266 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
168538658 + 0 mapped (98.54% : N/A)
168538658 + 0 primary mapped (98.54% : N/A)
171031266 + 0 paired in sequencing
85515633 + 0 read1
85515633 + 0 read2
166233644 + 0 properly paired (97.19% : N/A)
167189182 + 0 with itself and mate mapped
1349476 + 0 singletons (0.79% : N/A)
101336 + 0 with mate mapped to a different chr
9484 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除PCR重复+低质量+chrM后数据
63564228 + 0 in total (QC-passed reads + QC-failed reads)
63564228 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
63564228 + 0 mapped (100.00% : N/A)
63564228 + 0 primary mapped (100.00% : N/A)
63564228 + 0 paired in sequencing
31782114 + 0 read1
31782114 + 0 read2
63564228 + 0 properly paired (100.00% : N/A)
63564228 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除blacklist后数据
63547428 + 0 in total (QC-passed reads + QC-failed reads)
63547428 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
63547428 + 0 mapped (100.00% : N/A)
63547428 + 0 primary mapped (100.00% : N/A)
63547428 + 0 paired in sequencing
31773708 + 0 read1
31773720 + 0 read2
63547428 + 0 properly paired (100.00% : N/A)
63547428 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
到这一步，比对文件已经过滤完成。     


## 5.5 bamtobed
把处理好的bam比对文件转化为bed格式

```bash
# bam to bed
mkdir -p /mnt/xuruizhi/brain/bed/mouse
cd /mnt/xuruizhi/brain/final/mouse

parallel -j 6 "
   bedtools bamtobed -i ./{1} >../../bed/mouse/{1}.bed
" ::: $( ls *.final.bam)


# 此处只进行bed 

```

* 结果：
```bash
# bed
$ cat SRR11179779.final.bam.bed | head -n 10
chr1    3000665 3000748 SRR11179779.29048935/2  42      +
chr1    3001025 3001175 SRR11179779.29048935/1  42      -
chr1    3001182 3001332 SRR11179779.13631674/2  42      +
chr1    3001547 3001697 SRR11179779.13631674/1  42      -
chr1    3001663 3001813 SRR11179779.34885572/1  37      +
chr1    3001709 3001859 SRR11179779.34885572/2  37      -
chr1    3003183 3003326 SRR11179779.15040555/1  42      +
chr1    3003208 3003326 SRR11179779.15040555/2  42      -
chr1    3003705 3003785 SRR11179779.1674968/1   42      +
chr1    3003716 3003785 SRR11179779.1674968/2   42      -
```

```bash
# 必选的三列：
chrom - 染色体的名称（例如chr3，chrY，chr2_random）或支架（例如scaffold10671）。
chromStart- 染色体或scanfold中特征的起始位置。染色体中的第一个碱基编号为0。
chromEnd- 染色体或scanfold中特征的结束位置。所述 chromEnd碱没有包括在特征的显示。\
例如，染色体的前100个碱基定义为chromStart = 0，chromEnd = 100，并跨越编号为0-99的碱基。

# 9个可选的BED字段：
name - 定义BED行的名称。当轨道打开到完全显示模式时，此标签显示在Genome浏览器窗口中BED行的左侧，或者在打包模式下直接显示在项目的左侧。
score - 得分在0到1000之间。如果此注释数据集的轨迹线useScore属性设置为1，则得分值将确定显示此要素的灰度级别（较高的数字=较深的灰色）。
strand - 定义strand。要么“。” （=无绞线）或“+”或“ - ”。
```

# 6. shift reads

由于Tn5酶是以二聚体的形式结合到染色体上的，其跨度大致是9bp，在第一篇ATAC-seq出来的时候，作者就考虑到了这个问题，在分析的时候，需要回补这个9个bp的碱基差。具体做法就是将正链正向移动4bp，将负链负向移动5个bp。一般用alignmentSieve 一步到位。注意，不做reads shift 对单碱基分辨高的分析会有影响，例如TF motif footprinting，但也不是所有TF footprinting分析软件需要shifted reads，很多可以自己转换，e.g. NucleoATAC。   

```bash
mkdir -p /mnt/xuruizhi/brain/Tn5_shift/mouse
cp /mnt/d/ATAC/rmdup/config.raw /mnt/d/ATAC/bedpe/config.raw


# bed转化
cd /mnt/xuruizhi/brain/bed/mouse
cat >name_new.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049361
SRR14362276
SRR14362281
SRR14362282
SRR3595213
SRR13443549
SRR13443553
SRR13443554
EOF

cat name_new.list | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  cat ${sample}.final.bam.bed | awk -v \
  OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4} \
   else if($6=="-"){print $1,$2-5,$3-5}}' \
    > /mnt/xuruizhi/brain/Tn5_shift/mouse/${sample}.Tn5.bed
done
```
4. 结果解读：  


！注意，后续callpeak不可直接使用bedtools转化的bedpe文件，只能包含三行信息：chr,chrom_start,chrom_end
```bash
cd /mnt/xuruizhi/brain/Tn5_shift/mouse/
$ cat SRR13049361.Tn5.bed | head -n 5
chr1    3000777 3000877
chr1    3000779 3000879
chr1    3000797 3000897
chr1    3000877 3000973
chr1    3000922 3001022
$ wc -l SRR13049361.Tn5.bed
# 47997002
```




# 7. Call peaks 
1. 目的： 下一步需要在统计学上判断真实的peak，因为Tn5在染色体上结合是个概率事件，如何判断这个位置的reads足够为一个peak，这就需要用到统计检测。ATAC-seq 数据分析的第二个主要步骤是识别开放区域（也称为 Peak），后续高级分析以此为基础。  

```bash
mkdir -p /mnt/d/ATAC/macs2_peaks/
cd /mnt/d/ATAC/Tn5_shift/


# 如果用的不是专门双端测序的bedpe，而是bed文件，采用下面代码
# 单个样本
mkdir -p /mnt/xuruizhi/brain/macs2_peaks/mouse
cd /mnt/xuruizhi/brain/Tn5_shift/mouse
cat >name_new.list <<EOF
SRR11179779
SRR11179780
SRR11179781
SRR13049359
SRR13049361
SRR14362276
SRR14362281
SRR14362282
SRR3595213
SRR13443549
SRR13443553
SRR13443554
EOF

# 循环
cat name_new.list | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  macs2 callpeak  -g mm --nomodel \
  --shift -100 --extsize 200 -n ${sample} -t ./${sample}.Tn5.bed \
  --outdir /mnt/xuruizhi/brain/macs2_peaks/mouse 
done
mkdir -p /mnt/d/brain/brain/macs2_peaks/mouse
cp /mnt/xuruizhi/brain/macs2_peaks/mouse/* /mnt/d/brain/brain/macs2_peaks/mouse/
# 此参数peak 7万

# 使用新的参数call peak试试
mkdir -p /mnt/xuruizhi/brain/macs2_peaks_new/mouse
cd /mnt/xuruizhi/brain/Tn5_shift/mouse
cat name_new.list | while read id;
do 
  echo $id 
  arr=($id)
  sample=${arr[0]}

  macs2 callpeak  -g mm \
  --nomodel --call-summits --nolambda --keep-dup all -p 0.01 \
  --shift -75 --extsize 150 -n ${sample} -t ./${sample}.Tn5.bed \
  --outdir /mnt/d/ATAC/macs2_peaks2_new/ 
  # --keep-dup all会保留所有开头结尾相同的peak，不加
done
mkdir -p /mnt/d/brain/brain/macs2_peaks_new/mouse
cp /mnt/d/ATAC/macs2_peaks2_new/* /mnt/d/brain/brain/macs2_peaks_new/mouse/
cp /mnt/d/ATAC/macs2_peaks2_new/* /mnt/xuruizhi/brain/macs2_peaks_new/mouse/
# 此参数peak 24万



# 以此peak为准
mkdir -p /mnt/xuruizhi/brain/macs2_peaks_final/mouse
cd /mnt/xuruizhi/brain/Tn5_shift/mouse
cat name_new.list | while read id;
do 
  echo $id 
  arr=($id)
  sample=${arr[0]}

  macs2 callpeak  -g mm \
  --shift -75 --extsize 150 --nomodel \
  --nolambda --keep-dup all \
  -n ${sample} -t ./${sample}.Tn5.bed \
  --outdir /mnt/xuruizhi/brain/macs2_peaks_final_rm/mouse 
done
mkdir -p /mnt/d/brain/brain/macs2_peaks_final/mouse
cp /mnt/xuruizhi/brain/macs2_peaks_final/mouse/*  /mnt/d/brain/brain/macs2_peaks_final/mouse/
```

# 8. Visualization    
1. 目的： 将上文产生的文件放在`IGV`中可视化  


## 8.1 filterbam2Bw    

1. 目的： bw文件是用于方便可视化peak的文件，因为上游处理完的bam文件通常都较大，不方便于快速展示，而将其转变成bw(bigwig)或者wig就会方便的多，而bigWig文件的显示性能又相较wig文件快得多，故bw是更常用的。而相较于bed文件相说，它不只提供了peak的位置，还有peak的高低。 

* bam转bw: 因为此处不看细节位置，不看共同peak，所以使用final.bam文件  

```bash 
mkdir -p  /mnt/xuruizhi/brain/bw/mouse/
cd /mnt/xuruizhi/brain/final/mouse #该目录下需要包含最终过滤后的bam文件和其bai索引
ls *.bam | while read id; 
do 
  bamCoverage -p 6  -b $id \
  -o ../../bw/mouse/${id%%.*}.bw \
  --binSize 20 \
  --smoothLength 60 \
  --normalizeUsing RPKM \
  --centerReads 
  # 1 > ../../bw/mouse/${id%%.*}_bamCoverage.log
done
mkdir -p /mnt/d/brain/brain/bw/mouse
cp /mnt/xuruizhi/brain/bw/mouse/*  /mnt/d/brain/brain/bw/mouse/
# bamCoverage注意大小写
# --binSize Size of the bins, in bases, for the output of the bigwig/bedgraph file. (Default: 50)
# --smoothLength The smooth length defines a window, larger than the binSize, to average the number of reads.
# 可选--blackListFileName BED file  A BED or GTF file containing regions that should be excluded from all analyses.  
# --normalizeUsing {RPKM,CPM,BPM,RPGC,None} Use one of the entered methods to normalize the number of reads per bin. 
# （bw文件夹中last.bam文件使用CPM标准化）--normalizeTo1x: 按照1x测序深度(reads per genome coverage, RPGC)进行标准化
# --centerReads         By adding this option, reads are centered with respect to the fragment length. For paired-end
#                         data, the read is centered at the fragment length defined by the two ends of the fragment. For
#                         single-end data, the given fragment length is used. This option is useful to get a sharper
#                         signal around enriched regions. (default: False)
```

## 8.2 TSS enrichment  执行报错，还未完成

目的：通过观察 peaks 围绕 TSS 的分布情况，判断数据与理论推理是否一致；若一致则证明测序正常。  


来自 NFR（没有核小体的区域） 的片段预计会在基因的转录起始位点 (transcription start site, TSS) 附近富集，而来自核小体结合区域的片段预计会在 TSS 附近被耗尽，在 TSS 附近的侧翼区域会有少量富集 。  [(Fig. 1c)](https://github.com/outcastaaa/ATAC/blob/main/pictures/1c.png)  

① make dir
```bash
mkdir -p /mnt/xuruizhi/brain/TSS/mouse/
```
② 下载TSS注释文件：the BED file which contains the coordinates for all genes [下载地址](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=6884883_WoMR8YyIAAVII92Rr1Am3Kd0jr5H&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr12%3A56703576-56703740&hgta_outputType=primaryTable&hgta_outFileName=)   
[参数选择](https://www.jianshu.com/p/d6cb795af22a)   

genome:mouse --> assemble:mm10 --> gruop:genes and gene predictions --> track:UCSC genes or NCBI RefSeq --> table:如果track选择NCBI RefSeq，这里就选择RefSeq；如果track选择UCSC gene，这里就选knownGene --> output format根据自己的需求选择 --> file type returned这里选gzip compressed，这样就可以下载到压缩包格式的输出文件，选text则下载文本格式 --> output file一定要写上一个文件名字，如果为空则后面无法下载，而只能在浏览器上查看 --> 最后点击get output即可  

将`mm10.reseq.bed`保存在 /mnt/d/ATAC/TSS 文件夹内。 
```bash
cp /mnt/d/ATAC/TSS/mm10.refseq.bed /mnt/xuruizhi/brain/TSS/mouse/
``` 

③ 对比对后的bam文件转化为`bw文件`，保存在  /mnt/xuruizhi/brain/bw/mouse/ 文件夹内  

④ 绘图  
`computeMatrix`根据所提供的refseq.bed文件计算bw文件中在TSS附近左右信号强度，选取的左右可以直接调；若某些转录本附近没有reads，不会计算该位点的信号强度，也可以做自己得到的peaks附近的信号强度。    

用`plotHeatmap`以热图的方式对覆盖进行可视化，用`plotProfile`以折线图的方式展示覆盖情况，该图本质上是一个密度图，用于评估所有转录起始位点的reads密度。  

computeMatrix具有两个模式: `scale-region` 和 `reference-point`。前者用来信号在一个区域内分布，后者查看信号相对于某一个点的分布情况。无论是那个模式，都有两个参数是必须的，-S是提供bigwig文件，-R是提供基因的注释信息。还有更多个性化的可视化选项。  

* 每个样本单独画图  
```bash
pip install "numpy<1.24"

cd /mnt/xuruizhi/brain/bw/mouse/
mkdir -p /mnt/xuruizhi/brain/TSS/mouse
# 循环
ls *.bw | while read id; 
do 
  computeMatrix reference-point --referencePoint TSS -p 6 \
    -b 1000  -a 1000 \
    -R /mnt/xuruizhi/brain/TSS/mouse/mm10.refseq.bed \
    -S $id \
    --skipZeros \
    -o /mnt/xuruizhi/brain/TSS/mouse/$id_matrix.gz \
    --outFileSortedRegions /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_regions.bed
    1 > /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}.log
done
# --referencePoint Possible choices: TSS, TES, center
# -b, --upstream Distance upstream of the reference-point selected. (Default: 500)
# -a, --downstream Distance downstream of the reference-point selected. (Default: 1500)
# --missingDataAsZero  If set, missing data (NAs) will be treated as zeros. The default is to ignore such cases, which will be depicted as black areas in a heatmap.
# --skipZeros Whether regions with only scores of zero should be included or not. Default is to include them.  
# --binSize Length, in bases, of the non-overlapping bins for averaging the score over the regions length. (Default: 10)  
# --blackListFileName, -bl A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry.
# Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered.
# --binSize BINSIZE 几个bp分数取平均，默认:10bp  

# profile plot
ls *.log | while read id; 
do 
  plotProfile -m /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_matrix.gz \
    -out /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_profile.png \
    --perGroup \
    --colors green \
    --plotTitle "" \
    --refPointLabel "TSS" \
    -T "${id%%.*} read density" \
    -z ""
done
#--perGroup            The default is to plot all groups of regions by sample. Using this option instead plots all
                        # samples by group of regions. Note that this is only useful if you have multiple groups of
                        # regions. by sample rather than group. (default: False)


# heatmap and profile plot
ls *.log | while read id; 
do 
  plotHeatmap -m /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_matrix.gz \
    -out /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_heatmap.png \
    --colorMap RdBu \
    --zMin -12 --zMax 12
done


#单独heatmap
ls *.log | while read id; 
do
plotHeatmap -m /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_matrix.gz \
-out /mnt/xuruizhi/brain/TSS/mouse/${id%%.*}_heatmap2.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar' \
--zMin -8 --zMax 8  
done
```



* 画 `gene body` 区，使用 `scale-regions`  
```bash
cd /mnt/xuruizhi/brain/bw/mouse
mkdir -p /mnt/xuruizhi/brain/genebody/mouse/
# create a matrix 
computeMatrix scale-regions -p 6 \
    -b 10000  -a 10000 \
    -R /mnt/xuruizhi/brain/TSS/mouse/mm10.refseq.bed \
    -S SRR11539111.bw \
    --skipZeros \
    -o /mnt/d/ATAC/genebody/SRR11539111_matrix.gz 
  
cd /mnt/d/ATAC/genebody
plotHeatmap -m /mnt/d/ATAC/genebody/SRR11539111_matrix.gz \
    -out /mnt/d/ATAC/genebody/SRR11539111_heatmap.png 

plotProfile -m /mnt/d/ATAC/genebody/SRR11539111_matrix.gz \
    -out /mnt/d/ATAC/genebody/SRR11539111_profile.png 
    #不太好看，还需要调整参数
```
# 

# 9. 使用diffbind做主成分分析


① read in a set of peaksets and associated metadata  


* 原理：DiffBind 所需的输入是数据集中的所有样本以及每个样本的所有峰（不仅仅是高置信度峰），合并函数会查找 `overlap peak` 的基因组区间，如果某区间出现在两个及以上的样本中，定义为`consensus peakset`；具有rep需要单独使用，不可合并（因此在寻找差异peak时，不可使用IDR找到的consensus peak）。  

* [具体参数](https://rdrr.io/bioc/DiffBind/man/dba.html)  

* 输入文件：[参考官网man](https://rdrr.io/bioc/DiffBind/man/dba.html)    
文件格式：CSV表（，分隔）；表格.xls/xlsx  
sample sheet是一个列表，需要包括以下几列:"SamplelD"，"Tissue"，"Factor"，"Condition"， "Treatment"，"Replicate"， "bamReads"，"ControllD"，"bamControl"，"Peaks"和"PeakCaller"   

* 内容格式：  

| header  | detial  |
|:---:|:---:|
|  SampleID |  样本ID，给你输入的数据起个名  |
| Tissue，Factor，Condition，Treatment  |  都为数据的备注，包括组织来源/细胞系，状态，处理等，可不填，但是会影响后面分析的聚类。factor不是很重要；Treatment就是分组，对照或者不同处理，也可以是对照和过表达/KO等 |
| Replicate  |  第几次重复 |
| bamReads  |  ChIP-seq得到的bam文件，bam文件的绝对路径|
| ControlID  | Call peak时使用的input数据的ID，ATAC不需要  |
|  bamControl |  input对应的bam文件，ATAC不需要 |
|  Peaks | 峰文件，这里有多种数据格式可作为输入：1. macs2 输出的.narrowPeak等峰文件 2. 包括所有call peak 得到的peak位置信息的.bed 文件，不是Macs2直接得到的bed文件 3. 以上两种格式得到的.gz文件 |
| PeakCaller  |  用何种方式做的peak calling，默认峰值格式：narrowPeaks 文件 |   



* 输入：注意！！！一定把对照组放前面，把实验组放在后面  
```bash
mkdir -p /mnt/d/brain/brain/final/mouse/
cp /mnt/xuruizhi/brain/final/mouse/* /mnt/d/brain/brain/final/mouse/
mkdir -p /mnt/d/brain/brain/macs2_peaks/mouse/
cp /mnt/xuruizhi/brain/macs2_peaks/mouse/* /mnt/d/brain/brain/macs2_peaks/mouse/

```

| SampleID | Tissue         | Factor              | Condition | Treatment | Replicate | bamReads                                | ControlID | bamControl | Peaks                                               | PeakCaller |
|:--------:|:--------------:|:-------------------:|:---------:|:---------:|:---------:|:---------------------------------------:|:---------:|:----------:|:---------------------------------------------------:|:----------:|
| PFC1    | PFC | accessible\_regions | PFC      | PFC      | 1         | D:/brain/brain/final/mouse/SRR14362276\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR14362276\_peaks\.narrowPeak | narrowPeak |
| PFC2    | PFC | accessible\_regions | PFC      | PFC      | 2         | D:/brain/brain/final/mouse/SRR14362281\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR14362281\_peaks\.narrowPeak | narrowPeak |
| PFC3      | PFC    | accessible\_regions | PFC        | PFC        | 3         | D:/brain/brain/final/mouse/SRR14362282\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR14362282\_peaks\.narrowPeak | narrowPeak |
| cortex1      | cortex    | accessible\_regions | cortex        | cortex        | 1         |  D:/brain/brain/final/mouse/SRR13049359\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR13049359\_peaks\.narrowPeak | narrowPeak |  
| cortex2      | cortex    | accessible\_regions | cortex        | cortex        | 2         |  D:/brain/brain/final/mouse/SRR13049361\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR13049361\_peaks\.narrowPeak | narrowPeak |
| HIPP1      | HIPP    | accessible\_regions | HIPP        | HIPP        | 1         |  D:/brain/brain/final/mouse/SRR11179779\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR11179779\_peaks\.narrowPeak | narrowPeak |
| HIPP2      | HIPP    | accessible\_regions | HIPP        | HIPP        | 2         |  D:/brain/brain/final/mouse/SRR11179780\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR11179780\_peaks\.narrowPeak | narrowPeak |
| HIPP3      | HIPP    | accessible\_regions | HIPP        | HIPP        | 3         |  D:/brain/brain/final/mouse/SRR11179781\.final\.bam |           |            | D:/brain/brain/macs2_peaks/mouse/SRR11179781\_peaks\.narrowPeak | narrowPeak |


将上面表格写入文件`/mnt/d/brain/brain/R_analyse/mouse/sample_sheet.csv`，学会使用[格式转换器](https://tableconvert.com/zh-cn/csv-to-excel)，注意csv文件最后一行加一行空格，否则报错。注意删掉斜杠，写入csv不能有斜杠。




* 代码：

```r
# 在 R.studio 中进行操作
# 下载R包
BiocManager::install("DiffBind", force = TRUE)
library(DiffBind)
# DiffBind 3.8.4
getwd()
# [1] "D:/brain/brain/R_analyse"


# 导入数据
> samples <- read.csv("./sample_sheet.csv")
> names(samples)
#  [1] "SampleID"   "Tissue"     "Factor"     "Condition"  "Treatment" 
#  [6] "Replicate"  "bamReads"   "ControlID"  "bamControl" "Peaks"     
# [11] "PeakCaller"

#找到样本间共有peaks，比较相似性
> dbObj <- dba(sampleSheet = samples)  
> dbObj
# 8 Samples, 180339 sites in matrix (265042 total):
#        ID Tissue             Factor Condition Treatment Replicate Intervals
# 1    PFC1    PFC accessible_regions       PFC       PFC         1    143526
# 2    PFC2    PFC accessible_regions       PFC       PFC         2    146365
# 3    PFC3    PFC accessible_regions       PFC       PFC         3    119528
# 4 cortex1 cortex accessible_regions    cortex    cortex         1    144230
# 5 cortex2 cortex accessible_regions    cortex    cortex         2    223038
# 6   HIPP1   HIPP accessible_regions      HIPP      HIPP         1     75243
# 7   HIPP2   HIPP accessible_regions      HIPP      HIPP         2     84122
# 8   HIPP3   HIPP accessible_regions      HIPP      HIPP         3     66794
``` 
* 结果解读：    

This shows how many peaks are in each peakset, as well as (in the first line) the total number of unique peaks after merging overlapping ones (`265042`), and the dimensions of the default binding matrix of `8` samples by the `180339` sites that overlap in at least two of the samples.

* heatmap: 生成一个相关热图，利用矩阵的每一行的互相关联cross-correlations来给出样本的初始聚类  
```r
> plot(dbObj)
```


② Counting reads and creating a binding affinity matrix    

* 原理：   
The next step is to calculate a binding matrix with scores based on read counts for every sample (affinity scores), rather than confidence scores for only those peaks called in a specific sample (occupancy scores). 一旦一个 `consensus peak` 被推导出来，DiffBind可以使用提供的测序read文件来计算每个样本的每个区间有多少reads重叠。默认情况下，为了提供更多标准化的峰值区间，consensus peak中的峰会根据其峰值(最大读重叠点)重新调整中心点和trimmed。计数的最终结果是一个结合亲和矩阵，其中包含每个样本在每个共识结合位点的read count.  


* [具体参数](https://rdrr.io/bioc/DiffBind/man/dba.count.html)  

* 代码：  

```r
> db_count <- dba.count(dbObj)  #this step will take you a couple of minutes, be patient.
# 8 Samples, 180205 sites in matrix:
#        ID Tissue             Factor Condition Treatment Replicate    Reads FRiP
# 1    PFC1    PFC accessible_regions       PFC       PFC         1 39421908 0.25
# 2    PFC2    PFC accessible_regions       PFC       PFC         2 31444995 0.28
# 3    PFC3    PFC accessible_regions       PFC       PFC         3 30715211 0.22
# 4 cortex1 cortex accessible_regions    cortex    cortex         1  5826216 0.51
# 5 cortex2 cortex accessible_regions    cortex    cortex         2 31773714 0.43
# 6   HIPP1   HIPP accessible_regions      HIPP      HIPP         1 11534423 0.20
# 7   HIPP2   HIPP accessible_regions      HIPP      HIPP         2 27587078 0.16
# 8   HIPP3   HIPP accessible_regions      HIPP      HIPP         3 17589852 0.15
```

* 添加参数：  
`bUseSummarizeOverlaps`，这个参数会使得运行比较缓慢，但是是一个更标准的计算功能。如果你把它设置为TRUE，所有的read文件必须是bam，并且必须有其自己的索引文件 (.bam.bai) 。另外fragmentSize参数必须是缺省值。
```r
> db_count2 <- dba.count(dbObj,bUseSummarizeOverlaps=TRUE)
> db_count2
# 8 Samples, 180205 sites in matrix:
#        ID Tissue             Factor Condition Treatment Replicate    Reads FRiP
# 1    PFC1    PFC accessible_regions       PFC       PFC         1 39421908 0.25
# 2    PFC2    PFC accessible_regions       PFC       PFC         2 31444995 0.28
# 3    PFC3    PFC accessible_regions       PFC       PFC         3 30715211 0.22
# 4 cortex1 cortex accessible_regions    cortex    cortex         1  5826216 0.51
# 5 cortex2 cortex accessible_regions    cortex    cortex         2 31773714 0.43
# 6   HIPP1   HIPP accessible_regions      HIPP      HIPP         1 11534423 0.20
# 7   HIPP2   HIPP accessible_regions      HIPP      HIPP         2 27587078 0.16
# 8   HIPP3   HIPP accessible_regions      HIPP      HIPP         3 17589852 0.15 

> dba.plotPCA(db_count2, attributes=DBA_TREATMENT, label=DBA_ID)
> plot(db_count2)
```



* 通过 `dba.show` 命令整合：  
```r
> info <- dba.show(db_count2)
> libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
> rownames(libsizes) <- info$ID
> libsizes
#         LibReads FRiP PeakReads
# PFC1    39421908 0.25   9855477
# PFC2    31444995 0.28   8804598
# PFC3    30715211 0.22   6757346
# cortex1  5826216 0.51   2971370
# cortex2 31773714 0.43  13662697
# HIPP1   11534423 0.20   2306885
# HIPP2   27587078 0.16   4413932
# HIPP3   17589852 0.15   2638478
```

# 10. peak annotation

1. 目的： 

 获得 Peak 后，Peak 的注释可将染色质的可及性与基因调控联系起来。通常，Peak 由最接近的基因或调控元件进行注释。通常，来自 ATAC-seq 的 Peak 将代表不同的顺式调节元件的混合物，包括增强子和启动子 。在获得基因组特征列表之后，还可以使用 GO,KEGG和 Reactome等数据库进行功能富集分析。通常，Peak 注释会产生生物学和功能上有意义的结果，以供进一步研究。

2. 使用软件:   

 `ChIPseeker`和 `ChIPpeakAnno`被广泛用于为最接近或重叠的基因、外显子、内含子、启动子、5'UTR、3'UTR 和其他基因组特征分配 Peak。ChIPseeker 和 ChIPpeakAnno 还具有丰富的可视化功能，可用于解释注释结果，例如带有注释的基因组特征的饼图。  

3. 输入文件：  

一个就是要注释的peaks的文件，需满足BED格式。另一个就是注释参考文件，即需要一个包含注释信息的TxDb对象。  


4. 步骤： 

① 导入文件做基础分析  

```r
# Load libraries
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicFeatures")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", force = TRUE)
BiocManager::install("org.Mm.eg.db", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)

library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)


# 导入bed文件
> getwd()
# [1] "D:/brain/brain/R_analyse"
# 单个导入文件
> peak <- readPeakFile("D:/brain/brain/macs2_peaks/mouse/SRR3595213_summits.bed",sep ="") 

## 循环导入文件
# 创建一个存储文件名的向量
file_names <- c("SRR3595213_summits.bed", "SRR11179779_summits.bed", "SRR11179780_summits.bed", "SRR11179781_summits.bed", "SRR13049359_summits.bed","SRR13049361_summits.bed","SRR14362276_summits.bed","SRR14362281_summits.bed","SRR14362282_summits.bed")
# 创建一个空的列表，用于存储读取的数据
peak_list <- list()
# 使用循环遍历文件名向量
for (file_name in file_names) {
  # 构建文件路径
  file_path <- paste0("D:/brain/brain/macs2_peaks/mouse/", file_name)  # 将路径替换为您的实际路径
  # 使用适当的读取函数（readPeakFile）读取数据
  peak_data <- readPeakFile(file_path, sep = "")
  # 将数据添加到列表中
  peak_list[[file_name]] <- peak_data
}

# 现在，您可以通过列表中的名称访问读取的数据
# 例如，访问第一个文件的数据：
peak_list[["SRR3595213_summits.bed"]]



# peak在染色体上的分布
> covplot(peak)

# peak 在TSS位点附件的分布
> txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
  tagMatrix <- getTagMatrix(peak, windows=promoter)
> tagHeatmap(tagMatrix, xlim=c(-1000, 1000), color="red")
> plotAvgProf(
  tagMatrix,
  xlim=c(-1000, 1000),
  xlab="Genomic Region (5'->3')",
  ylab = "Peak Frequency")
```

* 结果注释:    

默认情况下，annotatePeak 函数使用 TSS 方法，并提供参数来指定最大距离截止。还可以报告此距离内的所有基因，无论是否与TSS重叠。对于注释基因组区域，当基因组区域为外显子或内含子时，annotatePeak 函数报告详细信息。例如，“外显子（uc002sbe.3/9736, exon 69 of 80）”意味着峰与转录本uc69sbe.80拥有的80个外显子的第69个外显子重叠，相应的Entrez基因ID为9736。   
![peak_chrs](./pictures/peak_chrs.png)  
![tss_heatmap](./pictures/tss_heatmap.png)   
![tss_line](./pictures/tss_line.png)    


② peak关联基因注释  
给出了关联的基因以及对应的基因组区域的类别，根据这个结果，可以提取关联基因进行下游的功能富集分析，比如提取geneid这一列，用clusterProfiler进行GO/KEGG等功能富集分析。  

```r
> peakAnnolist <- annotatePeak(
    peak,
    tssRegion = c(-1000, 1000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db")
# Annotated peaks generated by ChIPseeker
# 8265/8265  peaks were annotated
# Genomic Annotation Summary:
#              Feature   Frequency
# 9           Promoter 82.92800968
# 4             5' UTR  0.24198427
# 3             3' UTR  0.54446461
# 1           1st Exon  1.12522686
# 7         Other Exon  0.96793708
# 2         1st Intron  3.08529946
# 8       Other Intron  4.69449486
# 6 Downstream (<=300)  0.07259528
# 5  Distal Intergenic  6.33998790


> write.table(
    as.data.frame(peakAnnolist),
    "diff_DESeq2peak.annotation.tsv",
    sep="\t",
    row.names = F,
    quote = F)   

#可视化
> plotAnnoPie(peakAnnolist)
> plotAnnoBar(peakAnnolist) 
> plotDistToTSS(peakAnnolist,title="Distribution of accessible regions relative to TSS")
```  

* 结果注释:   

![peakpie](./pictures/peakpie.png)   

可以看到染色质可及区域绝大部分分布在TSS区域。 结果： annotatePeak传入annoDb参数,即可进行基因ID转换，将Entrez ID转化为ENSEMBL，SYMBOL，GENENAME.  


```bash
seqnames  start    end     width   strand   V4      V5            annotation          geneChr   geneStart    geneEnd geneLength    geneStrand   geneId  transcriptId         distanceToTSS   
ENSEMBL                SYMBOL          GENENAME
chr9    58807687 58808086   400     *       * -1.93607772246243   Distal Intergenic       9      58823412    58863175   39764           1       330953  ENSMUST00000034889.9    -15326 
 ENSMUSG00000032338      Hcn4    hyperpolarization-activated, cyclic nucleotide-gated K+ 4
```



















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

6. 预印本：灵长类动物精细运动系统的调节演化

数据分析:我们使用编码ATAC-seq管道(Landt等人，2012年)处理了ATAC-seg实验的原始FASTQ文件，访问网址为https://github.com/ENCODE-DCC/atac-seg-管道中。为了补充我们的猕猴和大鼠的数据，我们从(Srinivasan等，2020) 获得了小鼠大脑皮层和纹状体的ATAC-seq数据。我们还处理了来自人类死后大脑(Fullardetal..2018)的公开可用的NeuN排序的ATAC-seg数据，这些数据来自大致对应于从猕猴收集的区域，即:初级运动皮层(PMC)、腹外侧前额叶皮层 (VLPFC)背外侧前额叶皮层(dIPFC)、题上回(STC)、伏隔核(NAc)和壳核(PUT)。我们通过基因表达综合(GEO)登录IDGSE96949从序列读取存档(SRA)下载了这些数据。

我们使用猕猴的 rheMac8 组件、人类的 hg38 组件、大鼠的 m6 组件和小鼠的 mm10 组件运行 ENCODE 管道。 我们使用除“atac.multimapping”之外的默认参数运行管道：0、“atac.cap num peak”300.000。 “atac.smooth win”：150，“atac.enable idr”：true，以及“atacidr thresh”：0.1。 我们在处理数据时结合了技术复制。 我们为每个复制生成过滤的 bam 文件、峰值文件和信号轨道，并为每个组织生成复制池，每个物种。我们删除了 ENCODE 质量控制指标指示的低周期性样本，并重新处理了剩余的复制。 由于我们的重复通常在测序深度上存在显着差异，我们将可重现的峰定义为具有不可重现的发现率 irreproducible discovery rate (IDR) < 0.1 across pooled pseudo-replicates,（IDR，（Li 等人，2011））< 0.1 跨合并伪重复的峰，并将这些峰用于所有下游分析。 对于只有一个高质量生物复制的组织，我们使用了根据 IDR < 0.1 在自我伪复制中可重现的峰。

除了确定单个组织的峰集，对于每个物种，我们确定了一组峰作为全基因组背景集代表从所有加工组织确定的可重复的开放染色质峰的联盟。该背景集是使用bedtools (Ouinlan和Hall，2010年)与-wa和-u选项相交以结合每个物种的所有可再现峰集获得的。采取了许多步骤来准备OCR峰集用于下游分析。使用bedtools合并在50bp内的峰。我们使用bedtools subtract with option-A从任何注释的编码或非编码外显子中去除2 kb以内的峰，使我们能够从背景中排除启动子、编码序列和非编码RNA。与Y染色体对齐的峰被删除以控制性别偏见的影响。为了确定称猴的外显子排除区的完整集合，我们使用了完整的rheMac8 RefSeg注释(0’Leary等人，2016)，并补充了从他们的轨道浏览器获得的USC的“xenoRefSeg”注释，这代表了来自几十个其他物种通过liftOver与猕猴对齐(Karolchik等人，2004年)基因的RefSeg注释。用于人、大鼠和小鼠;我们分别使用了hq38、RN6和MM10RefSeg注释集(0“Leary等人，2016)。

为了识别不同组织对比中不同活性的OCR峰，我们首先使用特征数量化了每个组织与该物种的一致峰值相一致的读取次数(廖等人，2014年)。然后，我们使用DESeq2R软件包中的负二项分布模型，对比了各组织之间的重新计数(Love等人，2014年)。我们认为峰差在组织对比中呈对数倍差>1.5，调整后的p值<0.05。

为了确定跨物种的直系同源OCRS，我们将每个物种的每个组织的开放染色质数据与本研究中考虑的所有其他物种进行了比对。使用hallift over (Hickev等人，2013)和默认参数在物种之间映射OCRS，使用ZoonomiaCactus多个全基因组序列比对进行基于图形的基因组坐标映射(Armstrong等人，2019)。halLiftover的原始输出随后被过滤，并使用用于调控元件进化的halLiftover后处理(HALPER)工具组装成连续的OCR(Zhang等人，2020)，参数为-max frac12-minlen50和-protect_dist5。

使用 Genomic Regions Enrichment ofAnnotations Tool (GREAT) 4.0.4 版（McLean 等人，2010）进行基因本体论分析。 差异 OCR 峰集的基因组坐标用作前景区域。 对于背景区域，我们使用了从每个物种的所有处理过的组织中识别出的所有可重现的开放染色质峰的联合。 显着过度代表的本体类别按超几何错误发现率 q 值排序，并且仅考虑由至少 5 个基因组成的 GO 术语。
为了鉴定相对于打乱序列的感兴趣的差异OCR峰集中富集的转录因子结合基序，我们在MEME套件中使用了AME，对总优势分数进行费舍尔的精确测试（该序列的位置权重矩阵（PWM）基序得分总分），所有参数为默认。对于我们的PWM集，我们使用了JASPAR2018 CORE非冗余脊椎动物集motif。
