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
--very-sensitive -X 2000 -1  -U {}_trimmed.fq.gz \
2> ~/xuruizhi/brain/brain/alignment_new/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment_new/mouse/{}.sam

cat single.list  | while read id; do sed "s/{}/${id}/g" 222.sh > ${id}.sh; done

cat single.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/alignment_new/mouse/ "bash ${id}.sh"
done
# Job <8601568> is submitted to queue <mpi>.
# Job <8601569> is submitted to queue <mpi>.
# Job <8601570> is submitted to queue <mpi>.
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
mkdir -p ~/xuruizhi/brain/brain/sort_bam_new/mouse
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
samtools index -@ 48 ~/xuruizhi/brain/brain/sort_bam/mouse/{}.sort.bam
samtools flagstat  -@ 48 ~/xuruizhi/brain/brain/sort_bam/mouse/{}.sort.bam > ~/xuruizhi/brain/brain/sort_bam/mouse/{}.raw.stat

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



cd ~/xuruizhi/brain/brain/sort_bam_new/mouse
cat rmdup.sh
#!/usr/bin bash
parallel -k -j 24 '
picard MarkDuplicates -I {}.sort.bam -O ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam \
-REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT \
-METRICS_FILE ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.log'
samtools index -@ 48 ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam
samtools flagstat -@ 48 ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.bam > ~/xuruizhi/brain/brain/rmdup_new/mouse/{}.rmdup.stat

cat name.list  | while read id; do sed "s/{}/${id}/g" rmdup.sh > ${id}.sh; done

cat name.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/rmdup_new/mouse "bash ${id}.sh"
done
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
cp ~/xuruizhi/brain/brain/sort_bam_new/mouse/name.list ~/xuruizhi/brain/brain/filter_new/mouse/



cd ~/xuruizhi/brain/brain/rmdup_new/mouse
cat filter.sh
#!/usr/bin bash
parallel -k -j 8 ' samtools view -h -f 2 -F 1804 -q 30 {}.rmdup.bam | grep -v  chrM | samtools sort -@ 6 -O bam  -o ~/xuruizhi/brain/brain/filter_new/mouse/{}.filter.bam'
samtools index -@ 48 ~/xuruizhi/brain/brain/filter_new/mouse/{}.rmdup.bam
samtools flagstat -@ 48 ~/xuruizhi/brain/brain/filter_new/mouse/{}.rmdup.bam > ~/xuruizhi/brain/brain/filter_new/mouse/{}.rmdup.stat

cat name.list  | while read id; do sed "s/{}/${id}/g" filter.sh > ${id}.sh; done

cat name.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/filter_new/mouse "bash ${id}.sh"
done
```

## 5.4 Blacklist filtering

1. 目的：去除ENCODE blacklisted 区域，通过blacklist的过滤，可以进一步降低peak calling的假阳性。    

本流程采用的方法是：在peak calling之前去除，比对后的reads 去除PCR重复等后单独去除 blacklist region，再 call peak.  
   

```bash
# 下载软件
cd ~/xuruizhi/biosoft/biosoft/bedtools
tar -zxvf bedtools-2.30.0.tar.gz
cd bedtools2
make



# 下载对应物种的 blacklist.bed文件
mkdir -p /mnt/d/ATAC/blklist
cd /mnt/d/ATAC/blklist
wget https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
gzip -dc mm10.blacklist.bed.gz > mm10.blacklist.bed
rm *.gz
wc -l  mm10.blacklist.bed #164

cd /mnt/d/ATAC/filter
cat config.raw | while read id;
do 
  echo $id 
  arr=($id)
  sample=${arr[0]}

  echo ${sample}.filter.bam

  # 取交集看bam文件和blacklist有多少重合部分
  bedtools intersect -wa -a ${sample}.filter.bam  -b ../blklist/mm10.blacklist.bed | wc -l  
  # 16559
  # 15119
  # 15304
  # 20212

  # 凡是bam中含有blacklist都删除
  bedtools intersect -v -a ${sample}.filter.bam -b ../blklist/mm10.blacklist.bed > ../blklist/${sample}.final.bam
  samtools index  -@ 7 ../blklist/${sample}.final.bam
  samtools flagstat  -@ 7 ../blklist/${sample}.final.bam > ../blklist/${sample}.final.stat
done


cat config.raw | while read id;
do 
  echo $id 
  arr=($id)
  sample=${arr[0]}

  samtools index  -@ 7 ../blklist/${sample}.final.bam 
	samtools flagstat  -@ 7 ../blklist/${sample}.final.bam > ../blklist/${sample}.final.stat
done
```
6. 结果解读：  
```bash
# 原比对文件数据，以SRR11539111为例
98013300 + 0 in total (QC-passed reads + QC-failed reads)
98013300 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
96440057 + 0 mapped (98.39% : N/A)
96440057 + 0 primary mapped (98.39% : N/A)
98013300 + 0 paired in sequencing
49006650 + 0 read1
49006650 + 0 read2
94727152 + 0 properly paired (96.65% : N/A)
95584080 + 0 with itself and mate mapped
855977 + 0 singletons (0.87% : N/A)
160994 + 0 with mate mapped to a different chr
89323 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除PCR重复+低质量+chrM后数据
48111744 + 0 in total (QC-passed reads + QC-failed reads)
48111744 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
48111744 + 0 mapped (100.00% : N/A)
48111744 + 0 primary mapped (100.00% : N/A)
48111744 + 0 paired in sequencing
24055872 + 0 read1
24055872 + 0 read2
48111744 + 0 properly paired (100.00% : N/A)
48111744 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

# 删除blacklist后数据
47997002 + 0 in total (QC-passed reads + QC-failed reads)
47997002 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
47997002 + 0 mapped (100.00% : N/A)
47997002 + 0 primary mapped (100.00% : N/A)
47997002 + 0 paired in sequencing
23998484 + 0 read1
23998518 + 0 read2
47997002 + 0 properly paired (100.00% : N/A)
47997002 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
到这一步，比对文件已经过滤完成。     


## 5.5 bamtobed
1. 目的：后续需要用到 `bed bedpe` 文件，把处理好的bam比对文件转化为bed格式
2. 使用软件：`bedtools`,[参考文章](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html)  
3. 代码：
```bash
# bam to bed
mkdir -p /mnt/d/ATAC/bed
cd /mnt/d/ATAC/blklist

parallel -j 6 "
   bedtools bamtobed -i ./{1} >../bed/{1}.bed
" ::: $( ls *.final.bam)

# bam to bedpe 
mkdir -p /mnt/d/ATAC/bedpe
cd /mnt/d/ATAC/blklist
# the BAM file should be sorted by read name beforehand
parallel -j 6 "
  samtools sort -n -o ../bedpe/{1}.named {1}
" ::: $( ls *.final.bam)

cd /mnt/d/ATAC/bedpe
cat config.raw | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}
  samtools flagstat  -@ 7 ${sample}.final.bam.named > ${sample}.final.bam.named.stat
done
  


cd /mnt/d/ATAC/bedpe
# bedtools should extract the paired-end alignments as bedpe format, then awk should shift the fragments as needed
parallel -j 6 "
  bedtools bamtobed -i {1} -bedpe > {1}.bedpe
" ::: $( ls *.final.bam.named)
```
注：bedpe转化一定要按照name排序，把双端reads放一起；因为去除blacklist后有些reads被去除无法组成一个pair被skip
* 结果：
```bash
# bed
$ cat SRR11539111.final.bam.bed | head -n 5
chr1    3000773 3000873 SRR11539111.41226980/2  32      +
chr1    3000784 3000884 SRR11539111.41226980/1  32      -
chr1    3000793 3000893 SRR11539111.46953273/1  34      +
chr1    3000873 3000969 SRR11539111.16779100/1  36      +

# bedpe
$ cat SRR11539111.final.bam.named.bedpe | head -n 5
chr16   79178081        79178149        chr16   79178181        79178281        SRR11539111.1   42      +       -
chr2    64769626        64769726        chr2    64769944        64770041        SRR11539111.3   40      +       -
chr13   31981784        31981881        chr13   31981802        31981902        SRR11539111.6   42      +       -
chr7    45794613        45794710        chr7    45794641        45794740        SRR11539111.12  42      +       -
chr14   122435898       122435949       chr14   122435898       122435949       SRR11539111.15  42      +       -

```
* bedpe文件格式  [bed文件格式](https://www.cnblogs.com/djx571/p/9499795.html#:~:text=BED%20%E6%96%87%E4%BB%B6%28Browser%20Extensible%20Data%29%E6%A0%BC%E5%BC%8F%E6%98%AFucsc,%E7%9A%84genome%20browser%E7%9A%84%E4%B8%80%E4%B8%AA%E6%A0%BC%E5%BC%8F%20%2C%E6%8F%90%E4%BE%9B%E4%BA%86%E4%B8%80%E7%A7%8D%E7%81%B5%E6%B4%BB%E7%9A%84%E6%96%B9%E5%BC%8F%E6%9D%A5%E5%AE%9A%E4%B9%89%E7%9A%84%E6%95%B0%E6%8D%AE%E8%A1%8C%EF%BC%8C%E4%BB%A5%E7%94%A8%E6%9D%A5%E6%8F%8F%E8%BF%B0%E6%B3%A8%E9%87%8A%E4%BF%A1%E6%81%AF%E3%80%82%20BED%E8%A1%8C%E6%9C%893%E4%B8%AA%E5%BF%85%E9%A1%BB%E7%9A%84%E5%88%97%E5%92%8C9%E4%B8%AA%E9%A2%9D%E5%A4%96%E5%8F%AF%E9%80%89%E7%9A%84%E5%88%97%E3%80%82)  

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
thickStart- 绘制特征的起始位置（例如，基因显示中的起始密码子）。当没有厚部分时，thickStart和thickEnd通常设置为chromStart位置。
thickEnd - 绘制特征的结束位置（例如基因显示中的终止密码子）。
itemRgb- R，G，B形式的RGB值（例如255,0,0）。如果轨道行 itemRgb属性设置为“On”，则此RBG值将确定此BED行中包含的数据的显示颜色。\
注意：建议使用此属性的简单颜色方案（八种颜色或更少颜色），以避免压倒Genome浏览器和Internet浏览器的颜色资源。
blockCount- BED行中的块（外显子）数。
blockSizes- 块大小的逗号分隔列表。此列表中的项目数应与blockCount相对应。
blockStarts - 以逗号分隔的块开始列表。应该相对于chromStart计算所有 blockStart位置。此列表中的项目数应与blockCount相对应。

链接：https://www.jianshu.com/p/9208c3b89e44
```
* [bed bedpe格式的区别](https://www.jianshu.com/p/c73c1dc81c61)  
BEDPE 格式类似于 BED 格式，可用于描述成对的基因组区域。
由于bed文件原则上不能表示跨染色体的信息，因此，对于结构变异，一般采用的一种基于bed文件的变种文件bedpe格式进行存储。其格式与bed最大的区别在于，对于必须列即chrom、chromStart、chromEnd三列分别记录两次。  

# 6. shift reads
1. 目的：  

由于Tn5酶是以二聚体的形式结合到染色体上的，其跨度大致是9bp，在第一篇ATAC-seq出来的时候，作者就考虑到了这个问题，在分析的时候，需要回补这个9个bp的碱基差。具体做法就是将正链正向移动4bp，将负链负向移动5个bp。一般用alignmentSieve 一步到位。注意，不做reads shift 对单碱基分辨高的分析会有影响，例如TF motif footprinting，但也不是所有TF footprinting分析软件需要shifted reads，很多可以自己转换，e.g. NucleoATAC。   

方法：
分别对正链和负链的 reads 进行 + 4bp 和 -5bp 的移位（这个长度近似于一个完整的DNA螺旋[why参考文章](https://www.jianshu.com/p/13779b89e76b)），以解释 Tn5 转座酶修复损伤 DNA 所产生的 9bp 的重复，并实现 TF footprint 和 motif 相关分析的碱基对分辨率。  


2. 使用软件：该步有很多种[方法](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/)，本流程采用 `bedtools` and `awk`.

3. 代码：
```bash
mkdir -p /mnt/d/ATAC/Tn5_shift
cp /mnt/d/ATAC/rmdup/config.raw /mnt/d/ATAC/bedpe/config.raw

# bed转化
cd /mnt/d/ATAC/bed/
cat config.raw | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  cat ${sample}.final.bam.bed | awk -v \
  OFS="\t" '{if($6=="+"){print $1,$2+4,$3+4} \
   else if($6=="-"){print $1,$2-5,$3-5}}' \
    > ../Tn5_shift/${sample}.Tn5.bed
done


# bedpe转化
cd /mnt/d/ATAC/bedpe
cat config.raw | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  cat ${sample}.final.bam.named.bedpe | awk -v \
  OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4} \
   else if($9=="-"){print $1,$2-5,$6-5}}' \
    > ../Tn5_shift/${sample}.Tn5.bedpe
done
```
4. 结果解读：  


！注意，后续callpeak不可直接使用bedtools转化的bedpe文件，只能包含三行信息：chr,chrom_start,chrom_end
```bash
cd /mnt/d/ATAC/Tn5_shift
$ cat SRR11539111.Tn5.bed | head -n 5
chr1    3000777 3000877
chr1    3000779 3000879
chr1    3000797 3000897
chr1    3000877 3000973
chr1    3000922 3001022
$ wc -l SRR11539111.Tn5.bed
# 47997002

$ cat SRR11539111.Tn5.bedpe | head -n 5
chr16   79178085        79178285
chr2    64769630        64770045
chr13   31981788        31981906
chr7    45794617        45794744
chr14   122435902       122435953
$ wc -l SRR11539111.Tn5.bedpe
# 23998114
# bedpe文件行数应该是对应bed文件的一半，但是384对被blacklist去除了
```




# 7. Call peaks 
1. 目的： 下一步需要在统计学上判断真实的peak，因为Tn5在染色体上结合是个概率事件，如何判断这个位置的reads足够为一个peak，这就需要用到统计检测。ATAC-seq 数据分析的第二个主要步骤是识别开放区域（也称为 Peak），后续高级分析以此为基础。  

2. 使用软件：目前，`MACS2` 是 ENCODE ATAC-seq 流程的默认 Peak caller 程序。  

3. !!!重要：关于是否使用[-f BEDPE的讨论](https://github.com/macs3-project/MACS/issues/331)，可根据需要选择合适的callpeak参数。  


4. 其他： 


* ATAC-seq关心的是在哪里切断，断点才是peak的中心，所以使用shift模型，--shift -75或-100.   

* 这里选用固定宽度（fixed-width）的peaks,优点有：   
1）对大量的peaks进行counts和motif分析时可以减小误差；  
2）对于大量数据集的可以合并峰得到一致性的peaks;   

* 一个样本的overlaps他们是通过迭代移除的方法，首先保留最显著的peak，然后任何与最显著peak有直接overlap的peaks都被移除；接着对另一个最显著性的peak进行相同的操作，最终保留所有更显著的peaks，移除与其有直接overlaps的peaks  
* 注：后续分析过程需要用到IDR提取consensus peak，建议MACS2 callpeaks的步骤参数设置不要过于严格，以便鉴定出更多的peaks。

4. 代码：
```bash
mkdir -p /mnt/d/ATAC/macs2_peaks/
cd /mnt/d/ATAC/Tn5_shift/

# 注：本流程使用的是经过转化的bedpe
# 单个样本
macs2 callpeak  -g mm -f BEDPE --nomodel --keep-dup all \
  -n SRR11539111 -t ./SRR11539111.Tn5.bedpe \
  --outdir /mnt/d/ATAC/macs2_peaks/

# 循环
cp /mnt/d/ATAC/rmdup/config.raw /mnt/d/ATAC/Tn5_shift/config.raw
cat config.raw | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  macs2 callpeak  -g mm -f BEDPE --nomodel --keep-dup all \
   --cutoff-analysis -n ${sample} -t ./${sample}.Tn5.bedpe \
  --outdir ../macs2_peaks/
done

# 如果用的不是专门双端测序的bedpe，而是bed文件，采用下面代码
# 单个样本
mkdir -p /mnt/d/ATAC/macs2_peaks2/
cd /mnt/d/ATAC/Tn5_shift/
macs2 callpeak  -g mm --nomodel \
  --shift -100 --extsize 200 -n SRR11539111 -t ./SRR11539111.Tn5.bed \
  --outdir /mnt/d/ATAC/macs2_peaks2/

# 循环
cp /mnt/d/ATAC/rmdup/config.raw /mnt/d/ATAC/Tn5_shift/config.raw
cat config.raw | while read id;
do echo $id 
  arr=($id)
  sample=${arr[0]}

  macs2 callpeak  -g mm --nomodel \
  --shift -100 --extsize 200 -n ${sample} -t ./${sample}.Tn5.bed \
  --outdir /mnt/d/ATAC/macs2_peaks2/ 
done
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

6. 
