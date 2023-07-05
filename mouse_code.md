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
SRR14614715
SRR3595211
SRR3595212
SRR3595213
SRR3595214
SRR13443549
SRR13443553
SRR13443554
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


mkdir -p /mnt/xuruizhi/ATAC_brain/mouse/sequence
# 在超算中转换格式
rsync -av /mnt/xuruizhi/ATAC_brain/ \
wangq@202.119.37.251:/share/home/wangq/xuruizhi/brain/


fqdir= ~/xuruizhi/brain/fastq/物种
ls *.sra | while read id
do
 echo "fastq-dump --gzip --split-3 -O ${fqdir} ${id}"
done >sra2fq.sh
# 提交后台运行命令，脚本文件后缀为.sh，日志文件后缀为.log，运行脚本的命令为sh
sh sra2fq.sh > sra2fq.log 




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
bsub -q mpi -n 24 -J sra2fq_mouse -o ~/xuruizhi/brain/brain/fastq/mouse " \
cat MOUSE.list |  parallel --no-run-if-empty --linebuffer -k -j 10 \
'/share/home/wangq/bin/fastq-dump.3.0.0 --gzip --split-3 \
-O ~/xuruizhi/brain/brain/fastq/mouse ~/xuruizhi/brain/brain/sequence/{}'"

```
2. 基因组数据

# 3. 比对前质控


























```bash
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


# 写单端批量
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