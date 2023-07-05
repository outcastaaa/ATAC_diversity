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
SRR11179779.    0e6ac756b030cc84f42de5bce499f5fd
SRR11179780.    d1766b149fe0182ae872bf79a3decd5a
SRR11179781.    b93250c0ce5affeb603fa071c78fcce0
SRR13049359.    3f280eaf565605981db8d3062bf6a6a3
SRR13049360.    70677096eb974366da7152f36fa0c5b0
SRR13049361.    119bfddf4301e49835661bf24693199c
SRR13049362.    674996a273c74bbb3d52cfc025e155da
SRR13049363.    a54bb8dad247d4e7c114f12fd1b8791a
SRR13049364.    bf1dc5b5a62011235fdf854f860c98b8
SRR14362271.    8667509648f566de1e8b3b3e3638b3e0
SRR14362272.    bb509a5db71ff478426c773aa2c9770d
SRR14362275.    5eb17b41b9ac7c2caa4967752c416c98
SRR14362276.    075c3926a6a49c15ce85ac9194cc5370
SRR14362281.    c2cee7a6f408c5db7c8a81119e2ee43b
SRR14362282.    33a21b0f50030658b5a2d48f6573c1c3
SRR14614715×
SRR3595211.    fd2862221c839f1da0bf667704eb63c1
SRR3595212.    21076454a9730edede42d98f09aa10bc
SRR3595213.    c277a3891a38d31a0beb94105937e45d
SRR3595214.    b2ca2eee6d312f08831d9e27657dda37
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


cd /scratch/wangq/xrz/ATAC_brain/mouse/sequence



```
2. 基因组数据

# 3. 比对前质控
