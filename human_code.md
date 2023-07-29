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
# human 70个
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
SRR21163214
SRR21163215
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


mkdir -p ./sra
cd ./sra
cp /mnt/d/perl/perl_scripts/download_srr.pl ./
cp /mnt/xuruizhi/ATAC_brain/human/HUMAN.list ./

# 批量下载
cat HUMAN.list | parallel -k -j 6 "
  echo {} >> ./download.log
  perl download_srr.pl --output-dir . --srr {} >> ./download.log 2>&1
"
cat download.log | grep " downloaded successfully"
```
