# Figure 1
## F1A 将ATAC-seq、RNA-seq样本聚类

!!! 单样本没有去除chrY染色体，可能会影响聚类

### RNA-seq 
* 下载软件
```bash
cd /mnt/d/biosoft
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
tar -zxvf subread-2.0.6-Linux-x86_64.tar.gz
cd subread-2.0.6-Linux-x86_64/bin
chmod +x featureCounts
/mnt/d/biosoft/subread-2.0.6-Linux-x86_64/bin/featureCounts --help
```
1. 样本准备
```bash
# 为了保证样本间更高的相似性，去除Y染色体上的基因统计
# 使用TSS.bed文件
mkdir -p /mnt/xuruizhi/RNA_brain/human/featurecount
cd /mnt/xuruizhi/RNA_brain/human/featurecount
cp /mnt/xuruizhi/ATAC_brain/human/TSS/hg38.TSS.bed  ./

cat hg38.TSS.bed | grep -v "chrY"  >  hg38_rmchrY.bed
wc -l *.bed
#   274031 hg38.TSS.bed
#   273010 hg38_rmchrY.bed

# gtf文件
cd /mnt/xuruizhi/RNA_brain/human/annotation
cat hg38_rmchrY.gtf | tsv-filter --str-eq 3:transcript  > hg38_rmchrY_transcript.gtf
```
2. feature count计算每个样本的基因raw count数，太麻烦，使用bedtools coverage试试看，最终用htseq-count
```bash
# 先用feature count计算每个样本的基因raw count数，得到merge.csv
cd /mnt/xuruizhi/RNA_brain/human/sort_bam
parallel -j 2 "
    featureCounts -T 4 
    -a ../featurecount/hg38_rmchrY.bed 
    -F bed -o ../featurecount/{1}.count {1}.sort.bam 
    2>../featurecount/{1}.log

    featureCounts -T 10 -a $gtf -o read.count -p -B -C -f -t exon -g gene_id  *.bam
-O  Assign reads to all their overlapping meta-features #RNA-seq不用，ATAC-seq可以
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

# bedtools coverage
mkdir -p /mnt/xuruizhi/RNA_brain/human/coverage
cd /mnt/xuruizhi/RNA_brain/human/featurecount
sort -k1,1 -k2,2n hg38_rmchrY.bed > hg38_rmchrY_sort.bed
cp hg38_rmchrY_sort.bed ../coverage
cd ../sort_bam
parallel -j 6 "
   bedtools coverage -sorted -a ../coverage/hg38_rmchrY_sort.bed -b {1}.sort.bam  -counts  -f 1 > ../coverage/{1}.coverage.txt
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')


# htseq-count，统计exon
mkdir -p /mnt/xuruizhi/RNA_brain/human/HTseq_rmchrY_exon
cd /mnt/xuruizhi/RNA_brain/human/sort_bam
parallel -j 6 "
    htseq-count -t exon -s no -r pos -f bam {1}.sort.bam ../annotation/hg38_rmchrY.gtf > ../HTseq_rmchrY_exon/{1}.count  2>../HTseq_rmchrY_exon/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

mkdir -p /mnt/d/RNA_brain/human/HTseq_rmchrY_exon
```

2. TPM归一化
```bash
# 统计非冗余外显子(EXON)长度之和
BiocManager::install('GenomicFeatures')
library(GenomicFeatures)

setwd("D:/RNA_brain/human/HTseq_rmchrY_exon")
txdb <- makeTxDbFromGFF("../annotation/hg38_rmchrY.gtf",format="gtf")

exons_gene <- exonsBy(txdb, by = "gene")
gene_len <- list()
for(i in names(exons_gene)){
    range_info = reduce(exons_gene[[i]])
    width_info = width(range_info)
    sum_len    = sum(width_info)
    gene_len[[i]] = sum_len
}
# gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

class(gene_len)
length(gene_len)

gene_lengths <- t(as.data.frame(gene_len))
write.table(gene_lengths , file = "hg38_genelen.tsv", row.names = TRUE, sep="\t", quote = FALSE, col.names = FALSE)

# =========== 计算TPM ============ 
write.table(count, "TPM_normalize.count", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)
```

3. 计算dist
```r
setwd("D:/RNA_brain/human/TPM")
data <- read.csv("./terminal.csv", header=TRUE, row.names = 1)
log_data <- log10(data+1)
write.csv(log_data, "TPM_normalize.csv")
colnames(log_data) <- c("nonPSM","PSM","nonVLPFC","VLPFC","nonPSM","PSM","nonCRBLM","CRBLM","nonVLPFC","VLPFC","nonVLPFC","VLPFC","nonVLPFC","VLPFC","nonCRBLM","CRBLM","nonCRBLM","CRBLM","nonPSM","PSM","nonVLPFC","VLPFC","nonPSM","PSM","nonVLPFC","VLPFC","nonCRBLM","CRBLM")

sampleDists <- dist(t(log_data))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(8, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```






4. rld归一化做dist

* coldata
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/Deseq2_rmchrY_exon
mkdir -p /mnt/d/RNA_brain/human/Deseq2_rmchrY_exon
cd /mnt/d/RNA_brain/human/Deseq2_rmchrY_exon

vim coldata.csv
"ids","state","condition","treatment"
"SRR21161731","WT","neuron","PSM"
"SRR21161739","WT","neuron","PSM"
"SRR21161882","WT","neuron","PSM"
"SRR21161915","WT","neuron","PSM"
"SRR21161730","WT","non-neuron","nonPSM"
"SRR21161738","WT","non-neuron","nonPSM"
"SRR21161881","WT","non-neuron","nonPSM"
"SRR21161914","WT","non-neuron","nonPSM"
"SRR21161735","WT","neuron","VLPFC"
"SRR21161751","WT","neuron","VLPFC"
"SRR21161760","WT","neuron","VLPFC"
"SRR21161766","WT","neuron","VLPFC"
"SRR21161910","WT","neuron","VLPFC"
"SRR21161932","WT","neuron","VLPFC"
"SRR21161734","WT","non-neuron","nonVLPFC"
"SRR21161750","WT","non-neuron","nonVLPFC"
"SRR21161759","WT","non-neuron","nonVLPFC"
"SRR21161765","WT","non-neuron","nonVLPFC"
"SRR21161909","WT","non-neuron","nonVLPFC"
"SRR21161931","WT","non-neuron","nonVLPFC"
"SRR21161743","WT","neuron","CRBLM"
"SRR21161768","WT","neuron","CRBLM"
"SRR21161781","WT","neuron","CRBLM"
"SRR21161962","WT","neuron","CRBLM"
"SRR21161742","WT","non-neuron","nonCRBLM"
"SRR21161767","WT","non-neuron","nonCRBLM"
"SRR21161780","WT","non-neuron","nonCRBLM"
"SRR21161961","WT","non-neuron","nonCRBLM"

library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

setwd("D:/RNA_brain/human/Deseq2_rmchrY_exon")
dataframe <- read.csv("../HTseq_rmchrY_exon/merge.csv", header=TRUE, row.names = 1)
countdata <- dataframe[-(1:5),]
countdata <- countdata[rowSums(countdata) > 0,]
head(countdata)

coldata <- read.table("coldata.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds



# PCA分析 
# 归一化
rld <- rlog(dds, blind=FALSE)

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















































### ATAC-seq

1. 样本准备
```bash
# 为了保证样本间更高的相似性，去除Y染色体上的基因统计
# gtf文件
cd /mnt/xuruizhi/RNA_brain/human/annotation
cat hg38_rmchrY.gtf | tsv-filter --str-eq 3:transcript  > hg38_rmchrY_transcript.gtf
mkdir -p /mnt/xuruizhi/ATAC_brain/human/annotation
cp ./hg38_rmchrY_transcript.gtf /mnt/xuruizhi/ATAC_brain/human/annotation/
# 根据genebody图，做up+3kb,down+1kb,统计reads情况
awk 'BEGIN{OFS="\t"}{$4=$4-3000; $5=$5+1000; print}' hg38_rmchrY_transcript.gtf > hg38_rmchrY_transcript_adj.gtf
```
2. htseq-count
```bash
# htseq-count
mkdir -p /mnt/xuruizhi/ATAC_brain/human/HTseq_rmchrY_transcript
cd /mnt/xuruizhi/ATAC_brain/human/final
cat 1_all.list | while read id 
do
    htseq-count -n 4 -f bam -r pos --minaqual 30 -t transcript -i gene_id  -s no ${id}.final.bam ../annotation/hg38_rmchrY_transcript.gtf > ../HTseq_rmchrY_transcript/${id}.count  2>../HTseq_rmchrY_transcript/${id}.log
done
# --nonunique all 先不考虑两个都比对上的情况
```



3. rld归一化做dist —— 结果很不好

* coldata
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/Deseq2_rmchrY_transcript
mkdir -p /mnt/d/ATAC_brain/human/Deseq2_rmchrY_transcript
cd /mnt/d/ATAC_brain/human/Deseq2_rmchrY_transcript

vim coldata.csv
"ids","state","condition","treatment"
"SRR21163181","WT","neuron","PSM"
"SRR21163187","WT","neuron","PSM"
"SRR21163294","WT","neuron","PSM"
"SRR21163180","WT","non-neuron","nonPSM"
"SRR21163186","WT","non-neuron","nonPSM"
"SRR21163293","WT","non-neuron","nonPSM"
"SRR21163185","WT","neuron","VLPFC"
"SRR21163208","WT","neuron","VLPFC"
"SRR21163321","WT","neuron","VLPFC"
"SRR21163338","WT","neuron","VLPFC"
"SRR21163184","WT","non-neuron","nonVLPFC"
"SRR21163207","WT","non-neuron","nonVLPFC"
"SRR21163320","WT","non-neuron","nonVLPFC"
"SRR21163337","WT","non-neuron","nonVLPFC"
"SRR21163191","WT","neuron","CRBLM"
"SRR21163217","WT","neuron","CRBLM"
"SRR21163366","WT","neuron","CRBLM"
"SRR21163190","WT","non-neuron","nonCRBLM"
"SRR21163216","WT","non-neuron","nonCRBLM"
"SRR21163365","WT","non-neuron","nonCRBLM"
```
```r
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)

setwd("D:/ATAC_brain/human/Deseq2_rmchrY_transcript")
dataframe <- read.csv("../HTseq_rmchrY_transcript/merge.csv", header=TRUE, row.names = 1)
countdata <- dataframe[-(1:5),]
countdata <- countdata[rowSums(countdata) > 0,]
head(countdata) #39473

coldata <- read.table("coldata.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds



# PCA分析 
# 归一化
rld <- rlog(dds, blind=FALSE)

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





















## F1B 计算sample之间的correlation

软件：deeptools multiBamSummary，bins是给定bin size在全基因组范围内进行coverage的统计。先计算总体情况。  
功能：计算两个以上（含两个）BAM文件的基因组区域的覆盖度。    
小猪文章bin size定为2kb。  
这种相关性也可以通过基础R算法实现。

1. 统计reads在全基因组范围的情况
```bash
# 统计reads在全基因组范围的情况，如果结果不好，统计在基因固定位置的情况
# multiBamSummary bins --bamfiles file1.bam file2.bam -o results.npz -bs <bin size> --outRawCounts readCounts.tab
multiBamSummary bins \
--bamfiles  *.bam
--labels name \
-out counts_per_bin.npz --binSize 2000 \
--minMappingQuality 10 -p 6 --outRawCounts counts_per_bin.tab

# 散点图
plotCorrelation -in treat_results.npz -o treat_results.png --corMethod spearman -p scatterplot
# or
plotCorrelation -in number.of.bins.npz 
--corMethod pearson 
--skipZeros 
--plotTitle "Pearson Correlation of Average Scores Per Transcript" 
--plotFileFormat pdf 
--whatToPlot scatterplot 
-o scatterplot_PearsonCorr_bigwigScores.pdf 
--outFileCorMatrix PearsonCorr_bigwigScores.tab

# 热图
plotCorrelation -in treat_results.npz -o treat_results_heatmap.png --corMethod spearman -p heatmap
# or
plotCorrelation -in number.of.bins.npz 
--corMethod pearson 
--skipZeros 
--plotTitle "Pearson Correlation of Average Scores Per Transcript" 
--plotFileFormat pdf 
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers 
-o heatmap_PearsonCorr_bigwigScores.pdf 
--outFileCorMatrix PearsonCorr_bigwigScores.tab

# 主成分分析
plotPCA -in treat_results.npz  -o pca.png
```
2. 统计reads在功能区域内情况  





### RNA-seq 
```bash
mkdir -p /mnt/xuruizhi/RNA_brain/human/Correlation
cd /mnt/xuruizhi/RNA_brain/human/sort_bam
conda activate py3.8
multiBamSummary bins \
--bamfiles  SRR21161731.sort.bam SRR21161739.sort.bam SRR21161882.sort.bam SRR21161915.sort.bam SRR21161735.sort.bam SRR21161751.sort.bam SRR21161760.sort.bam SRR21161766.sort.bam SRR21161910.sort.bam SRR21161932.sort.bam SRR21161743.sort.bam SRR21161768.sort.bam SRR21161781.sort.bam SRR21161962.sort.bam SRR21161730.sort.bam SRR21161738.sort.bam SRR21161881.sort.bam SRR21161914.sort.bam SRR21161734.sort.bam SRR21161750.sort.bam SRR21161759.sort.bam SRR21161765.sort.bam SRR21161909.sort.bam SRR21161931.sort.bam SRR21161742.sort.bam SRR21161767.sort.bam SRR21161780.sort.bam SRR21161961.sort.bam \
--labels neu_PMC1 neu_PMC2 neu_PMC3 neu_PMC4 neu_VLPFC1 neu_VLPFC2 neu_VLPFC3 neu_VLPFC4 neu_VLPFC5 neu_VLPFC6 neu_CRBLM1 neu_CRBLM2 neu_CRBLM3 neu_CRBLM4 non_PMC1 non_PMC2 non_PMC3 non_PMC4 non_VLPFC1 non_VLPFC2 non_VLPFC3 non_VLPFC4 non_VLPFC5 non_VLPFC6 non_CRBLM1 non_CRBLM2 non_CRBLM3 non_CRBLM4 \
-out counts_per_bin.npz --binSize 2000 --verbose \
--minMappingQuality 10 -p 6 --outRawCounts counts_per_bin.tab

# 是否需要提前去除chrY？性别是否影响？
# 先不去除试试

# 热图
plotCorrelation -in counts_per_bin.npz \
--corMethod pearson \
--plotTitle "Pearson Correlation of Average Scores" \
--plotFileFormat pdf --whatToPlot heatmap \
--colorMap RdYlBu --skipZeros \
-o ../Correlation/heatmap_PearsonCorr_bamScores.pdf \
--outFileCorMatrix ../Correlation/PearsonCorr_bamScores.tab
#--plotNumbers



# 去除质量30以下的
multiBamSummary bins \
--bamfiles  SRR21161731.sort.bam SRR21161739.sort.bam SRR21161882.sort.bam SRR21161915.sort.bam SRR21161735.sort.bam SRR21161751.sort.bam SRR21161760.sort.bam SRR21161766.sort.bam SRR21161910.sort.bam SRR21161932.sort.bam SRR21161743.sort.bam SRR21161768.sort.bam SRR21161781.sort.bam SRR21161962.sort.bam SRR21161730.sort.bam SRR21161738.sort.bam SRR21161881.sort.bam SRR21161914.sort.bam SRR21161734.sort.bam SRR21161750.sort.bam SRR21161759.sort.bam SRR21161765.sort.bam SRR21161909.sort.bam SRR21161931.sort.bam SRR21161742.sort.bam SRR21161767.sort.bam SRR21161780.sort.bam SRR21161961.sort.bam \
--labels neu_PMC1 neu_PMC2 neu_PMC3 neu_PMC4 neu_VLPFC1 neu_VLPFC2 neu_VLPFC3 neu_VLPFC4 neu_VLPFC5 neu_VLPFC6 neu_CRBLM1 neu_CRBLM2 neu_CRBLM3 neu_CRBLM4 non_PMC1 non_PMC2 non_PMC3 non_PMC4 non_VLPFC1 non_VLPFC2 non_VLPFC3 non_VLPFC4 non_VLPFC5 non_VLPFC6 non_CRBLM1 non_CRBLM2 non_CRBLM3 non_CRBLM4 \
-out ../Correlation/counts_per_bin30.npz --binSize 2000 --verbose \
--minMappingQuality 30 -p 6 --outRawCounts ../Correlation/counts_per_bin30.tab

cd ../Correlation
plotCorrelation -in counts_per_bin30.npz \
--corMethod pearson \
--plotTitle "Pearson Correlation of Average Scores" \
--plotFileFormat pdf --whatToPlot heatmap \
--colorMap RdYlBu --skipZeros \
-o heatmap_PearsonCorr_bamScores30.pdf \
--outFileCorMatrix PearsonCorr_bamScores30.tab
```
```r
# 利用计算的correlation值画heatmap图
```
### ATAC-seq
* 质量值10或者30影响不大
```bash
mkdir -p /mnt/xuruizhi/ATAC_brain/human/Correlation
cd /mnt/xuruizhi/ATAC_brain/human/final
conda activate py3.8
multiBamSummary bins \
--bamfiles SRR21163181.final.bam SRR21163187.final.bam SRR21163294.final.bam SRR21163185.final.bam SRR21163208.final.bam SRR21163321.final.bam SRR21163338.final.bam SRR21163191.final.bam SRR21163217.final.bam SRR21163366.final.bam SRR21163180.final.bam SRR21163186.final.bam SRR21163293.final.bam SRR21163184.final.bam SRR21163207.final.bam SRR21163320.final.bam SRR21163337.final.bam SRR21163190.final.bam SRR21163216.final.bam SRR21163365.final.bam \
--labels neu_PMC1 neu_PMC2 neu_PMC3 neu_VLPFC1 neu_VLPFC2 neu_VLPFC3 neu_VLPFC4  neu_CRBLM1 neu_CRBLM2 neu_CRBLM3 non_PMC1 non_PMC2 non_PMC3 non_VLPFC1 non_VLPFC2 non_VLPFC3 non_VLPFC4 non_CRBLM1 non_CRBLM2 non_CRBLM3 \
-out ../Correlation/counts_per_bin.npz --binSize 2000 --verbose \
--minMappingQuality 10 -p 6 --outRawCounts ../Correlation/counts_per_bin.tab

# 是否需要提前去除chrY？性别是否影响？
# 先不去除试试

# 热图
cd ../Correlation
plotCorrelation -in counts_per_bin.npz \
--corMethod pearson \
--plotTitle "Pearson Correlation of Average Scores" \
--plotFileFormat pdf --whatToPlot heatmap \
--colorMap RdYlBu --skipZeros \
-o heatmap_PearsonCorr_bamScores.pdf \
--outFileCorMatrix PearsonCorr_bamScores.tab
#--plotNumbers


# 去除质量30以下的
multiBamSummary bins \
--bamfiles SRR21163181.final.bam SRR21163187.final.bam SRR21163294.final.bam SRR21163185.final.bam SRR21163208.final.bam SRR21163321.final.bam SRR21163338.final.bam SRR21163191.final.bam SRR21163217.final.bam SRR21163366.final.bam SRR21163180.final.bam SRR21163186.final.bam SRR21163293.final.bam SRR21163184.final.bam SRR21163207.final.bam SRR21163320.final.bam SRR21163337.final.bam SRR21163190.final.bam SRR21163216.final.bam SRR21163365.final.bam \
--labels neu_PMC1 neu_PMC2 neu_PMC3 neu_VLPFC1 neu_VLPFC2 neu_VLPFC3 neu_VLPFC4  neu_CRBLM1 neu_CRBLM2 neu_CRBLM3 non_PMC1 non_PMC2 non_PMC3 non_VLPFC1 non_VLPFC2 non_VLPFC3 non_VLPFC4 non_CRBLM1 non_CRBLM2 non_CRBLM3 \
-out ../Correlation/counts_per_bin30.npz --binSize 2000 --verbose \
--minMappingQuality 30 -p 6 --outRawCounts ../Correlation/counts_per_bin30.tab

cd ../Correlation
plotCorrelation -in counts_per_bin30.npz \
--corMethod pearson \
--plotTitle "Pearson Correlation of Average Scores" \
--plotFileFormat pdf --whatToPlot heatmap \
--colorMap RdYlBu --skipZeros \
-o heatmap_PearsonCorr_bamScores30.pdf \
--outFileCorMatrix PearsonCorr_bamScores30.tab
```
```r
# 利用计算的correlation值画heatmap图
```


## 统计RNA-seq基因表达分布情况
```bash
ls *.count | while read id
do
cat ${id} | cut -f 2 | head -n -5 > ${id%%.*}.txt
done
```
```r
setwd("D:/RNA_brain/human/HTseq_rmchrY_exon")
# 画peak length 直方图
summary_list <- list()
file_list <- list.files(path = "./", pattern = "\\d+\\.txt", full.names = TRUE)

for (file in file_list) {
  data <- read.table(file, encoding = "UTF-8")
  dim_data <- dim(data)
  summary_list[[file]] <- dim_data

  # 生成直方图并保存为png文件
  png(paste0(file, "_hist.png"))
  # 设置横坐标轴范围，例如从0到最大数据值
#   xlim_range <- c(0, quantile(data[,1], 0.99))
  xlim_range <- c(0, 5000)
  # 设置合适的 breaks 值，根据数据范围和希望的箱子数量进行调整
  breaks_value <- 1  # 根据需要调整
  hist(abs(as.numeric(data[,1])), breaks = breaks_value, 
       xlab = "count", ylab = "Frequency", main = file, xlim = xlim_range)
  dev.off()
}
```