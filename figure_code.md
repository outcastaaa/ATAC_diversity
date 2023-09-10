# Figure 1
## F1A 将ATAC-seq、RNA-seq样本聚类

!!! 单样本没有去除chrY染色体，会影响聚类

### RNA-seq 
1. 样本准备
```bash
cd /mnt/d/RNA_brain/human/Deseq2
# 先用htseq-count计算每个样本的基因raw count数，得到merge1.csv
while read -r i
do
  cp ../HTseq/${i}.count ./
done < 1.list
```

2. 初步过滤表达量低的基因
```r
library("RColorBrewer")
setwd("D:/RNA_brain/human/Deseq2")
dataframe <- read.csv("merge1.csv", header=TRUE, row.names = 1)
data <- dataframe[-(1:5),]
data <- data[rowSums(data) > 0,]
head(data)
```
```r
# 导入coltdata文件
coldata <- read.table("coldata.csv", row.names = 1, header = TRUE, sep = "," ) 
countdata <- countdata[row.names(coldata)]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)
dds

# 归一化
rld <- rlog(dds, blind=FALSE)

# 聚类热图

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



## F1B 计算sample之间的correlation

软件：deeptools multiBamSummary，bins是给定bin size在全基因组范围内进行coverage的统计。  
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
--corMethod pearson  #spearman,pearson  \
--plotTitle "Pearson Correlation of Average Scores Per Bin" \
--plotFileFormat pdf --whatToPlot heatmap \
--colorMap RdYlBu --plotNumbers --skipZeros \
-o heatmap_PearsonCorr_bamScores.pdf \
--outFileCorMatrix PearsonCorr_bamScores.tab

```
```r
# 利用计算的correlation值画heatmap图
```
### ATAC-seq