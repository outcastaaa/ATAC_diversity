#!/usr/bin R
# 人类GO分析代码
library(biomaRt)
library(ChIPseeker)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

region <- c("xxx")

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

  mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
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