setwd("D:/ATAC_brain/peaks")
getwd()
file_list <- list.files(pattern = "SRR\\d+_peaks_length\\.txt")

summary_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE)
  summary_list[[file]] <- summary(data)
}

for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}

summary_list <- list()
file_list <- list.files(path = "./", pattern = "SRR\\d+_peaks_length\\.txt", full.names = TRUE)

for (file in file_list) {
  data <- read.table(file)
  dim_data <- dim(data)
  summary_list[[file]] <- dim_data
  
  # 生成直方图并保存为PDF文件
  pdf(paste0(file, "_hist.pdf"))
  hist(abs(as.numeric(data[, 1])), breaks = 500, xlab = "Fragment length (bp)", ylab = "Frequency", main = file)
  dev.off()
}

for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}