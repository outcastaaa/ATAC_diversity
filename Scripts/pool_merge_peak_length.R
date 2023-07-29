setwd("D:/ATAC_brain/mouse/")
file_list <- list.files(path = "./IDR", pattern = ".*_pool_merge_length\\.txt$", full.names = TRUE)
summary_list <- list()
for (file in file_list) {
  data <- read.table(file, header = TRUE)
  summary_list[[file]] <- summary(data)
}

for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}

# 画直方图
summary_list <- list()
file_list <- list.files(path = "./IDR", pattern = ".*_pool_merge_length\\.txt$", full.names = TRUE)
for (file in file_list) {
  data <- read.table(file)
  dim_data <- dim(data)
  summary_list[[file]] <- dim_data
  
  png(paste0(file, "_hist.png"))
  hist(abs(as.numeric(data[, 1])), breaks = 500, xlab = "Fragment length (bp)", ylab = "Frequency", main = file)
  dev.off()
}
for (file in file_list) {
  cat("Summary for", file, ":\n")
  print(summary_list[[file]])
}